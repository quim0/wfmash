/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap 
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP 
#define COMPUTE_ALIGNMENTS_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>
#include <thread>
#include <memory>
#include <omp.h>
#include <htslib/faidx.h>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/wflign/src/wflign.hpp"
#include "common/atomic_queue/atomic_queue.h"
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "common/utils.hpp"

namespace align
{

long double float2phred(long double prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

struct seq_record_t {
    MappingBoundaryRow currentRecord;
    std::string mappingRecordLine;
    std::string refSequence;  
    std::string querySequence;
    uint64_t refStartPos;
    uint64_t refLen;
    uint64_t refTotalLength;
    uint64_t queryStartPos;
    uint64_t queryLen;
    uint64_t queryTotalLength;

    seq_record_t(const MappingBoundaryRow& c, const std::string& r, 
                 const std::string& ref, uint64_t refStart, uint64_t refLength, uint64_t refTotalLength,
                 const std::string& query, uint64_t queryStart, uint64_t queryLength, uint64_t queryTotalLength)
        : currentRecord(c)
        , mappingRecordLine(r)
        , refSequence(ref)
        , querySequence(query)
        , refStartPos(refStart)
        , refLen(refLength)
        , refTotalLength(refTotalLength)
        , queryStartPos(queryStart)
        , queryLen(queryLength)
        , queryTotalLength(queryTotalLength)
        { }
};

/**
 * @brief A single-producer, multi-consumer (SPMC) atomic queue for storing pointers to seq_record_t objects.
 *
 * This queue is designed for a setup where there is a single producer and multiple consumers.
 * The producer enqueues pointers to seq_record_t objects, which represent sequences to be processed.
 * Multiple consumers dequeue these pointers and perform long-running alignment processes on the sequences.
 *
 * The queue has the following characteristics:
 * - Capacity: 1024 elements
 * - Default value for empty elements: nullptr
 * - MINIMIZE_CONTENTION: true (minimizes contention among consumers)
 * - MAXIMIZE_THROUGHPUT: true (optimized for high throughput)
 * - TOTAL_ORDER: false (relaxed memory ordering for better performance)
 * - SPSC: false (single-producer, multi-consumer mode)
 */
typedef atomic_queue::AtomicQueue<seq_record_t*, 1024, nullptr, true, true, false, false> seq_atomic_queue_t;

/**
 * @brief A multi-producer, single-consumer (MPSC) atomic queue for storing pointers to std::string objects.
 *
 * This queue is designed for a setup where there are multiple producers and a single consumer.
 * Multiple producers enqueue pointers to std::string objects, which represent PAF (Pairwise Alignment Format) strings.
 * The single consumer dequeues these pointers and writes out the PAF strings.
 *
 * The queue has the following characteristics:
 * - Capacity: 1024 elements
 * - Default value for empty elements: nullptr
 * - MINIMIZE_CONTENTION: true (minimizes contention among producers)
 * - MAXIMIZE_THROUGHPUT: true (optimized for high throughput)
 * - TOTAL_ORDER: false (relaxed memory ordering for better performance)
 * - SPSC: false (multi-producer, single-consumer mode)
 */
typedef atomic_queue::AtomicQueue<std::string*, 1024, nullptr, true, true, false, false> paf_atomic_queue_t;


  /**
   * @class     align::Aligner
   * @brief     compute alignments and generate sam output
   *            from mashmap mappings
   */
  class Aligner
  {
    private:

      //algorithm parameters
      const align::Parameters &param;

      faidx_t* ref_faidx;
      faidx_t* query_faidx;

    public:

      explicit Aligner(const align::Parameters &p) : param(p) {
          assert(param.refSequences.size() == 1);
          assert(param.querySequences.size() == 1);
          ref_faidx = fai_load(param.refSequences.front().c_str());
          query_faidx = fai_load(param.querySequences.front().c_str());
      }

      ~Aligner() {
          fai_destroy(ref_faidx);
          fai_destroy(query_faidx);
      }

      /**
       * @brief                 compute alignments
       */

      void compute()
      {
        this->computeAlignments();
      }

      /**
       * @brief       parse mashmap row sequence
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline static void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord, const uint64_t target_padding) {
          std::stringstream ss(mappingRecordLine); // Insert the string into a stream
          std::string word; // Have a buffer string

          vector<std::string> tokens; // Create vector to hold our words

          while (ss >> word) {
              tokens.push_back(word);
          }

          // Check if the number of tokens is at least 13
          if (tokens.size() < 13) {
              throw std::runtime_error("[wfmash::align::parseMashmapRow] Error! Invalid mashmap mapping record: " + mappingRecordLine);
          }

          // Extract the mashmap identity from the string
          const vector<string> mm_id_vec = skch::CommonFunc::split(tokens[12], ':');
          // if the estimated identity is missing, avoid assuming too low values
          const float mm_id = wfmash::is_a_number(mm_id_vec.back()) ? std::stof(mm_id_vec.back()) : skch::fixed::percentage_identity;

          // Parse chain info if present (expecting format "chain:i:id.pos.len" in tokens[14])
          int32_t chain_id = -1;
          int32_t chain_length = 1;
          int32_t chain_pos = 1;
          if (tokens.size() > 14) {
              const vector<string> chain_vec = skch::CommonFunc::split(tokens[14], ':');
              if (chain_vec.size() == 3 && chain_vec[0] == "chain" && chain_vec[1] == "i") {
                  // Split the id.pos.len format
                  const vector<string> chain_parts = skch::CommonFunc::split(chain_vec[2], '.');
                  if (chain_parts.size() == 3) {
                      chain_id = std::stoi(chain_parts[0]);
                      chain_pos = std::stoi(chain_parts[1]);
                      chain_length = std::stoi(chain_parts[2]);
                  }
              }
          }

          // Save words into currentRecord
          {
              currentRecord.qId = tokens[0];
              currentRecord.qStartPos = std::stoi(tokens[2]);
              currentRecord.qEndPos = std::stoi(tokens[3]);
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = tokens[5];
              const uint64_t ref_len = std::stoi(tokens[6]);
              currentRecord.chain_id = chain_id;
              currentRecord.chain_length = chain_length;
              currentRecord.chain_pos = chain_pos;

              // Apply target padding while ensuring we don't go below 0 or above reference length
              uint64_t rStartPos = std::stoi(tokens[7]);
              uint64_t rEndPos = std::stoi(tokens[8]);

              // Always apply target padding
              if (target_padding > 0) {
                  if (rStartPos >= target_padding) {
                      rStartPos -= target_padding;
                  } else {
                      rStartPos = 0;
                  }
                  if (rEndPos + target_padding <= ref_len) {
                      rEndPos += target_padding;
                  } else {
                      rEndPos = ref_len;
                  }
              }

              // Validate coordinates against reference length
              if (rStartPos >= ref_len || rEndPos > ref_len) {
                  std::cerr << "[parse-debug] ERROR: Coordinates exceed reference length!" << std::endl;
                  throw std::runtime_error("[wfmash::align::parseMashmapRow] Error! Coordinates exceed reference length: " 
                                         + std::to_string(rStartPos) + "-" + std::to_string(rEndPos) 
                                         + " (ref_len=" + std::to_string(ref_len) + ")");
              }

              currentRecord.rStartPos = rStartPos;
              currentRecord.rEndPos = rEndPos;
              currentRecord.mashmap_estimated_identity = mm_id;
          }
      }

  private:

seq_record_t* createSeqRecord(const MappingBoundaryRow& currentRecord, 
                              const std::string& mappingRecordLine,
                              faidx_t* ref_faidx,
                              faidx_t* query_faidx) {
    // Get the reference sequence length
    const int64_t ref_size = faidx_seq_len(ref_faidx, currentRecord.refId.c_str());
    // Get the query sequence length
    const int64_t query_size = faidx_seq_len(query_faidx, currentRecord.qId.c_str());

    // Compute padding for sequence extraction
    const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor
        ? param.wflign_max_len_minor : currentRecord.rStartPos;
    const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor
        ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

    // Extract reference sequence
    int64_t ref_len;
    char* ref_seq = faidx_fetch_seq64(ref_faidx, currentRecord.refId.c_str(),
                                      currentRecord.rStartPos - head_padding, 
                                      currentRecord.rEndPos - 1 + tail_padding, &ref_len);

    // Extract query sequence
    int64_t query_len;
    char* query_seq = faidx_fetch_seq64(query_faidx, currentRecord.qId.c_str(),
                                        currentRecord.qStartPos, currentRecord.qEndPos - 1, &query_len);

    // Create a new seq_record_t object for the alignment
    seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                         std::string(ref_seq, ref_len), 
                                         currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                         std::string(query_seq, query_len), 
                                         currentRecord.qStartPos, query_len, query_size);

    // Clean up
    free(ref_seq);
    free(query_seq);

    return rec;
}

std::string processAlignment(seq_record_t* rec,
                             wfa::WFAlignerGapAffine2Pieces& wf_aligner) {
    std::string& ref_seq = rec->refSequence;
    std::string& query_seq = rec->querySequence;

    skch::CommonFunc::makeUpperCaseAndValidDNA(ref_seq.data(), ref_seq.length());
    skch::CommonFunc::makeUpperCaseAndValidDNA(query_seq.data(), query_seq.length());

    // Adjust the reference sequence to start from the original start position
    char* ref_seq_ptr = &ref_seq[rec->currentRecord.rStartPos - rec->refStartPos];

    std::vector<char> queryRegionStrand(query_seq.size() + 1);

    if(rec->currentRecord.strand == skch::strnd::FWD) {
        std::copy(query_seq.begin(), query_seq.end(), queryRegionStrand.begin());
    } else {
        skch::CommonFunc::reverseComplement(query_seq.data(), queryRegionStrand.data(), query_seq.size());
    }

    // Set up penalties for biWFA
    wflign_penalties_t wfa_penalties;
    wfa_penalties.match = 0;
    wfa_penalties.mismatch = param.wfa_patching_mismatch_score;
    wfa_penalties.gap_opening1 = param.wfa_patching_gap_opening_score1;
    wfa_penalties.gap_extension1 = param.wfa_patching_gap_extension_score1;
    wfa_penalties.gap_opening2 = param.wfa_patching_gap_opening_score2;
    wfa_penalties.gap_extension2 = param.wfa_patching_gap_extension_score2;

    std::stringstream output;

    // Do direct biWFA alignment
    wflign::wavefront::do_biwfa_alignment_reuse_aligner(
        rec->currentRecord.qId,
        queryRegionStrand.data(),
        rec->queryTotalLength,
        rec->queryStartPos,
        rec->queryLen,
        rec->currentRecord.strand != skch::strnd::FWD,
        rec->currentRecord.refId,
        ref_seq_ptr,
        rec->refTotalLength,
        rec->currentRecord.rStartPos,
        rec->currentRecord.rEndPos - rec->currentRecord.rStartPos,
        output,
        wfa_penalties,
        param.emit_md_tag,
        !param.sam_format,
        param.no_seq_in_sam,
        param.min_identity,
        param.wflign_max_len_minor,
        rec->currentRecord.mashmap_estimated_identity,
        rec->currentRecord.chain_id,
        rec->currentRecord.chain_length,
        rec->currentRecord.chain_pos,
        wf_aligner);

    return output.str();
}

void write_sam_header(std::ofstream& outstream) {
    for(const auto &fileName : param.refSequences) {
        // check if there is a .fai
        std::string fai_name = fileName + ".fai";
        if (fs::exists(fai_name)) {
            // if so, process the .fai to determine our sequence length
            std::string line;
            std::ifstream in(fai_name.c_str());
            while (std::getline(in, line)) {
                auto line_split = skch::CommonFunc::split(line, '\t');
                const std::string seq_name = line_split[0];
                const uint64_t seq_len = std::stoull(line_split[1]);
                outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq_len << "\n";
            }
        } else {
            // if not, warn that this is expensive
            std::cerr << "[wfmash::align] WARNING, no .fai index found for " << fileName << ", reading the file to prepare SAM header (slow)" << std::endl;
            seqiter::for_each_seq_in_file(
                fileName, {}, "",
                [&](const std::string& seq_name,
                    const std::string& seq) {
                    outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq.length() << "\n";
                });
        }
    }
    outstream << "@PG\tID:wfmash\tPN:wfmash\tVN:" << WFMASH_GIT_VERSION << "\tCL:wfmash\n";
}

void writer_thread(const std::string& output_file,
                   paf_atomic_queue_t& paf_queue,
                   std::atomic<bool>& reader_done,
                   std::atomic<bool>& processor_done,
                   const std::vector<std::atomic<bool>>& worker_working) {
    std::ofstream outstream(output_file);
    // if the output file is SAM, we write the header
    if (param.sam_format) {
        write_sam_header(outstream);
    }

    if (!outstream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open output file: " + output_file);
    }

    auto all_workers_done = [&]() {
        return std::all_of(worker_working.begin(), worker_working.end(),
                           [](const std::atomic<bool>& w) { return !w.load(); });
    };

    while (true) {
        std::string* paf_output = nullptr;
        if (paf_queue.try_pop(paf_output)) {
            outstream << *paf_output;
            outstream.flush();
            delete paf_output;
        } else if (reader_done.load() && processor_done.load() && paf_queue.was_empty() && all_workers_done()) {
            break;
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

    outstream.close();
}

void processLine(const std::string& line,
                 faidx_t* local_ref_faidx,
                 faidx_t* local_query_faidx,
                 std::vector<std::string>& paf_lines,
                 progress_meter::ProgressMeter& progress,
                 wfa::WFAlignerGapAffine2Pieces& wf_aligner,
                 const size_t linenum) {
    if (!line.empty()) {
        // (2) Processor thread
        MappingBoundaryRow currentRecord;
        parseMashmapRow(line, currentRecord, param.target_padding);
        seq_record_t* rec = createSeqRecord(currentRecord, line, local_ref_faidx, local_query_faidx);

        // (3) worker_thread

        std::string alignment_output = processAlignment(rec, wf_aligner);

        // Parse the alignment output to find CIGAR string and coordinates
        std::stringstream ss(alignment_output);
        std::string line;
        while (std::getline(ss, line)) {
            if (line.empty()) continue;

            std::vector<std::string> fields;
            std::stringstream field_ss(line);
            std::string field;
            while (field_ss >> field) {
                fields.push_back(field);
            }

            // Find the CIGAR string field (should be after cg:Z:)
            auto cigar_it = std::find_if(fields.begin(), fields.end(),
                [](const std::string& s) { return s.substr(0, 5) == "cg:Z:"; });

            std::string paf_str;
            if (cigar_it != fields.end()) {
                std::string cigar = cigar_it->substr(5); // Remove cg:Z: prefix
                uint64_t ref_start = std::stoull(fields[7]);
                uint64_t ref_end = std::stoull(fields[8]);


                // Just pass through the CIGAR string and coordinates unchanged
                // The trimming is now handled in wflign namespace

                // Reconstruct the line
                std::string new_line;
                for (const auto& f : fields) {
                    if (!new_line.empty()) new_line += '\t';
                    new_line += f;
                }
                new_line += '\n';

                // Push the modified alignment output to the paf_queue
                paf_str = std::string(std::move(new_line));
            } else {
                // If no CIGAR string found, output the line unchanged
                paf_str = std::string(line + '\n');
            }
            paf_lines[linenum] = paf_str;
        }
        uint64_t alignment_length = rec->currentRecord.qEndPos - rec->currentRecord.qStartPos;
        progress.increment(alignment_length);
        delete rec;
    }
}

void computeAlignments() {
    /*
    std::atomic<size_t> total_alignments_queued(0);
    std::atomic<bool> reader_done(false);
    std::atomic<bool> processor_done(false);

    // Create queues
    atomic_queue::AtomicQueue<std::string*, 1024> line_queue;
    seq_atomic_queue_t seq_queue;
    paf_atomic_queue_t paf_queue;  // Add this line
    */

    // Calculate max_processors based on the number of worker threads
    size_t max_processors = std::max(1UL, static_cast<unsigned long>(param.threads));

    // Calculate total alignment length
    // TODO: Avoid reading the whole file twice
    uint64_t total_alignment_length = 0;
    {
        std::ifstream mappingListStream(param.mashmapPafFile);
        if (!mappingListStream.is_open()) {
            throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: " + param.mashmapPafFile);
        }
        std::string mappingRecordLine;
        MappingBoundaryRow currentRecord;

        while(std::getline(mappingListStream, mappingRecordLine)) {
            if (!mappingRecordLine.empty()) {
                parseMashmapRow(mappingRecordLine, currentRecord, param.target_padding);
                total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
            }
        }
    }

    // Create progress meter
    progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align] aligned");

    // (1) Single reader thread
    std::ifstream mappingListStream(param.mashmapPafFile);
    if (!mappingListStream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: " + param.mashmapPafFile);
    }

    std::ofstream outstream(param.pafOutputFile);
    // if the output file is SAM, we write the header
    if (param.sam_format) {
        write_sam_header(outstream);
    }

    if (!outstream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open output file: " + param.pafOutputFile);
    }

    std::string line;
    std::vector<std::string> lines;
    while (std::getline(mappingListStream, line)) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }

    std::vector<std::string> paf_lines;
    paf_lines.resize(lines.size());

    omp_set_num_threads(max_processors);

    #pragma omp parallel
    {
    // Create a single aligner for each thread to be reused
    wflign_penalties_t penalties;
    penalties.match = 0;
    penalties.mismatch = param.wfa_patching_mismatch_score;
    penalties.gap_opening1 = param.wfa_patching_gap_opening_score1;
    penalties.gap_extension1 = param.wfa_patching_gap_extension_score1;
    penalties.gap_opening2 = param.wfa_patching_gap_opening_score2;
    penalties.gap_extension2 = param.wfa_patching_gap_extension_score2;

    wfa::WFAlignerGapAffine2Pieces wf_aligner(
        0,  // match
        penalties.mismatch,
        penalties.gap_opening1,
        penalties.gap_extension1,
        penalties.gap_opening2,
        penalties.gap_extension2,
        wfa::WFAligner::Alignment,
        wfa::WFAligner::MemoryUltralow);
    wf_aligner.setHeuristicNone();

    faidx_t* local_ref_faidx = fai_load(param.refSequences.front().c_str());
    faidx_t* local_query_faidx = fai_load(param.querySequences.front().c_str());

    // Process the lines in parallel
    #pragma omp for
    for (size_t i=0; i<lines.size(); i++) {
        processLine(lines[i], local_ref_faidx, local_query_faidx, paf_lines, progress, wf_aligner, i);
    }
    fai_destroy(local_ref_faidx);
    fai_destroy(local_query_faidx);
    }

    for (const auto& paf_line : paf_lines) {
        if (!paf_line.empty())
            outstream << paf_line;
    }

    mappingListStream.close();
    progress.finish();

} // computeAlignments
}; // Aligner class
} // namespace align


#endif
