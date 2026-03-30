#pragma once

#include "types.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <string_view>
#include <cstdint>

// Simple 2-bit Encoding
typedef uint64_t encoded_seq_t;

struct MismatchPreencodedTargets {
    bool fixed_len_fast_path = false;
    size_t fixed_target_len = 0;
    std::vector<encoded_seq_t> target_codes;
    std::vector<int> target_lengths;
};

// Legacy interface
void process_mismatch_search(const std::string& file_a, 
                             const std::string& file_b, 
                             const std::string& output_filename, 
                             int max_mismatch, 
                             int num_threads);

// Main interface - ALWAYS uses simple 2-bit brute-force (no strategy switching)
void process_mismatch_search_loc(const std::vector<SpacerLocation>& query_locs,
                                 const std::vector<SpacerLocation>& target_locs,
                                 const GenomeReference& query_genome,
                                 const GenomeReference& target_genome,
                                 const std::string& output_filename,
                                 int max_mismatch,
                                 int num_threads);

std::vector<SearchHit> process_mismatch_search_loc_collect(const std::vector<SpacerLocation>& query_locs,
                                                           const std::vector<SpacerLocation>& target_locs,
                                                           const GenomeReference& query_genome,
                                                           const GenomeReference& target_genome,
                                                           int max_mismatch,
                                                           int num_threads);

MismatchPreencodedTargets preencode_targets_for_mismatch(const std::vector<SpacerLocation>& target_locs,
                                                         const GenomeReference& target_genome,
                                                         int num_threads);

std::vector<SearchHit> process_mismatch_search_loc_collect_preencoded(
    const std::vector<SpacerLocation>& query_locs,
    const std::vector<SpacerLocation>& target_locs,
    const MismatchPreencodedTargets& preencoded_targets,
    const GenomeReference& query_genome,
    const GenomeReference& target_genome,
    int max_mismatch,
    int num_threads);

// Dummy for compatibility
typedef std::vector<std::unordered_map<uint64_t, std::vector<uint32_t>>> MismatchIndexOptimized;
MismatchIndexOptimized build_mismatch_index_optimized(const SequenceList& sequences_b, int max_mismatch, int num_threads);
