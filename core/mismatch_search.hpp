/*
 * SpacerScope - gRNA Design and Off-Target Effect Analysis Tool
 * Copyright (C) 2026 charlesqu666
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful for CRISPR
 * research and therapeutic applications, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 * 
 * Source code is available at: https://github.com/charlesqu666/SpacerScope
 */
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
