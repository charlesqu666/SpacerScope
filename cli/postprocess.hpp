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
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

using CombinedResult = SearchHit;

struct AlignmentEvent {
    std::string op;
    int pos = 0;
    char query_base = '-';
    char target_base = '-';
};

struct ScoredResult {
    CombinedResult base;
    double score = 0.0;
    std::vector<AlignmentEvent> details;
};

struct RegionHeader {
    std::string chrom;
    uint32_t start_1based = 0;
    uint32_t end_1based = 0;
    bool is_reverse = false;
};

bool parse_region_header(const std::string& header, RegionHeader& region);

std::vector<std::vector<AlignmentEvent>> align_left_fixed_cpp(const std::string& target,
                                                              const std::string& query,
                                                              int max_indel,
                                                              int max_mismatch,
                                                              int pam_direction);

std::vector<AlignmentEvent> filter_alignment_events(const std::vector<AlignmentEvent>& events);

double predict_off_target_cpp(const std::vector<AlignmentEvent>& mismatch_positions);

std::string alignment_events_to_json(const std::vector<AlignmentEvent>& events);

std::vector<ScoredResult> score_results_cpp(const std::vector<CombinedResult>& combined_results,
                                            const GenomeReference& target_genome,
                                            int max_indel,
                                            int max_mismatch,
                                            int pam_direction,
                                            int num_threads,
                                            size_t& discarded_count);

void write_scored_results_tsv(const std::string& output_path,
                              const std::vector<ScoredResult>& results);

std::string format_score_cpp(double score);
