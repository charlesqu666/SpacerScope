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
#include "postprocess.hpp"

#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string_view>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double kAP = 0.0014;
constexpr double kBP = 0.0021;
constexpr double kCP = 0.3055;

char complement_base(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

std::string reverse_complement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc.push_back(complement_base(*it));
    }
    return rc;
}

double p_pos(int i) {
    return kAP * i * i + kBP * i + kCP;
}

std::string extract_padded_target_sequence(const GenomeReference& genome,
                                           const RegionHeader& region,
                                           int pam_direction,
                                           int pad) {
    const size_t chrom_idx = genome.get_chrom_idx(region.chrom);
    if (chrom_idx == static_cast<size_t>(-1)) return "";
    if (chrom_idx >= genome.chrom_sequences.size()) return "";

    const size_t chrom_len = genome.chrom_sequences[chrom_idx].size();
    size_t start_1based = static_cast<size_t>(region.start_1based);
    size_t end_1based = static_cast<size_t>(region.end_1based);
    if (start_1based == 0 || end_1based == 0 || start_1based > end_1based || end_1based > chrom_len) {
        return "";
    }

    if (pam_direction == -1) {
        if (!region.is_reverse) {
            start_1based = (start_1based > static_cast<size_t>(pad)) ? (start_1based - static_cast<size_t>(pad)) : 1;
        } else {
            end_1based = std::min(chrom_len, end_1based + static_cast<size_t>(pad));
        }
    } else {
        if (!region.is_reverse) {
            end_1based = std::min(chrom_len, end_1based + static_cast<size_t>(pad));
        } else {
            start_1based = (start_1based > static_cast<size_t>(pad)) ? (start_1based - static_cast<size_t>(pad)) : 1;
        }
    }

    std::string_view seq = genome.get_sequence(chrom_idx, start_1based - 1, end_1based - start_1based + 1);
    if (seq.empty()) return "";
    if (region.is_reverse) return reverse_complement(seq);
    return std::string(seq);
}

bool event_has_indel(const std::vector<AlignmentEvent>& events) {
    for (const auto& event : events) {
        if (event.op == "ins" || event.op == "del") return true;
    }
    return false;
}

bool event_tuple_less(const AlignmentEvent& lhs, const AlignmentEvent& rhs) {
    if (lhs.op != rhs.op) return lhs.op < rhs.op;
    if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
    if (lhs.query_base != rhs.query_base) return lhs.query_base < rhs.query_base;
    return lhs.target_base < rhs.target_base;
}

}  // namespace

bool parse_region_header(const std::string& header, RegionHeader& region) {
    const size_t last_colon = header.rfind(':');
    if (last_colon == std::string::npos || last_colon + 1 >= header.size()) return false;
    const size_t second_last_colon = header.rfind(':', last_colon - 1);
    if (second_last_colon == std::string::npos || second_last_colon + 1 >= last_colon) return false;

    region.chrom = header.substr(0, second_last_colon);
    const std::string range = header.substr(second_last_colon + 1, last_colon - second_last_colon - 1);
    const std::string strand = header.substr(last_colon + 1);
    const size_t dash = range.find('-');
    if (dash == std::string::npos) return false;

    region.start_1based = static_cast<uint32_t>(std::stoul(range.substr(0, dash)));
    region.end_1based = static_cast<uint32_t>(std::stoul(range.substr(dash + 1)));
    region.is_reverse = (strand == "-");
    return true;
}

std::vector<std::vector<AlignmentEvent>> align_left_fixed_cpp(const std::string& target,
                                                              const std::string& query,
                                                              int max_indel,
                                                              int max_mismatch,
                                                              int pam_direction) {
    std::vector<std::vector<AlignmentEvent>> results;

    std::string t_seq;
    std::string q_seq;
    if (pam_direction == -1) {
        t_seq.assign(target.rbegin(), target.rend());
        q_seq.assign(query.rbegin(), query.rend());
    } else {
        t_seq = target;
        q_seq = query;
    }

    const int q_len = static_cast<int>(q_seq.size());
    const int t_len = static_cast<int>(t_seq.size());
    std::vector<AlignmentEvent> path;
    path.reserve(static_cast<size_t>(q_len + max_indel));

    auto solve = [&](auto&& self, int q_idx, int t_idx, int mm_cnt, int id_cnt, bool has_indel) -> void {
        if (has_indel) {
            if (mm_cnt + id_cnt > max_indel) return;
        } else {
            if (mm_cnt > max_mismatch) return;
        }
        if (id_cnt > max_indel) return;

        if (q_idx == q_len) {
            std::vector<AlignmentEvent> final_path = path;
            if (pam_direction == -1) {
                std::reverse(final_path.begin(), final_path.end());
            }
            results.push_back(std::move(final_path));
            return;
        }

        if (t_idx == t_len && id_cnt >= max_indel) return;

        const int current_pos = (pam_direction == -1) ? -(q_idx + 1) : (q_idx + 1);

        if (t_idx < t_len) {
            const char q_char = q_seq[static_cast<size_t>(q_idx)];
            const char t_char = t_seq[static_cast<size_t>(t_idx)];
            const bool is_match = (q_char == t_char);
            const int new_mm = mm_cnt + (is_match ? 0 : 1);

            bool should_prune = false;
            if (has_indel) {
                if (new_mm + id_cnt > max_indel) should_prune = true;
            } else {
                if (new_mm > max_mismatch) should_prune = true;
            }

            if (!should_prune) {
                path.push_back({is_match ? "match" : "sub", current_pos, q_char, t_char});
                self(self, q_idx + 1, t_idx + 1, new_mm, id_cnt, has_indel);
                path.pop_back();
            }
        }

        if (id_cnt < max_indel) {
            const char q_char = q_seq[static_cast<size_t>(q_idx)];
            path.push_back({"del", current_pos, q_char, '-'});
            self(self, q_idx + 1, t_idx, mm_cnt, id_cnt + 1, true);
            path.pop_back();
        }

        if (t_idx < t_len && id_cnt < max_indel) {
            const char t_char = t_seq[static_cast<size_t>(t_idx)];
            path.push_back({"ins", current_pos, '-', t_char});
            self(self, q_idx, t_idx + 1, mm_cnt, id_cnt + 1, true);
            path.pop_back();
        }
    };

    solve(solve, 0, 0, 0, 0, false);
    return results;
}

std::vector<AlignmentEvent> filter_alignment_events(const std::vector<AlignmentEvent>& events) {
    std::vector<AlignmentEvent> filtered;
    filtered.reserve(events.size());
    for (const auto& event : events) {
        if (event.op != "match") filtered.push_back(event);
    }
    return filtered;
}

double predict_off_target_cpp(const std::vector<AlignmentEvent>& mismatch_positions) {
    double z = 1.0;
    int last = 99;
    std::vector<AlignmentEvent> sorted_positions = mismatch_positions;
    std::sort(sorted_positions.begin(), sorted_positions.end(), event_tuple_less);
    const int mis_num = static_cast<int>(sorted_positions.size());

    for (const auto& event : sorted_positions) {
        const int i = std::abs(event.pos);
        if (event.op != "sub") {
            z = z * (p_pos(i) - 0.23 - i * 0.01) / static_cast<double>(mis_num);
        } else {
            if (mis_num == 3) {
                z = z * p_pos(i) / 1.8;
            } else if (mis_num == 4) {
                z = z * p_pos(i) / 4.0;
            } else if (mis_num == 2) {
                z = z * p_pos(i) / 1.8;
            } else if (mis_num >= 5) {
                z = 0.0000001;
            } else {
                z = z * p_pos(i);
            }
        }

        if (std::abs(i - last) == 1) z *= 1.15;
        last = i;
    }
    return z;
}

std::string alignment_events_to_json(const std::vector<AlignmentEvent>& events) {
    std::string json;
    json.push_back('[');
    for (size_t i = 0; i < events.size(); ++i) {
        if (i > 0) json.push_back(',');
        const auto& event = events[i];
        json.append("[\"");
        json.append(event.op);
        json.append("\",");
        json.append(std::to_string(event.pos));
        json.append(",\"");
        json.push_back(event.query_base);
        json.append("\",\"");
        json.push_back(event.target_base);
        json.append("\"]");
    }
    json.push_back(']');
    return json;
}

std::string format_score_cpp(double score) {
    std::ostringstream oss;
    oss << std::setprecision(12) << std::defaultfloat << score;
    std::string out = oss.str();
    if (out.find_first_of(".eE") == std::string::npos) out += ".0";
    return out;
}

std::vector<ScoredResult> score_results_cpp(const std::vector<CombinedResult>& combined_results,
                                            const GenomeReference& target_genome,
                                            int max_indel,
                                            int max_mismatch,
                                            int pam_direction,
                                            int num_threads,
                                            size_t& discarded_count) {
    std::vector<ScoredResult> scored_slots(combined_results.size());
    std::vector<uint8_t> keep(combined_results.size(), 0);
    size_t discarded_local = 0;

    #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(guided, 1) reduction(+:discarded_local)
    #endif
    for (int idx = 0; idx < static_cast<int>(combined_results.size()); ++idx) {
        const auto& raw = combined_results[static_cast<size_t>(idx)];
        RegionHeader target_region;
        if (!parse_region_header(raw.target_name, target_region)) {
            ++discarded_local;
            continue;
        }

        const std::string target_seq = extract_padded_target_sequence(target_genome,
                                                                      target_region,
                                                                      pam_direction,
                                                                      max_indel);
        if (target_seq.empty()) {
            ++discarded_local;
            continue;
        }

        const auto all_paths = align_left_fixed_cpp(target_seq,
                                                    raw.query_seq,
                                                    max_indel,
                                                    max_mismatch,
                                                    pam_direction);
        if (all_paths.empty()) {
            ++discarded_local;
            continue;
        }

        bool has_best = false;
        double best_score = -1.0;
        std::vector<AlignmentEvent> best_details;
        for (const auto& path : all_paths) {
            std::vector<AlignmentEvent> filtered = filter_alignment_events(path);
            const double score = predict_off_target_cpp(filtered);
            if (!has_best || score > best_score) {
                best_score = score;
                best_details = std::move(filtered);
                has_best = true;
            }
            if (score == 1.0) break;
        }

        if (!has_best) {
            ++discarded_local;
            continue;
        }

        ScoredResult result;
        result.base = raw;
        result.score = best_score;
        result.details = std::move(best_details);
        if (result.details.empty()) {
            result.base.distance = 0;
            result.base.search_type = "Exact";
        } else {
            result.base.distance = static_cast<int>(result.details.size());
            result.base.search_type = event_has_indel(result.details) ? "Indel" : "Mismatch";
        }

        scored_slots[static_cast<size_t>(idx)] = std::move(result);
        keep[static_cast<size_t>(idx)] = 1;
    }

    discarded_count = discarded_local;
    std::vector<ScoredResult> compact_results;
    compact_results.reserve(combined_results.size() - discarded_count);
    for (size_t i = 0; i < scored_slots.size(); ++i) {
        if (keep[i]) compact_results.push_back(std::move(scored_slots[i]));
    }
    return compact_results;
}

void write_scored_results_tsv(const std::string& output_path,
                              const std::vector<ScoredResult>& results) {
    std::ofstream out_stream(output_path, std::ios::binary | std::ios::trunc);
    out_stream << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\tSearchType\tGenomicFrequency\tScore\tDetails\n";
    for (const auto& result : results) {
        out_stream << result.base.query_name << '\t'
                   << result.base.query_seq << '\t'
                   << result.base.target_name << '\t'
                   << result.base.target_seq << '\t'
                   << result.base.distance << '\t'
                   << result.base.search_type << '\t'
                   << result.base.genomic_frequency << '\t'
                   << format_score_cpp(result.score) << '\t'
                   << alignment_events_to_json(result.details) << '\n';
    }
    out_stream.flush();
}
