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
// indel_search.cpp - OPTIMIZED VERSION with Banded DP
#include "indel_search.hpp"
#include "fasta_parser.hpp"
#include "progress.hpp"
#include "utils.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <omp.h>
#include <iomanip>
#include <cstring>
#include <atomic>
#include <cstddef>
#include <stdexcept>

#ifdef _MSC_VER
    #include <intrin.h>
#endif

// ============================================================================
// OPTIMIZATION: Thread-local pre-allocated buffers for Banded DP
// ============================================================================
thread_local std::array<int, 128> dp_buffer_a;
thread_local std::array<int, 128> dp_buffer_b;

// Forward declarations
static char complement(char c);
static void append_result_line(std::string& buffer,
                               const std::string& query_name,
                               const std::string& query_seq,
                               const std::string& target_name,
                               const std::string& target_seq,
                               int distance);
static void get_sequence_from_loc(const SpacerLocation& loc,
                                  const GenomeReference& genome,
                                  std::string& out);
static void encode_sequence_channels(std::string_view seq,
                                     bool reverse_complement,
                                     uint64_t& ch_a,
                                     uint64_t& ch_c,
                                     uint64_t& ch_g,
                                     uint64_t& ch_t);

static unsigned count_trailing_zeros_16(uint16_t value) {
#ifdef _MSC_VER
    unsigned long idx = 0;
    _BitScanForward(&idx, static_cast<unsigned long>(value));
    return static_cast<unsigned>(idx);
#else
    return static_cast<unsigned>(__builtin_ctz(static_cast<unsigned>(value)));
#endif
}

// ============================================================================
// OPTIMIZATION: Banded DP implementation of right-anchored edit distance
// ============================================================================
int right_anchored_distance(const std::string_view& Q, const std::string_view& T, int max_dist) {
    int m = static_cast<int>(Q.size());
    int n = static_cast<int>(T.size());
    const int inf = max_dist + 1;
    
    if (m == 0) return n;
    if (n == 0) return m;
    
    // Early rejection: if length difference > max_dist, can't be within threshold
    if (std::abs(m - n) > max_dist) return max_dist + 1;
    
    // Fall back to standard allocation for very long sequences
    if (n >= 128) {
        std::vector<int> prev(n + 1, 0), curr(n + 1, inf);
        for (int j = 0; j <= n; ++j) prev[j] = 0;

        for (int i = 1; i <= m; ++i) {
            int j_min = std::max(1, i - max_dist);
            int j_max = std::min(n, i + max_dist);
            std::fill(curr.begin(), curr.end(), inf);
            curr[0] = 0;
            
            for (int j = j_min; j <= j_max; ++j) {
                int cost = (Q[i - 1] == T[j - 1]) ? 0 : 1;
                curr[j] = std::min({prev[j - 1] + cost, prev[j] + 1, curr[j - 1] + 1});
            }
            std::swap(prev, curr);
        }
        return prev[n];
    }

    // Fast path: use pre-allocated buffers
    int* prev = dp_buffer_a.data();
    int* curr = dp_buffer_b.data();
    
    // Initialize first row (right-anchored: free prefix deletion)
    for (int j = 0; j <= n; ++j) prev[j] = 0;

    for (int i = 1; i <= m; ++i) {
        // Compute valid column range for this row (banded)
        int j_min = std::max(1, i - max_dist);
        int j_max = std::min(n, i + max_dist);

        // Important: initialize all states so out-of-band cells are never read uninitialized.
        for (int j = 0; j <= n; ++j) curr[j] = inf;
        curr[0] = 0;  // First column: free prefix deletion
        
        for (int j = j_min; j <= j_max; ++j) {
            int cost = (Q[i - 1] == T[j - 1]) ? 0 : 1;
            int subst = prev[j - 1] + cost;
            int del = prev[j] + 1;
            int ins = curr[j - 1] + 1;
            curr[j] = std::min({subst, del, ins});
        }
        
        // Early termination: if entire row > max_dist, we can abort
        bool all_above = true;
        for (int j = j_min; j <= j_max && all_above; ++j) {
            if (curr[j] <= max_dist) all_above = false;
        }
        if (all_above && j_min <= j_max) return max_dist + 1;
        
        std::swap(prev, curr);
    }
    
    return prev[n];
}

// ============================================================================
// Binary Channel Filter Implementation
// ============================================================================

static void dna_to_binary_channels(const std::string_view& seq,
                                   uint64_t& channel_a, uint64_t& channel_c,
                                   uint64_t& channel_g, uint64_t& channel_t) {
    channel_a = 0;
    channel_c = 0;
    channel_g = 0;
    channel_t = 0;

    if (seq.length() > 64) {
        return;
    }

    for (char base : seq) {
        channel_a <<= 1;
        channel_c <<= 1;
        channel_g <<= 1;
        channel_t <<= 1;
        switch (std::toupper(base)) {
            case 'A': channel_a |= 1; break;
            case 'C': channel_c |= 1; break;
            case 'G': channel_g |= 1; break;
            case 'T': channel_t |= 1; break;
            default: break;
        }
    }
}

Channel BinaryChannelFilter::apply_substitution(const Channel& ch, size_t pos, char new_base) const {
    if (pos >= ch.len) return ch;
    size_t bit_pos = ch.len - 1 - pos;
    uint64_t mask = 1ULL << bit_pos;

    Channel new_ch = ch;
    new_ch.a &= ~mask;
    new_ch.c &= ~mask;
    new_ch.g &= ~mask;
    new_ch.t &= ~mask;
    
    switch (std::toupper(new_base)) {
        case 'A': new_ch.a |= mask; break;
        case 'C': new_ch.c |= mask; break;
        case 'G': new_ch.g |= mask; break;
        case 'T': new_ch.t |= mask; break;
        default: break;
    }
    return new_ch;
}

Channel BinaryChannelFilter::apply_deletion(const Channel& ch, size_t pos) const {
    if (ch.len <= 1 || pos >= ch.len) return ch;
    size_t bit_pos = ch.len - 1 - pos;
    uint64_t mask = 1ULL << bit_pos;

    Channel new_ch;
    new_ch.len = ch.len - 1;
    
    auto apply_del_to_ch = [&](uint64_t orig_ch) -> uint64_t {
        orig_ch &= ~mask;
        uint64_t high_mask = ~((1ULL << (bit_pos + 1)) - 1ULL);
        uint64_t high_part = (orig_ch & high_mask) >> 1;
        uint64_t low_part = orig_ch & ((1ULL << bit_pos) - 1ULL);
        return high_part | low_part;
    };
    
    new_ch.a = apply_del_to_ch(ch.a);
    new_ch.c = apply_del_to_ch(ch.c);
    new_ch.g = apply_del_to_ch(ch.g);
    new_ch.t = apply_del_to_ch(ch.t);
    return new_ch;
}

Channel BinaryChannelFilter::apply_insertion(const Channel& ch, size_t pos, char new_base) const {
    size_t prefix_len = pos;
    size_t suffix_len = ch.len - pos;
    size_t new_len = ch.len + 1;
    size_t new_bit_pos = suffix_len;
    uint64_t new_mask = 1ULL << new_bit_pos;

    Channel new_ch;
    new_ch.len = new_len;
    
    auto apply_ins_to_ch = [&](uint64_t orig_ch, bool set_new) -> uint64_t {
        uint64_t prefix_part = orig_ch >> suffix_len;
        uint64_t suffix_part = orig_ch & ((1ULL << suffix_len) - 1ULL);
        uint64_t result = (prefix_part << (suffix_len + 1));
        if (set_new) result |= new_mask;
        result |= suffix_part;
        return result;
    };
    
    new_ch.a = apply_ins_to_ch(ch.a, std::toupper(new_base) == 'A');
    new_ch.c = apply_ins_to_ch(ch.c, std::toupper(new_base) == 'C');
    new_ch.g = apply_ins_to_ch(ch.g, std::toupper(new_base) == 'G');
    new_ch.t = apply_ins_to_ch(ch.t, std::toupper(new_base) == 'T');
    return new_ch;
}

void BinaryChannelFilter::add_normalized(const Channel& ch) {
    if (ch.len == 0) return;

    if (ch.len >= fixed_len) {
        uint64_t norm_a = ch.a & mask;
        uint64_t norm_c = ch.c & mask;
        uint64_t norm_g = ch.g & mask;
        uint64_t norm_t = ch.t & mask;
        insert_variant(norm_a, norm_c, norm_g, norm_t);
    } else {
        size_t d = fixed_len - ch.len;
        uint64_t low_mask = (1ULL << ch.len) - 1;
        uint64_t low_a = ch.a & low_mask;
        uint64_t low_c = ch.c & low_mask;
        uint64_t low_g = ch.g & low_mask;
        uint64_t low_t = ch.t & low_mask;

        uint64_t num_pats = 1ULL << d;
        for (uint64_t pat = 0; pat < num_pats; ++pat) {
            uint64_t high_shifted = pat << ch.len;
            insert_variant(high_shifted | low_a,
                           high_shifted | low_c,
                           high_shifted | low_g,
                           high_shifted | low_t);
        }
    }
}

void BinaryChannelFilter::insert_variant(uint64_t a, uint64_t c, uint64_t g, uint64_t t) {
    variants_A.insert(a);
    variants_C.insert(c);
    variants_G.insert(g);
    variants_T.insert(t);
}

bool BinaryChannelFilter::contains_variant(uint64_t a, uint64_t c, uint64_t g, uint64_t t) const {
    if (!variants_A.contains(a)) return false;
    if (!variants_C.contains(c)) return false;
    if (!variants_G.contains(g)) return false;
    if (!variants_T.contains(t)) return false;
    return true;
}

void BinaryChannelFilter::accumulate_query_mask_u32(uint16_t query_bit,
                                                    std::vector<uint16_t>& a_map,
                                                    std::vector<uint16_t>& c_map,
                                                    std::vector<uint16_t>& g_map,
                                                    std::vector<uint16_t>& t_map) const {
    variants_A.for_each([&](uint64_t value) {
        a_map[static_cast<uint32_t>(value)] |= query_bit;
    });
    variants_C.for_each([&](uint64_t value) {
        c_map[static_cast<uint32_t>(value)] |= query_bit;
    });
    variants_G.for_each([&](uint64_t value) {
        g_map[static_cast<uint32_t>(value)] |= query_bit;
    });
    variants_T.for_each([&](uint64_t value) {
        t_map[static_cast<uint32_t>(value)] |= query_bit;
    });
}

BinaryChannelFilter::BinaryChannelFilter(const std::string_view& query_seq, int max_edits, size_t spacer_len)
    : fixed_len(spacer_len), mask((1ULL << spacer_len) - 1) {
    constexpr size_t variant_reserve = 4096;
    constexpr size_t visited_reserve = 8192;
    constexpr float max_load = 0.7f;

    variants_A.reserve(variant_reserve);
    variants_C.reserve(variant_reserve);
    variants_G.reserve(variant_reserve);
    variants_T.reserve(variant_reserve);
    
    size_t init_len = query_seq.length();
    min_len = std::max(1, static_cast<int>(init_len) - max_edits);
    max_len = static_cast<int>(init_len) + max_edits;

    uint64_t init_a, init_c, init_g, init_t;
    dna_to_binary_channels(query_seq, init_a, init_c, init_g, init_t);

    Channel start_ch{init_a, init_c, init_g, init_t, init_len};
    std::unordered_set<Channel, ChannelHash> visited;
    visited.max_load_factor(max_load);
    visited.reserve(visited_reserve);
    std::queue<std::pair<Channel, int>> q;

    q.push({start_ch, 0});
    visited.insert(start_ch);
    add_normalized(start_ch);

    const char bases[] = {'A', 'C', 'G', 'T'};

    while (!q.empty()) {
        auto [current_ch, edits] = q.front();
        q.pop();

        if (edits >= max_edits) continue;

        // Substitutions
        for (size_t p = 0; p < current_ch.len; ++p) {
            for (char b : bases) {
                Channel next_ch = apply_substitution(current_ch, p, b);
                if (visited.insert(next_ch).second) {
                    add_normalized(next_ch);
                    q.push({next_ch, edits + 1});
                }
            }
        }

        // Deletions
        if (current_ch.len > 1) {
            for (size_t p = 0; p < current_ch.len; ++p) {
                Channel next_ch = apply_deletion(current_ch, p);
                if (visited.insert(next_ch).second) {
                    add_normalized(next_ch);
                    q.push({next_ch, edits + 1});
                }
            }
        }

        // Insertions
        for (size_t p = 0; p <= current_ch.len; ++p) {
            for (char b : bases) {
                Channel next_ch = apply_insertion(current_ch, p, b);
                if (visited.insert(next_ch).second) {
                    add_normalized(next_ch);
                    q.push({next_ch, edits + 1});
                }
            }
        }
    }

}

bool BinaryChannelFilter::check(const std::string_view& target_seq) const {
    size_t target_len = target_seq.length();
    bool length_ok = (target_len >= static_cast<size_t>(min_len) && target_len <= static_cast<size_t>(max_len));

    if (!length_ok) {
        return false;
    }

    uint64_t t_ch_a, t_ch_c, t_ch_g, t_ch_t;
    dna_to_binary_channels(target_seq, t_ch_a, t_ch_c, t_ch_g, t_ch_t);

    uint64_t norm_a = t_ch_a & mask;
    uint64_t norm_c = t_ch_c & mask;
    uint64_t norm_g = t_ch_g & mask;
    uint64_t norm_t = t_ch_t & mask;

    return contains_variant(norm_a, norm_c, norm_g, norm_t);
}

// OPTIMIZATION: Check pre-encoded target without string conversion
bool BinaryChannelFilter::check_encoded(uint64_t t_a, uint64_t t_c, uint64_t t_g, uint64_t t_t, size_t target_len) const {
    bool length_ok = (target_len >= static_cast<size_t>(min_len) && target_len <= static_cast<size_t>(max_len));
    if (!length_ok) return false;

    return check_encoded_fixed(t_a, t_c, t_g, t_t);
}

bool BinaryChannelFilter::check_encoded_fixed(uint64_t t_a, uint64_t t_c, uint64_t t_g, uint64_t t_t) const {
    return contains_variant(t_a, t_c, t_g, t_t);
}

// Check pre-encoded target against filter (optimized - no string conversion)
bool check_preencoded(const BinaryChannelFilter& filter,
                      const PreencodedTargets& targets,
                      size_t idx) {
    if (targets.use_u32) {
        return filter.check_encoded_fixed(targets.a32[idx], targets.c32[idx], targets.g32[idx], targets.t32[idx]);
    }
    return filter.check_encoded_fixed(targets.a64[idx], targets.c64[idx], targets.g64[idx], targets.t64[idx]);
}

static void encode_sequence_channels(std::string_view seq,
                                     bool reverse_complement,
                                     uint64_t& ch_a,
                                     uint64_t& ch_c,
                                     uint64_t& ch_g,
                                     uint64_t& ch_t) {
    ch_a = 0;
    ch_c = 0;
    ch_g = 0;
    ch_t = 0;

    auto encode_base = [&](char base) {
        ch_a <<= 1;
        ch_c <<= 1;
        ch_g <<= 1;
        ch_t <<= 1;
        switch (std::toupper(base)) {
            case 'A': ch_a |= 1; break;
            case 'C': ch_c |= 1; break;
            case 'G': ch_g |= 1; break;
            case 'T': ch_t |= 1; break;
            default: break;
        }
    };

    if (reverse_complement) {
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            encode_base(complement(*it));
        }
        return;
    }

    for (char base : seq) {
        encode_base(base);
    }
}

static void append_result_line(std::string& buffer,
                               const std::string& query_name,
                               const std::string& query_seq,
                               const std::string& target_name,
                               const std::string& target_seq,
                               int distance) {
    buffer.append(query_name);
    buffer.push_back('\t');
    buffer.append(query_seq);
    buffer.push_back('\t');
    buffer.append(target_name);
    buffer.push_back('\t');
    buffer.append(target_seq);
    buffer.push_back('\t');
    buffer.append(std::to_string(distance));
    buffer.push_back('\n');
}

// OPTIMIZATION: Pre-encode all targets for indel search (parallel)
PreencodedTargets preencode_targets_for_indel(
    const std::vector<SpacerLocation>& target_locs,
    const GenomeReference& target_genome,
    int num_threads) {
    PreencodedTargets result;
    bool use_u32 = true;
    for (const auto& loc : target_locs) {
        if (loc.length > 32) {
            use_u32 = false;
            break;
        }
    }

    result.use_u32 = use_u32;
    result.count = target_locs.size();
    if (use_u32) {
        result.a32.reset(new uint32_t[target_locs.size()]);
        result.c32.reset(new uint32_t[target_locs.size()]);
        result.g32.reset(new uint32_t[target_locs.size()]);
        result.t32.reset(new uint32_t[target_locs.size()]);
    } else {
        result.a64.reset(new uint64_t[target_locs.size()]);
        result.c64.reset(new uint64_t[target_locs.size()]);
        result.g64.reset(new uint64_t[target_locs.size()]);
        result.t64.reset(new uint64_t[target_locs.size()]);
    }
    
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(target_locs.size()); ++i) {
            const size_t idx = static_cast<size_t>(i);
            const auto& loc = target_locs[idx];
            std::string_view seq = target_genome.get_sequence(loc);

            // Encode to 4 channels
            uint64_t ch_a, ch_c, ch_g, ch_t;
            encode_sequence_channels(seq, loc.is_reverse, ch_a, ch_c, ch_g, ch_t);
            if (use_u32) {
                result.a32[idx] = static_cast<uint32_t>(ch_a);
                result.c32[idx] = static_cast<uint32_t>(ch_c);
                result.g32[idx] = static_cast<uint32_t>(ch_g);
                result.t32[idx] = static_cast<uint32_t>(ch_t);
            } else {
                result.a64[idx] = ch_a;
                result.c64[idx] = ch_c;
                result.g64[idx] = ch_g;
                result.t64[idx] = ch_t;
            }
        }
    }
    
    return result;
}

// Helper: complement a single base
char complement(char c) {
    switch (std::toupper(c)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

// Helper: reverse complement
static std::string get_reverse_complement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.length());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (std::toupper(*it)) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            default: rc += 'N'; break;
        }
    }
    return rc;
}

static std::string get_sequence_from_loc(const SpacerLocation& loc, const GenomeReference& genome) {
    std::string_view seq = genome.get_sequence(loc);
    if (loc.is_reverse) {
        return get_reverse_complement(seq);
    }
    return std::string(seq);
}

static void get_sequence_from_loc(const SpacerLocation& loc,
                                  const GenomeReference& genome,
                                  std::string& out) {
    std::string_view seq = genome.get_sequence(loc);
    out.clear();
    if (seq.empty()) return;

    out.reserve(seq.length());
    if (loc.is_reverse) {
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            out.push_back(complement(*it));
        }
        return;
    }

    out.append(seq.data(), seq.size());
}

// Legacy interface
void process_indel_search(const std::string& file_a, const std::string& file_b,
                          const std::string& output_filename, int max_indel,
                          int num_threads, int spacer_len) {
    SequenceList sequences_b = parse_fasta(file_b);
    SequenceList sequences_a = parse_fasta(file_a);

    std::ofstream out_file(output_filename);
    if (!out_file.is_open()) return;
    std::vector<char> out_file_buffer(1 << 20);
    out_file.rdbuf()->pubsetbuf(out_file_buffer.data(), static_cast<std::streamsize>(out_file_buffer.size()));
    out_file << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\n";

    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(sequences_a.size()); ++i) {
        const size_t query_idx = static_cast<size_t>(i);
        const auto& query_seq = sequences_a.get(query_idx);
        BinaryChannelFilter filter(query_seq, max_indel, spacer_len);

        std::vector<SequenceID> filtered_candidates;
        for (size_t j = 0; j < sequences_b.size(); ++j) {
            if (filter.check(sequences_b.get(j))) {
                filtered_candidates.push_back(j);
            }
        }

        std::vector<MatchResult> local_results;
        for (const auto& candidate_id : filtered_candidates) {
            const auto& target_seq = sequences_b.get(candidate_id);
            int dist = right_anchored_distance(query_seq, target_seq, max_indel);
            if (dist <= max_indel) {
                const auto& target_name = sequences_b.get_header(candidate_id);
                local_results.emplace_back(target_name, target_seq, dist);
            }
        }

        if (!local_results.empty()) {
            std::stringstream task_buffer;
            const auto& query_name = sequences_a.headers[query_idx];
            for (const auto& res : local_results) {
                task_buffer << query_name << "\t" << query_seq << "\t"
                            << std::get<0>(res) << "\t"
                            << std::get<1>(res) << "\t"
                            << std::get<2>(res) << "\n";
            }
            #pragma omp critical
            {
                out_file << task_buffer.str();
            }
        }
    }
    out_file.flush();
}

// Location-based interface
std::vector<SearchHit> process_indel_search_loc_collect_preencoded(const std::vector<SpacerLocation>& query_locs,
                                                                   const std::vector<SpacerLocation>& target_locs,
                                                                   const GenomeReference& query_genome,
                                                                   const GenomeReference& target_genome,
                                                                   const PreencodedTargets& preencoded_targets,
                                                                   int max_indel,
                                                                   int num_threads,
                                                                   int spacer_len,
                                                                   size_t progress_base,
                                                                   size_t progress_total) {
    if (query_locs.empty() || target_locs.empty()) return {};

    const size_t query_batch_size = 16;
    const bool mask_fusion_enabled = preencoded_targets.use_u32 && spacer_len > 0 && spacer_len <= 20;
    const size_t total_queries = query_locs.size();
    const size_t progress_denominator = (progress_total == 0) ? total_queries : progress_total;

    std::vector<std::vector<SearchHit>> thread_hits(num_threads);
    std::vector<uint16_t> batch_mask_a;
    std::vector<uint16_t> batch_mask_c;
    std::vector<uint16_t> batch_mask_g;
    std::vector<uint16_t> batch_mask_t;
    if (mask_fusion_enabled) {
        const size_t value_space = 1ULL << spacer_len;
        batch_mask_a.resize(value_space);
        batch_mask_c.resize(value_space);
        batch_mask_g.resize(value_space);
        batch_mask_t.resize(value_space);
    }

    for (size_t batch_start = 0; batch_start < total_queries; batch_start += query_batch_size) {
        const size_t batch_count = std::min(query_batch_size, query_locs.size() - batch_start);

        std::vector<std::string> batch_query_seqs(batch_count);
        std::vector<std::string> batch_query_names(batch_count);
        std::vector<BinaryChannelFilter> batch_filters;
        batch_filters.reserve(batch_count);

        for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
            const auto& query_loc = query_locs[batch_start + local_idx];
            batch_query_seqs[local_idx] = get_sequence_from_loc(query_loc, query_genome);
            batch_query_names[local_idx] = query_genome.format_header(query_loc);
            batch_filters.emplace_back(batch_query_seqs[local_idx], max_indel, spacer_len);
        }

        std::vector<std::vector<std::vector<size_t>>> thread_candidates(
            num_threads, std::vector<std::vector<size_t>>(batch_count));
        const bool use_mask_fusion = mask_fusion_enabled && batch_count <= 16;
        if (use_mask_fusion) {
            std::fill(batch_mask_a.begin(), batch_mask_a.end(), 0);
            std::fill(batch_mask_c.begin(), batch_mask_c.end(), 0);
            std::fill(batch_mask_g.begin(), batch_mask_g.end(), 0);
            std::fill(batch_mask_t.begin(), batch_mask_t.end(), 0);

            for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
                const uint16_t query_bit = static_cast<uint16_t>(1u << local_idx);
                batch_filters[local_idx].accumulate_query_mask_u32(
                    query_bit,
                    batch_mask_a,
                    batch_mask_c,
                    batch_mask_g,
                    batch_mask_t);
            }
        }

        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            auto& local_candidates = thread_candidates[tid];
            const size_t num_targets_local = preencoded_targets.size();

            if (preencoded_targets.use_u32) {
                const uint32_t* target_a = preencoded_targets.a32.get();
                const uint32_t* target_c = preencoded_targets.c32.get();
                const uint32_t* target_g = preencoded_targets.g32.get();
                const uint32_t* target_t = preencoded_targets.t32.get();

                if (use_mask_fusion) {
                    #pragma omp for schedule(static)
                    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(num_targets_local); ++j) {
                        const size_t target_idx = static_cast<size_t>(j);
                        uint16_t remaining =
                            batch_mask_a[target_a[target_idx]] &
                            batch_mask_c[target_c[target_idx]] &
                            batch_mask_g[target_g[target_idx]] &
                            batch_mask_t[target_t[target_idx]];

                        while (remaining != 0) {
                            const unsigned bit = count_trailing_zeros_16(remaining);
                            local_candidates[bit].push_back(target_idx);
                            remaining &= static_cast<uint16_t>(remaining - 1);
                        }
                    }
                } else {
                    #pragma omp for schedule(static)
                    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(num_targets_local); ++j) {
                        const size_t target_idx = static_cast<size_t>(j);
                        const uint32_t a = target_a[target_idx];
                        const uint32_t c = target_c[target_idx];
                        const uint32_t g = target_g[target_idx];
                        const uint32_t t = target_t[target_idx];
                        for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
                            if (!batch_filters[local_idx].check_encoded_fixed(a, c, g, t)) continue;
                            local_candidates[local_idx].push_back(target_idx);
                        }
                    }
                }
            } else {
                const uint64_t* target_a = preencoded_targets.a64.get();
                const uint64_t* target_c = preencoded_targets.c64.get();
                const uint64_t* target_g = preencoded_targets.g64.get();
                const uint64_t* target_t = preencoded_targets.t64.get();

                #pragma omp for schedule(static)
                    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(num_targets_local); ++j) {
                        const size_t target_idx = static_cast<size_t>(j);
                        const uint64_t a = target_a[target_idx];
                        const uint64_t c = target_c[target_idx];
                        const uint64_t g = target_g[target_idx];
                        const uint64_t t = target_t[target_idx];
                        for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
                            if (!batch_filters[local_idx].check_encoded_fixed(a, c, g, t)) continue;
                            local_candidates[local_idx].push_back(target_idx);
                        }
                    }
                }
        }

        std::vector<std::vector<size_t>> candidate_ids(batch_count);
        for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
            size_t total_candidates = 0;
            for (int tid = 0; tid < num_threads; ++tid) {
                total_candidates += thread_candidates[tid][local_idx].size();
            }
            auto& merged = candidate_ids[local_idx];
            merged.reserve(total_candidates);
            for (int tid = 0; tid < num_threads; ++tid) {
                auto& local = thread_candidates[tid][local_idx];
                merged.insert(merged.end(), local.begin(), local.end());
            }
        }

        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            auto& local_hits = thread_hits[tid];
            std::string target_seq_buffer;
            target_seq_buffer.reserve(spacer_len + max_indel);

            #pragma omp for schedule(guided, 1)
            for (std::ptrdiff_t local_idx = 0; local_idx < static_cast<std::ptrdiff_t>(batch_count); ++local_idx) {
                const size_t batch_idx = static_cast<size_t>(local_idx);
                const std::string& query_seq = batch_query_seqs[batch_idx];
                const std::string& query_name = batch_query_names[batch_idx];
                const auto& candidates = candidate_ids[batch_idx];

                for (size_t target_idx : candidates) {
                    const auto& target_loc = target_locs[target_idx];
                    get_sequence_from_loc(target_loc, target_genome, target_seq_buffer);
                    if (target_seq_buffer.empty()) {
                        continue;
                    }

                    int dist = right_anchored_distance(query_seq, target_seq_buffer, max_indel);
                    if (dist <= max_indel) {
                        SearchHit hit;
                        hit.query_name = query_name;
                        hit.query_seq = query_seq;
                        hit.target_name = target_genome.format_header(target_loc);
                        hit.target_seq = target_seq_buffer;
                        hit.distance = dist;
                        local_hits.push_back(std::move(hit));
                    }
                }
            }
        }

        progress::progress_update("indel",
                                  progress_base + batch_start + batch_count,
                                  progress_denominator,
                                  "queries");
    }

    size_t total_hits = 0;
    for (const auto& local_hits : thread_hits) total_hits += local_hits.size();
    std::vector<SearchHit> collected_hits;
    collected_hits.reserve(total_hits);
    for (auto& local_hits : thread_hits) {
        collected_hits.insert(collected_hits.end(),
                              std::make_move_iterator(local_hits.begin()),
                              std::make_move_iterator(local_hits.end()));
    }
    return collected_hits;
}

std::vector<SearchHit> process_indel_search_loc_collect(const std::vector<SpacerLocation>& query_locs,
                                                        const std::vector<SpacerLocation>& target_locs,
                                                        const GenomeReference& query_genome,
                                                        const GenomeReference& target_genome,
                                                        int max_indel,
                                                        int num_threads,
                                                        int spacer_len) {
    progress::log("Pre-encoding " + std::to_string(target_locs.size()) + " target sequences for indel search...");
    const double preencode_start = omp_get_wtime();
    auto preencoded_targets = preencode_targets_for_indel(target_locs, target_genome, num_threads);
    const double preencode_time = omp_get_wtime() - preencode_start;
    progress::log("Target pre-encoding complete.");
    progress::log("  Using " + std::to_string(num_threads) +
                  " threads for indel search (query batch size 16)");

    const double search_start = omp_get_wtime();
    std::vector<SearchHit> collected_hits =
        process_indel_search_loc_collect_preencoded(query_locs,
                                                    target_locs,
                                                    query_genome,
                                                    target_genome,
                                                    preencoded_targets,
                                                    max_indel,
                                                    num_threads,
                                                    spacer_len,
                                                    0,
                                                    query_locs.size());
    const double search_time = omp_get_wtime() - search_start;

    {
        std::ostringstream oss;
        oss << "  Indel phase breakdown:";
        progress::log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "    Preencode wall time: " << preencode_time << "s";
        progress::log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "    Filter+scan + DP wall time: " << search_time << "s";
        progress::log(oss.str());
    }
    progress::log("    Writeback wall time: 0s (in-memory collect mode)");
    {
        std::ostringstream oss;
        oss << "    Collected hits: " << collected_hits.size();
        progress::log(oss.str());
    }
    return collected_hits;
}

void process_indel_search_loc(const std::vector<SpacerLocation>& query_locs,
                              const std::vector<SpacerLocation>& target_locs,
                              const GenomeReference& query_genome,
                              const GenomeReference& target_genome,
                              const std::string& output_filename,
                              int max_indel,
                              int num_threads,
                              int spacer_len) {
    std::vector<SearchHit> collected_hits =
        process_indel_search_loc_collect(query_locs,
                                         target_locs,
                                         query_genome,
                                         target_genome,
                                         max_indel,
                                         num_threads,
                                         spacer_len);

    std::ofstream out_file(output_filename, std::ios::binary | std::ios::trunc);
    if (!out_file.is_open()) return;
    std::vector<char> out_file_buffer(1 << 20);
    out_file.rdbuf()->pubsetbuf(out_file_buffer.data(), static_cast<std::streamsize>(out_file_buffer.size()));
    out_file << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\n";
    for (const auto& hit : collected_hits) {
        out_file << hit.query_name << '\t'
                 << hit.query_seq << '\t'
                 << hit.target_name << '\t'
                 << hit.target_seq << '\t'
                 << hit.distance << '\n';
    }
    out_file.flush();
}
