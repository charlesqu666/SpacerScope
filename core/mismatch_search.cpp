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
// Strategy A Simple ONLY: 2-bit Encoding Brute-force Search
// NO strategy switching - always uses brute-force for fair comparison testing

#include "mismatch_search.hpp"
#include "fasta_parser.hpp"
#include "progress.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#ifdef _MSC_VER
    #include <intrin.h>
    #define POPCOUNT64(x) __popcnt64(x)
#elif defined(__GNUC__) || defined(__clang__)
    #define POPCOUNT64(x) __builtin_popcountll(x)
#else
    inline int popcount64_fallback(uint64_t x) {
        x = x - ((x >> 1) & 0x5555555555555555ULL);
        x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
        x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
        return static_cast<int>((x * 0x0101010101010101ULL) >> 56);
    }
    #define POPCOUNT64(x) popcount64_fallback(x)
#endif

static char complement(char c) {
    switch (std::toupper(c)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

static std::string get_rc(const std::string_view& seq) {
    std::string rc;
    rc.reserve(seq.length());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) rc.push_back(complement(*it));
    return rc;
}

static void materialize_oriented_sequence(std::string_view seq, bool reverse_complement, std::string& out) {
    out.clear();
    out.reserve(seq.size());
    if (!reverse_complement) {
        out.append(seq.data(), seq.size());
        return;
    }

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) out.push_back(complement(*it));
}

// 2-bit Encoding: A=00, C=01, G=10, T=11
encoded_seq_t encode_simple(std::string_view seq) {
    encoded_seq_t res = 0;
    for (size_t i = 0; i < seq.length() && i < 32; ++i) {
        res <<= 2;
        switch (std::toupper(seq[i])) {
            case 'A': res |= 0; break;
            case 'C': res |= 1; break;
            case 'G': res |= 2; break;
            case 'T': res |= 3; break;
            default:  res |= 0; break;
        }
    }
    return res;
}

static encoded_seq_t encode_simple_oriented(std::string_view seq, bool reverse_complement) {
    if (!reverse_complement) return encode_simple(seq);

    encoded_seq_t res = 0;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        res <<= 2;
        switch (std::toupper(complement(*it))) {
            case 'A': res |= 0; break;
            case 'C': res |= 1; break;
            case 'G': res |= 2; break;
            case 'T': res |= 3; break;
            default:  res |= 0; break;
        }
    }
    return res;
}

inline int hamming_distance_simple(encoded_seq_t q, encoded_seq_t t, int len) {
    (void)len;
    const encoded_seq_t diff = q ^ t;
    const encoded_seq_t mask = 0x5555555555555555ULL;
    const encoded_seq_t has_diff = (diff | (diff >> 1)) & mask;
    return POPCOUNT64(has_diff);
}

static void append_mismatch_line(std::string& buffer,
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

struct CodeUniquenessEstimate {
    size_t sample_size = 0;
    size_t unique_count = 0;
    bool exact = false;
};

static CodeUniquenessEstimate estimate_code_uniqueness(const std::vector<encoded_seq_t>& codes) {
    CodeUniquenessEstimate estimate;
    if (codes.empty()) return estimate;

    static constexpr size_t EXACT_LIMIT = 2000000;
    static constexpr size_t SAMPLE_LIMIT = 1000000;

    std::vector<encoded_seq_t> sampled_codes;
    if (codes.size() <= EXACT_LIMIT) {
        sampled_codes = codes;
        estimate.exact = true;
    } else {
        const size_t sample_size = std::min(SAMPLE_LIMIT, codes.size());
        sampled_codes.reserve(sample_size);
        for (size_t i = 0; i < sample_size; ++i) {
            const size_t idx = (i * codes.size()) / sample_size;
            sampled_codes.push_back(codes[idx]);
        }
    }

    std::sort(sampled_codes.begin(), sampled_codes.end());
    auto unique_end = std::unique(sampled_codes.begin(), sampled_codes.end());
    estimate.sample_size = sampled_codes.size();
    estimate.unique_count = static_cast<size_t>(std::distance(sampled_codes.begin(), unique_end));
    return estimate;
}

std::vector<SearchHit> process_mismatch_search_loc_collect(
    const std::vector<SpacerLocation>& query_locs,
    const std::vector<SpacerLocation>& target_locs,
    const GenomeReference& query_genome,
    const GenomeReference& target_genome,
    int max_mismatch,
    int num_threads
) {
    if (target_locs.empty() || query_locs.empty()) return {};

    progress::log("[SIMPLE 2-BIT VERSION] Always using brute-force with 2-bit encoding (m=" +
                  std::to_string(max_mismatch) + ")");
    MismatchPreencodedTargets preencoded_targets =
        preencode_targets_for_mismatch(target_locs, target_genome, num_threads);
    return process_mismatch_search_loc_collect_preencoded(query_locs,
                                                          target_locs,
                                                          preencoded_targets,
                                                          query_genome,
                                                          target_genome,
                                                          max_mismatch,
                                                          num_threads);
}

MismatchPreencodedTargets preencode_targets_for_mismatch(const std::vector<SpacerLocation>& target_locs,
                                                         const GenomeReference& target_genome,
                                                         int num_threads) {
    MismatchPreencodedTargets preencoded;
    if (target_locs.empty()) return preencoded;

    preencoded.fixed_len_fast_path = true;
    preencoded.fixed_target_len = target_locs.front().length;
    if (preencoded.fixed_target_len == 0 || preencoded.fixed_target_len > 32) {
        preencoded.fixed_len_fast_path = false;
    } else {
        for (const auto& loc : target_locs) {
            if (loc.length != preencoded.fixed_target_len) {
                preencoded.fixed_len_fast_path = false;
                break;
            }
        }
    }

    progress::log("Pre-encoding " + std::to_string(target_locs.size()) + " target sequences with 2-bit...");
    preencoded.target_codes.resize(target_locs.size());
    if (!preencoded.fixed_len_fast_path) {
        preencoded.target_lengths.resize(target_locs.size());
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < static_cast<int>(target_locs.size()); ++idx) {
        const auto& loc = target_locs[static_cast<size_t>(idx)];
        const std::string_view seq_view = target_genome.get_sequence(loc);
        if (seq_view.empty() || seq_view.length() > 32) {
            preencoded.target_codes[static_cast<size_t>(idx)] = 0;
            if (!preencoded.fixed_len_fast_path) {
                preencoded.target_lengths[static_cast<size_t>(idx)] = 0;
            }
            continue;
        }

        preencoded.target_codes[static_cast<size_t>(idx)] = encode_simple_oriented(seq_view, loc.is_reverse);
        if (!preencoded.fixed_len_fast_path) {
            preencoded.target_lengths[static_cast<size_t>(idx)] = static_cast<int>(seq_view.length());
        }
    }

    if (preencoded.fixed_len_fast_path) {
        const CodeUniquenessEstimate uniqueness = estimate_code_uniqueness(preencoded.target_codes);
        if (uniqueness.sample_size > 0) {
            const double unique_ratio = static_cast<double>(uniqueness.unique_count) /
                                        static_cast<double>(uniqueness.sample_size);
            if (uniqueness.exact) {
                std::ostringstream oss;
                oss << "  Target code uniqueness: " << uniqueness.unique_count << "/"
                    << uniqueness.sample_size << " exact ("
                    << std::fixed << std::setprecision(2) << (unique_ratio * 100.0)
                    << "% unique)";
                progress::log(oss.str());
            } else {
                std::ostringstream oss;
                oss << "  Target code uniqueness estimate: " << uniqueness.unique_count << "/"
                    << uniqueness.sample_size << " sampled ("
                    << std::fixed << std::setprecision(2) << (unique_ratio * 100.0)
                    << "% unique)";
                progress::log(oss.str());
            }
        }
    }
    return preencoded;
}

std::vector<SearchHit> process_mismatch_search_loc_collect_preencoded(
    const std::vector<SpacerLocation>& query_locs,
    const std::vector<SpacerLocation>& target_locs,
    const MismatchPreencodedTargets& preencoded_targets,
    const GenomeReference& query_genome,
    const GenomeReference& target_genome,
    int max_mismatch,
    int num_threads
) {
    if (target_locs.empty() || query_locs.empty()) return {};

    bool fixed_len_fast_path = preencoded_targets.fixed_len_fast_path;
    const size_t fixed_target_len = preencoded_targets.fixed_target_len;
    const size_t fixed_query_len = query_locs.front().length;
    if (fixed_len_fast_path && fixed_query_len != fixed_target_len) {
        fixed_len_fast_path = false;
    }
    if (fixed_len_fast_path) {
        for (const auto& loc : query_locs) {
            if (loc.length != fixed_query_len) {
                fixed_len_fast_path = false;
                break;
            }
        }
    }

    const std::vector<encoded_seq_t>& target_codes = preencoded_targets.target_codes;
    const std::vector<int>& target_lengths = preencoded_targets.target_lengths;

    std::vector<std::string> query_seqs(query_locs.size());
    std::vector<std::string> query_names(query_locs.size());
    std::vector<encoded_seq_t> query_codes(query_locs.size());
    std::vector<int> query_lengths(query_locs.size(), 0);

    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(guided, 1)
    for (int i = 0; i < static_cast<int>(query_locs.size()); ++i) {
        const auto& loc = query_locs[static_cast<size_t>(i)];
        const std::string_view seq_view = query_genome.get_sequence(loc);
        query_names[static_cast<size_t>(i)] = query_genome.format_header(loc);
        if (seq_view.empty() || seq_view.length() > 32) continue;

        materialize_oriented_sequence(seq_view, loc.is_reverse, query_seqs[static_cast<size_t>(i)]);
        query_codes[static_cast<size_t>(i)] = encode_simple_oriented(seq_view, loc.is_reverse);
        query_lengths[static_cast<size_t>(i)] = static_cast<int>(seq_view.length());
    }

    std::vector<std::vector<SearchHit>> thread_hits(num_threads);
    const size_t query_batch_size = 16;

    if (fixed_len_fast_path) {
        for (size_t batch_start = 0; batch_start < query_locs.size(); batch_start += query_batch_size) {
            const size_t batch_count = std::min(query_batch_size, query_locs.size() - batch_start);
            const bool full_batch = (batch_count == query_batch_size);
            std::array<encoded_seq_t, query_batch_size> batch_query_codes{};
            std::array<uint8_t, query_batch_size> batch_query_valid{};

            for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
                const size_t query_idx = batch_start + local_idx;
                batch_query_codes[local_idx] = query_codes[query_idx];
                batch_query_valid[local_idx] = static_cast<uint8_t>(query_lengths[query_idx] != 0);
            }

            #pragma omp parallel
            {
                const int tid = omp_get_thread_num();
                auto& local_hits = thread_hits[tid];
                std::string target_seq_buffer;
                target_seq_buffer.reserve(fixed_target_len);
                std::string target_name_buffer;

                #pragma omp for schedule(static)
                for (int j = 0; j < static_cast<int>(target_locs.size()); ++j) {
                    const encoded_seq_t t_code = target_codes[static_cast<size_t>(j)];
                    bool target_materialized = false;
                    auto emit_hit = [&](size_t local_idx, int dist) {
                        if (!target_materialized) {
                            const auto& target_loc = target_locs[static_cast<size_t>(j)];
                            const std::string_view t_seq_view = target_genome.get_sequence(target_loc);
                            materialize_oriented_sequence(t_seq_view, target_loc.is_reverse, target_seq_buffer);
                            target_name_buffer = target_genome.format_header(target_loc);
                            target_materialized = true;
                        }

                        const size_t query_idx = batch_start + local_idx;
                        SearchHit hit;
                        hit.query_name = query_names[query_idx];
                        hit.query_seq = query_seqs[query_idx];
                        hit.target_name = target_name_buffer;
                        hit.target_seq = target_seq_buffer;
                        hit.distance = dist;
                        local_hits.push_back(std::move(hit));
                    };

                    if (full_batch) {
                        if (batch_query_valid[0]) {
                            const int dist = hamming_distance_simple(batch_query_codes[0], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(0, dist);
                        }
                        if (batch_query_valid[1]) {
                            const int dist = hamming_distance_simple(batch_query_codes[1], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(1, dist);
                        }
                        if (batch_query_valid[2]) {
                            const int dist = hamming_distance_simple(batch_query_codes[2], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(2, dist);
                        }
                        if (batch_query_valid[3]) {
                            const int dist = hamming_distance_simple(batch_query_codes[3], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(3, dist);
                        }
                        if (batch_query_valid[4]) {
                            const int dist = hamming_distance_simple(batch_query_codes[4], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(4, dist);
                        }
                        if (batch_query_valid[5]) {
                            const int dist = hamming_distance_simple(batch_query_codes[5], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(5, dist);
                        }
                        if (batch_query_valid[6]) {
                            const int dist = hamming_distance_simple(batch_query_codes[6], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(6, dist);
                        }
                        if (batch_query_valid[7]) {
                            const int dist = hamming_distance_simple(batch_query_codes[7], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(7, dist);
                        }
                        if (batch_query_valid[8]) {
                            const int dist = hamming_distance_simple(batch_query_codes[8], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(8, dist);
                        }
                        if (batch_query_valid[9]) {
                            const int dist = hamming_distance_simple(batch_query_codes[9], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(9, dist);
                        }
                        if (batch_query_valid[10]) {
                            const int dist = hamming_distance_simple(batch_query_codes[10], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(10, dist);
                        }
                        if (batch_query_valid[11]) {
                            const int dist = hamming_distance_simple(batch_query_codes[11], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(11, dist);
                        }
                        if (batch_query_valid[12]) {
                            const int dist = hamming_distance_simple(batch_query_codes[12], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(12, dist);
                        }
                        if (batch_query_valid[13]) {
                            const int dist = hamming_distance_simple(batch_query_codes[13], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(13, dist);
                        }
                        if (batch_query_valid[14]) {
                            const int dist = hamming_distance_simple(batch_query_codes[14], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(14, dist);
                        }
                        if (batch_query_valid[15]) {
                            const int dist = hamming_distance_simple(batch_query_codes[15], t_code, static_cast<int>(fixed_target_len));
                            if (dist <= max_mismatch) emit_hit(15, dist);
                        }
                        continue;
                    }

                    for (size_t local_idx = 0; local_idx < batch_count; ++local_idx) {
                        if (!batch_query_valid[local_idx]) continue;

                        const int dist = hamming_distance_simple(
                            batch_query_codes[local_idx], t_code, static_cast<int>(fixed_target_len));
                        if (dist > max_mismatch) continue;
                        emit_hit(local_idx, dist);
                    }
                }

            }
        }
    } else {
        omp_set_num_threads(num_threads);
        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            auto& local_hits = thread_hits[tid];
            std::string target_seq_buffer;
            std::string target_name_buffer;

            #pragma omp for schedule(dynamic, 64)
            for (int i = 0; i < static_cast<int>(query_locs.size()); ++i) {
                const int q_len = query_lengths[static_cast<size_t>(i)];
                if (q_len == 0) continue;

                for (size_t j = 0; j < target_locs.size(); ++j) {
                    if (target_lengths[j] != q_len) continue;
                    const int dist = hamming_distance_simple(query_codes[static_cast<size_t>(i)], target_codes[j], q_len);
                    if (dist > max_mismatch) continue;

                    const auto& target_loc = target_locs[j];
                    const std::string_view t_seq_view = target_genome.get_sequence(target_loc);
                    materialize_oriented_sequence(t_seq_view, target_loc.is_reverse, target_seq_buffer);
                    target_name_buffer = target_genome.format_header(target_loc);

                    SearchHit hit;
                    hit.query_name = query_names[static_cast<size_t>(i)];
                    hit.query_seq = query_seqs[static_cast<size_t>(i)];
                    hit.target_name = target_name_buffer;
                    hit.target_seq = target_seq_buffer;
                    hit.distance = dist;
                    local_hits.push_back(std::move(hit));
                }
            }
        }
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

void process_mismatch_search_loc(
    const std::vector<SpacerLocation>& query_locs,
    const std::vector<SpacerLocation>& target_locs,
    const GenomeReference& query_genome,
    const GenomeReference& target_genome,
    const std::string& output_filename,
    int max_mismatch,
    int num_threads
) {
    std::vector<SearchHit> collected_hits =
        process_mismatch_search_loc_collect(query_locs,
                                            target_locs,
                                            query_genome,
                                            target_genome,
                                            max_mismatch,
                                            num_threads);

    std::ofstream out_file(output_filename, std::ios::binary | std::ios::trunc);
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

void process_mismatch_search(const std::string& file_a,
                             const std::string& file_b,
                             const std::string& output_filename,
                             int max_mismatch,
                             int num_threads) {
    SequenceList sequences_b = parse_fasta(file_b);
    if (sequences_b.size() == 0) return;
    SequenceList sequences_a = parse_fasta(file_a);

    progress::log("[SIMPLE 2-BIT VERSION] Legacy mode - 2-bit brute-force (m=" +
                  std::to_string(max_mismatch) + ")");

    std::vector<encoded_seq_t> target_codes;
    std::vector<int> target_lengths;
    target_codes.reserve(sequences_b.size());
    target_lengths.reserve(sequences_b.size());
    for (const auto& seq : sequences_b.sequences) {
        target_codes.push_back(encode_simple(seq));
        target_lengths.push_back(static_cast<int>(seq.length()));
    }

    std::ofstream out_file(output_filename);
    out_file << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\n";

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 64)
    for (int i = 0; i < static_cast<int>(sequences_a.size()); ++i) {
        const encoded_seq_t q_code = encode_simple(sequences_a.sequences[static_cast<size_t>(i)]);
        const int q_len = static_cast<int>(sequences_a.sequences[static_cast<size_t>(i)].length());
        std::stringstream ss;

        for (size_t j = 0; j < sequences_b.size(); ++j) {
            if (target_lengths[j] != q_len) continue;
            const int dist = hamming_distance_simple(q_code, target_codes[j], q_len);
            if (dist <= max_mismatch) {
                ss << sequences_a.headers[static_cast<size_t>(i)] << "\t"
                   << sequences_a.sequences[static_cast<size_t>(i)] << "\t"
                   << sequences_b.headers[j] << "\t"
                   << sequences_b.sequences[j] << "\t"
                   << dist << "\n";
            }
        }

        if (ss.tellp() > 0) {
            #pragma omp critical
            out_file << ss.str();
        }
    }
}

MismatchIndexOptimized build_mismatch_index_optimized(const SequenceList& sequences_b, int max_mismatch, int num_threads) {
    (void)sequences_b;
    (void)max_mismatch;
    (void)num_threads;
    return MismatchIndexOptimized();
}
