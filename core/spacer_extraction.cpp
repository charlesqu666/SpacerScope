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
#include "spacer_extraction.hpp"
#include "fasta_parser.hpp"
#include <iostream>
#include <fstream>
#include <string_view>
#include <vector>
#include <cctype>
#include <algorithm>
#include <thread>
#include <mutex>
#include <future>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

// Global mutex for thread-safe extraction
std::mutex g_extraction_mutex;

// --- Helper functions ---

static char complement_base(char c) {
    switch (std::toupper(c)) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        default:  return 'N';
    }
}

static std::string reverse_complement(std::string_view seq) {
    std::string rc_seq;
    rc_seq.reserve(seq.length());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc_seq += complement_base(*it);
    }
    return rc_seq;
}

static std::string reverse_complement_pam(const std::string& pam) {
    std::string rc_pam;
    rc_pam.reserve(pam.length());
    for (auto it = pam.rbegin(); it != pam.rend(); ++it) {
        rc_pam += complement_base(*it);
    }
    return rc_pam;
}

static bool pam_matches(std::string_view seq_chunk, const std::string& pam) {
    if (seq_chunk.length() != pam.length()) return false;
    for (size_t i = 0; i < pam.length(); ++i) {
        if (std::toupper(pam[i]) != 'N' && std::toupper(pam[i]) != std::toupper(seq_chunk[i])) {
            return false;
        }
    }
    return true;
}

static bool pam_matches_ptr(const char* seq_chunk, const std::string& pam) {
    for (size_t i = 0; i < pam.length(); ++i) {
        const char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pam[i])));
        const char s = static_cast<char>(std::toupper(static_cast<unsigned char>(seq_chunk[i])));
        if (p != 'N' && p != s) {
            return false;
        }
    }
    return true;
}

static bool is_ngg_nag_pair(const std::string& pam1, const std::string& pam2) {
    auto upper = [](const std::string& s) {
        std::string out = s;
        std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        return out;
    };
    const std::string a = upper(pam1);
    const std::string b = upper(pam2);
    return (a == "NGG" && b == "NAG") || (a == "NAG" && b == "NGG");
}

static inline bool match_ngg_fast(const char* s) {
    return s[1] == 'G' && s[2] == 'G';
}

static inline bool match_nag_fast(const char* s) {
    return s[1] == 'A' && s[2] == 'G';
}

static inline bool match_ccn_fast(const char* s) {
    return s[0] == 'C' && s[1] == 'C';
}

static inline bool match_ctn_fast(const char* s) {
    return s[0] == 'C' && s[1] == 'T';
}

static inline int base_to_idx(char c) {
    switch (std::toupper(c)) {
        case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3;
        default: return -1;
    }
}

static bool is_complex(std::string_view seq, int threshold) {
    if (seq.length() < 3) return true;
    uint16_t counts[64] = {0};
    for (size_t j = 0; j <= seq.length() - 3; ++j) {
        int b1 = base_to_idx(seq[j]), b2 = base_to_idx(seq[j+1]), b3 = base_to_idx(seq[j+2]);
        if (b1 != -1 && b2 != -1 && b3 != -1) counts[(b1 << 4) | (b2 << 2) | b3]++;
    }
    int score = 0;
    for (int c : counts) if (c > 1) score += c * (c - 1) / 2;
    return score <= threshold;
}

// --- Multi-threading processors ---

static void process_segment_with_filter(
    const std::string& seq, const std::string& name, const std::string& pam, const std::string& rc_pam,
    int len, bool up, size_t start, size_t end, int dust_thru, SequenceList& out, size_t offset) 
{
    SequenceList local_list;
    size_t p_len = pam.length();

    for (size_t j = start; j < end && (j + p_len <= seq.length()); ++j) {
        // ������
        if (pam_matches(std::string_view(seq).substr(j, p_len), pam)) {
            size_t pos; bool ok = false;
            if (up) { if (j >= (size_t)len) { pos = j - len; ok = true; } }
            else { if (j + p_len + len <= seq.length()) { pos = j + p_len; ok = true; } }
            
            if (ok) {
                std::string_view s_view = std::string_view(seq).substr(pos, len);
                if (is_complex(s_view, dust_thru)) {
                    std::string header = name + ":" + std::to_string(pos + offset + 1) + "-" + 
                                        std::to_string(pos + offset + len) + ":+";
                    local_list.push_back(header, std::string(s_view));
                }
            }
        }
        // ������
        if (pam_matches(std::string_view(seq).substr(j, p_len), rc_pam)) {
            size_t pos; bool ok = false;
            if (up) { if (j + p_len + len <= seq.length()) { pos = j + p_len; ok = true; } }
            else { if (j >= (size_t)len) { pos = j - len; ok = true; } }

            if (ok) {
                std::string_view s_view = std::string_view(seq).substr(pos, len);
                if (is_complex(s_view, dust_thru)) {
                    std::string header = name + ":" + std::to_string(pos + offset + 1) + "-" + 
                                        std::to_string(pos + offset + len) + ":-";
                    local_list.push_back(header, reverse_complement(s_view));
                }
            }
        }
    }

    if (local_list.size() > 0) {
        std::lock_guard<std::mutex> lock(g_extraction_mutex);
        for (size_t i = 0; i < local_list.size(); ++i) 
            out.push_back(local_list.headers[i], local_list.sequences[i]);
    }
}

// --- Location-based processor (NEW: zero-copy extraction) ---

static void process_segment_loc_with_filter(
    const std::string& seq, uint32_t chrom_idx, const std::string& pam, const std::string& rc_pam,
    int len, bool up, size_t start, size_t end, int dust_thru, std::vector<SpacerLocation>& out, size_t offset)
{
    std::vector<SpacerLocation> local_list;
    size_t p_len = pam.length();

    for (size_t j = start; j < end && (j + p_len <= seq.length()); ++j) {
        // Forward strand (+)
        if (pam_matches(std::string_view(seq).substr(j, p_len), pam)) {
            size_t pos; bool ok = false;
            if (up) { if (j >= (size_t)len) { pos = j - len; ok = true; } }
            else { if (j + p_len + len <= seq.length()) { pos = j + p_len; ok = true; } }
            
            if (ok) {
                std::string_view s_view = std::string_view(seq).substr(pos, len);
                if (is_complex(s_view, dust_thru)) {
                    local_list.push_back({chrom_idx,
                                          static_cast<uint32_t>(pos + offset),
                                          static_cast<uint32_t>(len),
                                          false});
                }
            }
        }
        // Reverse strand (-)
        if (pam_matches(std::string_view(seq).substr(j, p_len), rc_pam)) {
            size_t pos; bool ok = false;
            if (up) { if (j + p_len + len <= seq.length()) { pos = j + p_len; ok = true; } }
            else { if (j >= (size_t)len) { pos = j - len; ok = true; } }

            if (ok) {
                std::string_view s_view = std::string_view(seq).substr(pos, len);
                if (is_complex(s_view, dust_thru)) {
                    local_list.push_back({chrom_idx,
                                          static_cast<uint32_t>(pos + offset),
                                          static_cast<uint32_t>(len),
                                          true});
                }
            }
        }
    }

    if (!local_list.empty()) {
        std::lock_guard<std::mutex> lock(g_extraction_mutex);
        out.insert(out.end(), local_list.begin(), local_list.end());
    }
}

static void process_segment_loc_with_dual_filter(
    const std::string& seq,
    uint32_t chrom_idx,
    const std::string& pam,
    const std::string& alt_pam,
    const std::string& rc_pam,
    const std::string& rc_alt_pam,
    int len,
    bool up,
    size_t start,
    size_t end,
    int dust_thru,
    std::vector<SpacerLocation>& out,
    size_t offset,
    bool use_ngg_nag_fast_path)
{
    std::vector<SpacerLocation> local_list;
    const size_t p_len = pam.length();
    const char* data = seq.data();
    const size_t scan_end = std::min(end, seq.length() - p_len + 1);

    auto add_forward = [&](size_t j) {
        size_t pos;
        bool ok = false;
        if (up) {
            if (j >= static_cast<size_t>(len)) {
                pos = j - len;
                ok = true;
            }
        } else {
            if (j + p_len + len <= seq.length()) {
                pos = j + p_len;
                ok = true;
            }
        }
        if (ok && is_complex(std::string_view(seq).substr(pos, len), dust_thru)) {
            local_list.push_back({chrom_idx,
                                  static_cast<uint32_t>(pos + offset),
                                  static_cast<uint32_t>(len),
                                  false});
        }
    };

    auto add_reverse = [&](size_t j) {
        size_t pos;
        bool ok = false;
        if (up) {
            if (j + p_len + len <= seq.length()) {
                pos = j + p_len;
                ok = true;
            }
        } else {
            if (j >= static_cast<size_t>(len)) {
                pos = j - len;
                ok = true;
            }
        }
        if (ok && is_complex(std::string_view(seq).substr(pos, len), dust_thru)) {
            local_list.push_back({chrom_idx,
                                  static_cast<uint32_t>(pos + offset),
                                  static_cast<uint32_t>(len),
                                  true});
        }
    };

    if (use_ngg_nag_fast_path) {
        for (size_t j = start; j < scan_end; ++j) {
            const char* s = data + j;
            if (match_ngg_fast(s)) {
                add_forward(j);
            }
            if (match_nag_fast(s)) {
                add_forward(j);
            }
            if (match_ccn_fast(s)) {
                add_reverse(j);
            }
            if (match_ctn_fast(s)) {
                add_reverse(j);
            }
        }
    } else {
        const bool same_forward = (pam == alt_pam);
        const bool same_reverse = (rc_pam == rc_alt_pam);
        for (size_t j = start; j < scan_end; ++j) {
            const char* s = data + j;
            if (pam_matches_ptr(s, pam)) {
                add_forward(j);
            }
            if (!same_forward && pam_matches_ptr(s, alt_pam)) {
                add_forward(j);
            }
            if (pam_matches_ptr(s, rc_pam)) {
                add_reverse(j);
            }
            if (!same_reverse && pam_matches_ptr(s, rc_alt_pam)) {
                add_reverse(j);
            }
        }
    }

    if (!local_list.empty()) {
        out.insert(out.end(), local_list.begin(), local_list.end());
    }
}

// --- Interface 1: File-based extraction (legacy) ---

bool spacerExtraction(const std::string& genome_file, const std::string& pam, int spacer_len, 
                      bool upstream, SequenceList& spacers, int num_threads, int score_threshold) {
    std::ifstream fasta(genome_file);
    if (!fasta.is_open()) return false;

    std::string rc_pam;
    rc_pam.reserve(pam.length());
    for (auto it = pam.rbegin(); it != pam.rend(); ++it) rc_pam += complement_base(*it);

    const size_t CHUNK_SIZE = 10 * 1024 * 1024; // 10MB
    const size_t OVERLAP = spacer_len + pam.length() + 5; 
    
    std::string line, cur_name, buffer;
    std::vector<std::future<void>> futures;
    size_t current_offset = 0;

    auto submit_chunk = [&](std::string seq_chunk, std::string name, size_t offset) {
        auto shared_seq = std::make_shared<std::string>(std::move(seq_chunk));
        size_t sub_chunk_size = (shared_seq->length() + num_threads - 1) / num_threads;
        
        for (int t = 0; t < num_threads; ++t) {
            size_t start = t * sub_chunk_size;
            if (start >= shared_seq->length()) break;
            size_t end = std::min(start + sub_chunk_size, shared_seq->length());

            futures.push_back(std::async(std::launch::async, [=, &spacers]() {
                process_segment_with_filter(*shared_seq, name, pam, rc_pam, spacer_len, 
                                          upstream, start, end, score_threshold, spacers, offset);
            }));
        }
        while (futures.size() > (size_t)num_threads * 4) {
            futures.begin()->get();
            futures.erase(futures.begin());
        }
    };

    while (std::getline(fasta, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!buffer.empty()) {
                submit_chunk(buffer, cur_name, current_offset);
                buffer.clear();
            }
            cur_name = line.substr(1);
            size_t first_space = cur_name.find_first_of(" \t\r");
            if (first_space != std::string::npos) cur_name = cur_name.substr(0, first_space);
            current_offset = 0;
        } else {
            buffer += line;
            if (buffer.length() >= CHUNK_SIZE + OVERLAP) {
                std::string to_process = buffer.substr(0, CHUNK_SIZE + OVERLAP);
                submit_chunk(to_process, cur_name, current_offset);
                buffer.erase(0, CHUNK_SIZE);
                current_offset += CHUNK_SIZE;
            }
        }
    }
    if (!buffer.empty()) submit_chunk(buffer, cur_name, current_offset);
    for (auto& f : futures) f.get();
    return true;
}

// --- Interface 2: In-memory extraction (legacy) ---

bool spacerExtraction(const SequenceList& sequences, const std::string& pam, int spacer_len, 
                      bool upstream, SequenceList& spacers, int num_threads, int score_threshold) {
    
    std::string rc_pam;
    rc_pam.reserve(pam.length());
    for (auto it = pam.rbegin(); it != pam.rend(); ++it) rc_pam += complement_base(*it);

    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < sequences.size(); ++i) {
        const std::string& seq = sequences.sequences[i];
        const std::string& name = sequences.headers[i];
        size_t total_len = seq.length();
        if (total_len < (size_t)pam.length()) continue;

        size_t sub_chunk = (total_len + num_threads - 1) / num_threads;
        for (int t = 0; t < num_threads; ++t) {
            size_t start = t * sub_chunk;
            if (start >= total_len) break;
            size_t end = std::min(start + sub_chunk, total_len);
            const std::string* seq_ptr = &seq;
            const std::string* name_ptr = &name;

            futures.push_back(std::async(std::launch::async, [&, seq_ptr, name_ptr, start, end]() {
                process_segment_with_filter(*seq_ptr, *name_ptr, pam, rc_pam, spacer_len, upstream, 
                                          start, end, score_threshold, spacers, 0);
            }));
        }
    }
    for (auto& f : futures) f.get();
    return true;
}
// --- NEW: Location-based extraction interfaces (zero-copy) ---

// Interface 3: File-based extraction returning locations
bool spacerExtractionLoc(const std::string& genome_file, const std::string& pam, int spacer_len,
                         bool upstream, std::vector<SpacerLocation>& locations, int num_threads, int score_threshold) {
    std::ifstream fasta(genome_file);
    if (!fasta.is_open()) return false;

    std::string rc_pam;
    rc_pam.reserve(pam.length());
    for (auto it = pam.rbegin(); it != pam.rend(); ++it) rc_pam += complement_base(*it);

    const size_t CHUNK_SIZE = 10 * 1024 * 1024; // 10MB
    const size_t OVERLAP = spacer_len + pam.length() + 5;
    
    std::string line, cur_name, buffer;
    std::vector<std::future<void>> futures;
    size_t current_offset = 0;
    uint32_t current_chrom_idx = 0;

    auto submit_chunk = [&](std::string seq_chunk, uint32_t chrom_idx, size_t offset) {
        auto shared_seq = std::make_shared<std::string>(std::move(seq_chunk));
        size_t sub_chunk_size = (shared_seq->length() + num_threads - 1) / num_threads;
        
        for (int t = 0; t < num_threads; ++t) {
            size_t start = t * sub_chunk_size;
            if (start >= shared_seq->length()) break;
            size_t end = std::min(start + sub_chunk_size, shared_seq->length());

            futures.push_back(std::async(std::launch::async, [=, &locations]() {
                process_segment_loc_with_filter(*shared_seq, chrom_idx, pam, rc_pam, spacer_len,
                                                upstream, start, end, score_threshold, locations, offset);
            }));
        }
        while (futures.size() > (size_t)num_threads * 4) {
            futures.begin()->get();
            futures.erase(futures.begin());
        }
    };

    while (std::getline(fasta, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!buffer.empty()) {
                submit_chunk(buffer, current_chrom_idx, current_offset);
                buffer.clear();
                ++current_chrom_idx;
            }
            cur_name = line.substr(1);
            size_t first_space = cur_name.find_first_of(" \t\r");
            if (first_space != std::string::npos) cur_name = cur_name.substr(0, first_space);
            current_offset = 0;
        } else {
            buffer += line;
            if (buffer.length() >= CHUNK_SIZE + OVERLAP) {
                std::string to_process = buffer.substr(0, CHUNK_SIZE + OVERLAP);
                submit_chunk(to_process, current_chrom_idx, current_offset);
                buffer.erase(0, CHUNK_SIZE);
                current_offset += CHUNK_SIZE;
            }
        }
    }
    if (!buffer.empty()) submit_chunk(buffer, current_chrom_idx, current_offset);
    for (auto& f : futures) f.get();
    return true;
}

// Interface 4: In-memory extraction returning locations
bool spacerExtractionLoc(const SequenceList& sequences, const std::string& pam, int spacer_len,
                         bool upstream, std::vector<SpacerLocation>& locations, int num_threads, int score_threshold) {
    
    std::string rc_pam;
    rc_pam.reserve(pam.length());
    for (auto it = pam.rbegin(); it != pam.rend(); ++it) rc_pam += complement_base(*it);

    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < sequences.size(); ++i) {
        const std::string& seq = sequences.sequences[i];
        size_t total_len = seq.length();
        if (total_len < (size_t)pam.length()) continue;

        size_t sub_chunk = (total_len + num_threads - 1) / num_threads;
        for (int t = 0; t < num_threads; ++t) {
            size_t start = t * sub_chunk;
            if (start >= total_len) break;
            size_t end = std::min(start + sub_chunk, total_len);
            const std::string* seq_ptr = &seq;
            const uint32_t chrom_idx = static_cast<uint32_t>(i);

            futures.push_back(std::async(std::launch::async, [&, seq_ptr, chrom_idx, start, end]() {
                process_segment_loc_with_filter(*seq_ptr, chrom_idx, pam, rc_pam, spacer_len, upstream,
                                                start, end, score_threshold, locations, 0);
            }));
        }
    }
    for (auto& f : futures) f.get();
    return true;
}

// Interface 5: GenomeReference-based extraction returning locations
bool spacerExtractionLoc(const GenomeReference& genome, const std::string& pam, int spacer_len,
                         bool upstream, std::vector<SpacerLocation>& locations, int num_threads, int score_threshold) {
    
    std::string rc_pam = reverse_complement_pam(pam);

    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < genome.size(); ++i) {
        const std::string& seq = genome.chrom_sequences[i];
        size_t total_len = seq.length();
        if (total_len < (size_t)pam.length()) continue;

        size_t sub_chunk = (total_len + num_threads - 1) / num_threads;
        for (int t = 0; t < num_threads; ++t) {
            size_t start = t * sub_chunk;
            if (start >= total_len) break;
            size_t end = std::min(start + sub_chunk, total_len);
            const std::string* seq_ptr = &seq;
            const uint32_t chrom_idx = static_cast<uint32_t>(i);

            futures.push_back(std::async(std::launch::async, [&, seq_ptr, chrom_idx, start, end]() {
                process_segment_loc_with_filter(*seq_ptr, chrom_idx, pam, rc_pam, spacer_len, upstream,
                                                start, end, score_threshold, locations, 0);
            }));
        }
    }
    for (auto& f : futures) f.get();
    return true;
}

bool spacerExtractionLoc(const GenomeReference& genome, const std::string& pam, const std::string& alt_pam,
                         int spacer_len, bool upstream, std::vector<SpacerLocation>& locations,
                         int num_threads, int score_threshold) {
    if (alt_pam.empty() || pam.length() != alt_pam.length()) {
        const bool ok_primary = spacerExtractionLoc(
            genome, pam, spacer_len, upstream, locations, num_threads, score_threshold);
        if (!ok_primary || alt_pam.empty()) {
            return ok_primary;
        }
        return spacerExtractionLoc(
            genome, alt_pam, spacer_len, upstream, locations, num_threads, score_threshold);
    }

    const std::string rc_pam = reverse_complement_pam(pam);
    const std::string rc_alt_pam = reverse_complement_pam(alt_pam);
    const bool use_ngg_nag_fast_path = is_ngg_nag_pair(pam, alt_pam);
    const int threads = std::max(1, num_threads);
    const size_t p_len = pam.length();

    for (size_t i = 0; i < genome.size(); ++i) {
        const std::string& seq = genome.chrom_sequences[i];
        const size_t total_len = seq.length();
        if (total_len < p_len || total_len < static_cast<size_t>(spacer_len)) {
            continue;
        }

        const size_t scan_positions = total_len - p_len + 1;
        const size_t worker_count = std::min(scan_positions, static_cast<size_t>(threads));
        std::vector<std::vector<SpacerLocation>> thread_lists(worker_count);

#ifdef _OPENMP
#pragma omp parallel num_threads(static_cast<int>(worker_count))
        {
            const int tid = omp_get_thread_num();
            const int team_size = omp_get_num_threads();
            const size_t chunk = (scan_positions + static_cast<size_t>(team_size) - 1) /
                                 static_cast<size_t>(team_size);
            const size_t start = static_cast<size_t>(tid) * chunk;
            const size_t end = std::min(start + chunk, scan_positions);
            if (start < end) {
                process_segment_loc_with_dual_filter(
                    seq,
                    static_cast<uint32_t>(i),
                    pam,
                    alt_pam,
                    rc_pam,
                    rc_alt_pam,
                    spacer_len,
                    upstream,
                    start,
                    end,
                    score_threshold,
                    thread_lists[static_cast<size_t>(tid)],
                    0,
                    use_ngg_nag_fast_path);
            }
        }
#else
        std::vector<std::future<void>> futures;
        futures.reserve(worker_count);
        const size_t chunk = (scan_positions + worker_count - 1) / worker_count;
        for (size_t t = 0; t < worker_count; ++t) {
            const size_t start = t * chunk;
            const size_t end = std::min(start + chunk, scan_positions);
            if (start >= end) {
                break;
            }
            futures.push_back(std::async(std::launch::async, [&, t, start, end]() {
                process_segment_loc_with_dual_filter(
                    seq,
                    static_cast<uint32_t>(i),
                    pam,
                    alt_pam,
                    rc_pam,
                    rc_alt_pam,
                    spacer_len,
                    upstream,
                    start,
                    end,
                    score_threshold,
                    thread_lists[t],
                    0,
                    use_ngg_nag_fast_path);
            }));
        }
        for (auto& f : futures) {
            f.get();
        }
#endif

        for (auto& local : thread_lists) {
            locations.insert(locations.end(), local.begin(), local.end());
        }
    }
    return true;
}
