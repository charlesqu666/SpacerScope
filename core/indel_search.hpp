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
// indel_search.hpp
#pragma once

#include "types.hpp"
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cstdint>
#include <functional>
#include <vector>
#include <algorithm>
#include <array>
#include <memory>

class FlatHashSet64 {
private:
    std::vector<uint64_t> keys_;
    std::vector<uint8_t> occupied_;
    size_t mask_ = 0;
    size_t size_ = 0;

    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    static size_t next_power_of_two(size_t n) {
        size_t cap = 8;
        while (cap < n) cap <<= 1;
        return cap;
    }

    void rehash(size_t new_capacity) {
        std::vector<uint64_t> old_keys = std::move(keys_);
        std::vector<uint8_t> old_occupied = std::move(occupied_);

        keys_.assign(new_capacity, 0);
        occupied_.assign(new_capacity, 0);
        mask_ = new_capacity - 1;
        size_ = 0;

        for (size_t i = 0; i < old_keys.size(); ++i) {
            if (old_occupied[i]) insert(old_keys[i]);
        }
    }

public:
    FlatHashSet64() = default;

    explicit FlatHashSet64(size_t expected_size) {
        reserve(expected_size);
    }

    void reserve(size_t expected_size) {
        size_t needed = next_power_of_two((expected_size * 10 + 6) / 7);
        if (needed <= keys_.size()) return;
        rehash(needed);
    }

    bool insert(uint64_t key) {
        if (keys_.empty()) reserve(1024);
        if ((size_ + 1) * 10 > keys_.size() * 7) rehash(keys_.size() << 1);

        size_t idx = static_cast<size_t>(mix(key)) & mask_;
        while (occupied_[idx]) {
            if (keys_[idx] == key) return false;
            idx = (idx + 1) & mask_;
        }
        occupied_[idx] = 1;
        keys_[idx] = key;
        ++size_;
        return true;
    }

    bool contains(uint64_t key) const {
        if (keys_.empty()) return false;
        size_t idx = static_cast<size_t>(mix(key)) & mask_;
        while (occupied_[idx]) {
            if (keys_[idx] == key) return true;
            idx = (idx + 1) & mask_;
        }
        return false;
    }

    size_t size() const {
        return size_;
    }

    template <typename Fn>
    void for_each(Fn&& fn) const {
        for (size_t i = 0; i < keys_.size(); ++i) {
            if (occupied_[i]) fn(keys_[i]);
        }
    }
};

struct Channel {
    uint64_t a, c, g, t;
    size_t len;
    bool operator==(const Channel& other) const {
        return len == other.len && a == other.a && c == other.c && g == other.g && t == other.t;
    }
};

struct ChannelHash {
    size_t operator()(const Channel& ch) const {
        std::hash<uint64_t> hasher;
        size_t h = std::hash<size_t>()(ch.len);
        h ^= hasher(ch.a) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hasher(ch.c) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hasher(ch.g) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hasher(ch.t) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

class BinaryChannelFilter {
private:
    size_t fixed_len;
    uint64_t mask;
    FlatHashSet64 variants_A;
    FlatHashSet64 variants_C;
    FlatHashSet64 variants_G;
    FlatHashSet64 variants_T;
    int min_len, max_len;

    void add_normalized(const Channel& ch);
    void insert_variant(uint64_t a, uint64_t c, uint64_t g, uint64_t t);
    bool contains_variant(uint64_t a, uint64_t c, uint64_t g, uint64_t t) const;

    Channel apply_substitution(const Channel& ch, size_t pos, char new_base) const;
    Channel apply_deletion(const Channel& ch, size_t pos) const;
    Channel apply_insertion(const Channel& ch, size_t pos, char new_base) const;

public:
    BinaryChannelFilter(const std::string_view& query_seq, int max_edits, size_t spacer_len);
    bool check(const std::string_view& target_seq) const;
    void accumulate_query_mask_u32(uint16_t query_bit,
                                   std::vector<uint16_t>& a_map,
                                   std::vector<uint16_t>& c_map,
                                   std::vector<uint16_t>& g_map,
                                   std::vector<uint16_t>& t_map) const;
    
    // OPTIMIZATION: Check pre-encoded target without string conversion
    bool check_encoded(uint64_t t_a, uint64_t t_c, uint64_t t_g, uint64_t t_t, size_t target_len) const;
    bool check_encoded_fixed(uint64_t t_a, uint64_t t_c, uint64_t t_g, uint64_t t_t) const;
};

// Optimized right-anchored edit distance using banded DP
int right_anchored_distance(const std::string_view& Q, const std::string_view& T, int max_dist);

// Legacy interface: file-based
void process_indel_search(const std::string& file_a,
                          const std::string& file_b,
                          const std::string& output_filename,
                          int max_indel,
                          int num_threads,
                          int spacer_len);

// Pre-encoded targets in structure-of-arrays form for scan-heavy indel search.
struct PreencodedTargets {
    bool use_u32 = false;
    size_t count = 0;
    std::unique_ptr<uint32_t[]> a32;
    std::unique_ptr<uint32_t[]> c32;
    std::unique_ptr<uint32_t[]> g32;
    std::unique_ptr<uint32_t[]> t32;
    std::unique_ptr<uint64_t[]> a64;
    std::unique_ptr<uint64_t[]> c64;
    std::unique_ptr<uint64_t[]> g64;
    std::unique_ptr<uint64_t[]> t64;

    size_t size() const {
        return count;
    }
};

// OPTIMIZATION: New location-based interface (zero-copy, no temporary files)
void process_indel_search_loc(const std::vector<SpacerLocation>& query_locs,
                              const std::vector<SpacerLocation>& target_locs,
                              const GenomeReference& query_genome,
                              const GenomeReference& target_genome,
                              const std::string& output_filename,
                              int max_indel,
                              int num_threads,
                              int spacer_len);

std::vector<SearchHit> process_indel_search_loc_collect(const std::vector<SpacerLocation>& query_locs,
                                                        const std::vector<SpacerLocation>& target_locs,
                                                        const GenomeReference& query_genome,
                                                        const GenomeReference& target_genome,
                                                        int max_indel,
                                                        int num_threads,
                                                        int spacer_len);

std::vector<SearchHit> process_indel_search_loc_collect_preencoded(const std::vector<SpacerLocation>& query_locs,
                                                                   const std::vector<SpacerLocation>& target_locs,
                                                                   const GenomeReference& query_genome,
                                                                   const GenomeReference& target_genome,
                                                                   const PreencodedTargets& preencoded_targets,
                                                                   int max_indel,
                                                                   int num_threads,
                                                                   int spacer_len,
                                                                   size_t progress_base = 0,
                                                                   size_t progress_total = 0);

// OPTIMIZATION: Pre-encode targets for indel filter
PreencodedTargets preencode_targets_for_indel(
    const std::vector<SpacerLocation>& target_locs,
    const GenomeReference& target_genome,
    int num_threads);

// Check pre-encoded target against filter
bool check_preencoded(const BinaryChannelFilter& filter,
                      const PreencodedTargets& targets,
                      size_t idx);
