#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <string_view>
#include <tuple>
#include <cstdint>
#include <limits>
#include <cctype>
#include <array>
#include <set>
#include <stdexcept>

// Forward declare to resolve circular dependency
struct SequenceList;
#include "utils.hpp"

// --- Core Type Aliases ---
using Sequence = std::string;
using SequenceID = size_t;
using MatchResult = std::tuple<std::string, std::string, int>;

// --- Correct SequenceList struct, with helper methods ---
struct SequenceList {
    std::vector<std::string> headers;
    std::vector<std::string> sequences;
    void push_back(const std::string& header, const std::string& seq) {
        headers.push_back(header);
        sequences.push_back(seq);
    }

    void reserve(size_t num_sequences) {
        headers.reserve(num_sequences);
        sequences.reserve(num_sequences);
    }

    std::string_view get(size_t idx) const {
        if (idx >= sequences.size()) return "";
        return sequences[idx];
    }

    std::string_view get_header(size_t idx) const {
        if (idx >= headers.size()) return "";
        return headers[idx];
    }

    size_t size() const { return sequences.size(); }
};

// --- OPTIMIZATION: Location-based spacer reference (compact index form) ---
// Stores only chromosome index and coordinates, keeping per-hit memory low.
struct SpacerLocation {
    uint32_t chrom_idx;         // Chromosome/scaffold index in GenomeReference
    uint32_t start_pos;         // 0-based start position
    uint32_t length;            // Spacer length
    bool is_reverse;            // true = reverse strand
};

struct SearchHit {
    std::string query_name;
    std::string query_seq;
    std::string target_name;
    std::string target_seq;
    int distance = 0;
    std::string search_type;
    int genomic_frequency = 0;
};

static_assert(sizeof(SpacerLocation) <= 16, "SpacerLocation should remain compact");

// Genome reference data - keeps original sequences in memory for fast lookup
struct GenomeReference {
    std::vector<std::string> chrom_names;
    std::vector<std::string> chrom_sequences;
    std::unordered_map<std::string, size_t> chrom_name_to_idx;
    
    void add_chromosome(const std::string& name, std::string&& seq) {
        if (chrom_names.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::runtime_error("Too many chromosomes for compact SpacerLocation");
        }
        if (seq.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::runtime_error("Chromosome length exceeds compact SpacerLocation limit");
        }
        chrom_name_to_idx[name] = chrom_names.size();
        chrom_names.push_back(name);
        chrom_sequences.push_back(std::move(seq));
    }
    
    size_t get_chrom_idx(const std::string& name) const {
        auto it = chrom_name_to_idx.find(name);
        return (it != chrom_name_to_idx.end()) ? it->second : static_cast<size_t>(-1);
    }
    
    // Get sequence view by location (zero-copy)
    std::string_view get_sequence(size_t chrom_idx, size_t start, size_t len) const {
        if (chrom_idx >= chrom_sequences.size()) return "";
        const auto& seq = chrom_sequences[chrom_idx];
        if (start + len > seq.length()) return "";
        return std::string_view(seq.data() + start, len);
    }
    
    // Get sequence by location
    std::string_view get_sequence(const SpacerLocation& loc) const {
        return get_sequence(static_cast<size_t>(loc.chrom_idx),
                            static_cast<size_t>(loc.start_pos),
                            static_cast<size_t>(loc.length));
    }

    std::string format_header(const SpacerLocation& loc) const {
        const size_t idx = static_cast<size_t>(loc.chrom_idx);
        if (idx >= chrom_names.size()) return "";
        return chrom_names[idx] + ":" + std::to_string(static_cast<size_t>(loc.start_pos) + 1) + "-" +
               std::to_string(static_cast<size_t>(loc.start_pos) + static_cast<size_t>(loc.length)) +
               (loc.is_reverse ? ":-" : ":+");
    }
    
    size_t size() const { return chrom_names.size(); }
};

// --- Indexing and Caching Types ---

// ** OLD MismatchIndex (Pigeonhole/Q-gram Index) **
using MismatchPartition = std::string;
using MismatchIndexMap = std::unordered_map<MismatchPartition, std::vector<SequenceID>>;
// 索引是一个向量，每个元素代表一个分区（q-gram）的哈希表
using MismatchIndex = std::vector<MismatchIndexMap>; 


// ** NEW MismatchIndex (Binary Channel based) **
// Key is the 64-bit integer channel value, Value is the list of Sequence IDs
// using ChannelMap = std::unordered_map<uint64_t, std::vector<SequenceID>>; // 不再需要
// using MismatchIndex = std::array<ChannelMap, 4>; // 不再用于索引，只用于距离计算

using QGramIndex = std::vector<std::vector<SequenceID>>;

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

struct DistanceCache {
    size_t num_locks;
    std::vector<std::mutex> mutexes;
    std::vector<std::unordered_map<std::pair<std::string, std::string>, int, pair_hash>> caches;

    explicit DistanceCache(size_t locks = 128) : num_locks(locks), mutexes(locks), caches(locks) {}

    int get(std::string_view s1, std::string_view s2, int max_indel) {
        auto key = std::make_pair(std::string(s1), std::string(s2));
        size_t h = pair_hash{}(key) % num_locks;
        {
            std::lock_guard<std::mutex> lock(mutexes[h]);
            auto it = caches[h].find(key);
            if (it != caches[h].end()) {
                return it->second;
            }
        }
        // NOTE: assumes levenshtein_distance_banded is defined in utils.hpp/cpp
        int dist = levenshtein_distance_banded(s1, s2, max_indel); 
        {
            std::lock_guard<std::mutex> lock(mutexes[h]);
            caches[h][key] = dist;
        }
        return dist;
    }
};

// OPTIMIZATION: Encodes a DNA sequence into a 64-bit integer for fast hashing and comparison.
inline uint64_t encode_spacer(std::string_view seq) {
    uint64_t val = 0;
    if (seq.length() > 32) return std::numeric_limits<uint64_t>::max(); // Safety check

    for (char c : seq) {
        switch (std::toupper(c)) {
            case 'A': val = (val << 2) | 0; break;
            case 'C': val = (val << 2) | 1; break;
            case 'G': val = (val << 2) | 2; break;
            case 'T': val = (val << 2) | 3; break;
            default: return std::numeric_limits<uint64_t>::max(); // Invalid character
        }
    }
    return val;
}
