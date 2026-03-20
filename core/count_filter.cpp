#include "count_filter.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>

// OPTIMIZATION: countFrequency now operates on uint64_t keys.
static bool countFrequency(const SequenceList& all_spacers,
                           std::unordered_map<uint64_t, int>& frequency_map) {
    frequency_map.clear();
    frequency_map.reserve(all_spacers.size() / 2);
    for (size_t i = 0; i < all_spacers.size(); ++i) {
        uint64_t key = encode_spacer(all_spacers.sequences[i]);
        if (key != std::numeric_limits<uint64_t>::max()) {
            frequency_map[key]++;
        }
    }
    return true;
}

// OPTIMIZATION: The main function now uses the uint64_t map for all filtering logic.
bool countFilter(const SequenceList& dusted_spacers,
                 SequenceList& candidate_spacers,
                 std::unordered_map<uint64_t, int>& frequency_map,
                 const std::string& over_represented_fasta_path,
                 int frequency_threshold) {

    countFrequency(dusted_spacers, frequency_map);

    std::ofstream out_file(over_represented_fasta_path);
    if (!out_file.is_open()) return false;

    std::unordered_set<uint64_t> written_over_rep;
    for (size_t i = 0; i < dusted_spacers.size(); ++i) {
        uint64_t key = encode_spacer(dusted_spacers.sequences[i]);
        if (key != std::numeric_limits<uint64_t>::max()) {
            auto it = frequency_map.find(key);
            if (it != frequency_map.end() && it->second > frequency_threshold) {
                if (written_over_rep.find(key) == written_over_rep.end()) {
                    out_file << ">count_" << it->second << "\n" << dusted_spacers.sequences[i] << "\n";
                    written_over_rep.insert(key);
                }
            }
        }
    }
    out_file.close();

    candidate_spacers.reserve(dusted_spacers.size());
    for (size_t i = 0; i < dusted_spacers.size(); ++i) {
        const std::string& seq_str = dusted_spacers.sequences[i];
        uint64_t key = encode_spacer(seq_str);
        if (key != std::numeric_limits<uint64_t>::max()) {
            auto it = frequency_map.find(key);
            if (it != frequency_map.end() && it->second <= frequency_threshold) {
                candidate_spacers.push_back(dusted_spacers.headers[i], seq_str);
            }
        }
    }
    return true;
}