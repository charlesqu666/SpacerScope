#pragma once

#include "types.hpp"
#include <string>
#include <unordered_map>

// OPTIMIZATION: The frequency map now uses uint64_t as its key for major performance gain.
bool countFilter(const SequenceList& dusted_spacers,
                 SequenceList& candidate_spacers,
                 std::unordered_map<uint64_t, int>& frequency_map, // Changed key type
                 const std::string& over_represented_fasta_path,
                 int frequency_threshold = 10);