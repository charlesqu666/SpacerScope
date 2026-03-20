#pragma once

#include "types.hpp"
#include <string>
#include <unordered_map>

// OPTIMIZATION: The frequency map now uses uint64_t as its key.
int onTargetFilter(const SequenceList& all_gene_spacers,
                   const std::unordered_map<uint64_t, int>& frequency_map, // Changed key type
                   SequenceList& final_on_target_spacers,
                   int frequency_threshold = 10);