#pragma once

#include <string>
#include <string_view>

// Calculates Hamming distance, accelerated with SIMD.
int hamming_distance(std::string_view s1, std::string_view s2);

// Calculates Levenshtein distance with a bandwidth optimization.
int levenshtein_distance_banded(std::string_view s1, std::string_view s2, int max_indel);

int base_to_int(char c);