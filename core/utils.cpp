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
#include "utils.hpp"
#include "types.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <cctype>
#include <string>

// --- A NEW, ROBUST IMPLEMENTATION to replace the fragile, transplanted code ---

/**
 * @brief Calculates the edit distance based on "events" (affine gap style).
 * Substitutions are counted one by one.
 * Consecutive blocks of indels (e.g., "---" or "AAA") are counted as a single event.
 * This function correctly implements your intended scoring logic in a robust way.
 *
 * @param s1 First sequence.
 * @param s2 Second sequence.
 * @param max_dist An optimization to stop early if the distance exceeds a threshold.
 * @return The calculated event-based distance.
 */
int calculate_distance_with_affine_gap(std::string_view s1, std::string_view s2, int max_dist) {
    int len1 = s1.length();
    int len2 = s2.length();

    // The DP table stores the standard Levenshtein distance
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));

    for (int i = 0; i <= len1; ++i) dp[i][0] = i;
    for (int j = 0; j <= len2; ++j) dp[0][j] = j;

    for (int i = 1; i <= len1; ++i) {
        for (int j = 1; j <= len2; ++j) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            dp[i][j] = std::min({dp[i - 1][j] + 1,        // Deletion
                                 dp[i][j - 1] + 1,        // Insertion
                                 dp[i - 1][j - 1] + cost}); // Match/Mismatch
        }
    }

    // If standard Levenshtein distance already exceeds the max, we can stop.
    if (dp[len1][len2] > max_dist) {
        return dp[len1][len2];
    }

    // Now, backtrack through the DP matrix to count "events" correctly.
    int events = 0;
    int i = len1, j = len2;
    bool in_indel = false;

    while (i > 0 || j > 0) {
        int current_dist = dp[i][j];
        int sub_cost = (i > 0 && j > 0 && s1[i - 1] != s2[j - 1]) ? 1 : 0;
        
        int match_diag = (i > 0 && j > 0) ? dp[i - 1][j - 1] + sub_cost : -1;
        int del_up = (i > 0) ? dp[i - 1][j] + 1 : -1;
        int ins_left = (j > 0) ? dp[i][j - 1] + 1 : -1;

        if (match_diag != -1 && current_dist == match_diag) {
            if (sub_cost == 1) {
                events++;
            }
            i--; j--;
            in_indel = false;
        } else if (del_up != -1 && current_dist == del_up) {
            if (!in_indel) {
                events++;
                in_indel = true;
            }
            i--;
        } else if (ins_left != -1 && current_dist == ins_left) {
            if (!in_indel) {
                events++;
                in_indel = true;
            }
            j--;
        } else {
             // Should not happen if DP is correct, but as a safeguard
             if (i > 0) i--;
             if (j > 0) j--;
        }
    }
    
    return events;
}


// --- ADAPTER FUNCTION: We now call the new robust implementation ---
int levenshtein_distance_banded(std::string_view s1, std::string_view s2, int max_indel) {
    // The old SW implementation is no longer used.
    // We now call the new, robust function that correctly implements your logic.
    return calculate_distance_with_affine_gap(s1, s2, max_indel);
}


// --- Other utility functions remain the same ---
int hamming_distance(std::string_view s1, std::string_view s2) {
    if (s1.length() != s2.length()) {
        return -1;
    }
    int diff = 0;
    for (size_t i = 0; i < s1.length(); ++i) {
        if (s1[i] != s2[i]) {
            diff++;
        }
    }
    return diff;
}

int base_to_int(char c) {
    switch (std::toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}