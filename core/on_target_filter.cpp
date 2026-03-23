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
#include "on_target_filter.hpp"
#include <string>
#include <unordered_set>

// OPTIMIZATION: The filtering logic now uses the fast uint64_t key.
int onTargetFilter(const SequenceList& all_gene_spacers,
                   const std::unordered_map<uint64_t, int>& frequency_map,
                   SequenceList& final_on_target_spacers,
                   int frequency_threshold) {

    int discarded_count = 0;
    std::unordered_set<uint64_t> unique_spacers_added;
    final_on_target_spacers.reserve(all_gene_spacers.size());

    for (size_t i = 0; i < all_gene_spacers.size(); ++i) {
        const std::string& seq_str = all_gene_spacers.sequences[i];
        uint64_t key = encode_spacer(seq_str);
        if(key == std::numeric_limits<uint64_t>::max()) continue;

        int count = 0;
        auto it = frequency_map.find(key);
        if (it != frequency_map.end()) {
            count = it->second;
        }

        if (count > frequency_threshold) {
            discarded_count++;
        } else {
            if (unique_spacers_added.find(key) == unique_spacers_added.end()) {
                final_on_target_spacers.push_back(all_gene_spacers.headers[i], seq_str);
                unique_spacers_added.insert(key);
            }
        }
    }
    return discarded_count;
}