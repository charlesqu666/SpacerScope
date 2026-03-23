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
#pragma once

#include "types.hpp"
#include <string>
#include <vector>

// Original interfaces (kept for compatibility)
bool spacerExtraction(const SequenceList& sequences,
                      const std::string& pam,
                      int spacer_len,
                      bool upstream,
                      SequenceList& spacers,
                      int num_threads = 1,
                      int score_threshold = 12);

bool spacerExtraction(const std::string& genome_file,
                      const std::string& pam,
                      int spacer_len,
                      bool upstream,
                      SequenceList& spacers,
                      int num_threads = 1,
                      int score_threshold = 12);

// OPTIMIZATION: New location-based interfaces (zero-copy extraction)
// Extract spacers from genome file, returning only location info
bool spacerExtractionLoc(const std::string& genome_file,
                         const std::string& pam,
                         int spacer_len,
                         bool upstream,
                         std::vector<SpacerLocation>& locations,
                         int num_threads = 1,
                         int score_threshold = 12);

// Extract spacers from in-memory sequences, returning location info
bool spacerExtractionLoc(const SequenceList& sequences,
                         const std::string& pam,
                         int spacer_len,
                         bool upstream,
                         std::vector<SpacerLocation>& locations,
                         int num_threads = 1,
                         int score_threshold = 12);

// Extract spacers from genome reference (used after loading genome once)
bool spacerExtractionLoc(const GenomeReference& genome,
                         const std::string& pam,
                         int spacer_len,
                         bool upstream,
                         std::vector<SpacerLocation>& locations,
                         int num_threads = 1,
                         int score_threshold = 12);

// Extract spacers from genome reference using primary + alternative PAM in one pass.
bool spacerExtractionLoc(const GenomeReference& genome,
                         const std::string& pam,
                         const std::string& alt_pam,
                         int spacer_len,
                         bool upstream,
                         std::vector<SpacerLocation>& locations,
                         int num_threads = 1,
                         int score_threshold = 12);
