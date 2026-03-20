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
