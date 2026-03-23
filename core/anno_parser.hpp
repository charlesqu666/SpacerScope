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

#include <string>
#include <vector>
#include <cstdint>

#include "types.hpp"

struct AnnotationPlanSegment {
    uint32_t start_1based = 0;
    uint32_t end_1based = 0;
    uint32_t phase = 0;
};

struct AnnotationTranscriptPlan {
    std::string id;
    char strand = '+';
    std::vector<AnnotationPlanSegment> exons;
    std::vector<AnnotationPlanSegment> cds;
};

struct AnnotationQueryPlan {
    std::string gene_id;
    std::string chrom;
    uint32_t gene_start_1based = 0;
    uint32_t gene_end_1based = 0;
    char strand = '+';
    std::vector<AnnotationTranscriptPlan> transcripts;
};

bool build_annotation_query_plan(const std::string& gene_id,
                                 const std::string& annotation_file,
                                 AnnotationQueryPlan& out_plan,
                                 std::string* error_message = nullptr);

bool materialize_annotation_sequences(const AnnotationQueryPlan& plan,
                                      std::string_view chrom_sequence,
                                      const std::string& seq_type,
                                      SequenceList& out_sequences,
                                      std::string* error_message = nullptr);

bool extract_annotation_sequences(const std::string& gene_id,
                                  const std::string& annotation_file,
                                  const GenomeReference& genome,
                                  const std::string& seq_type,
                                  SequenceList& out_sequences,
                                  std::string* error_message = nullptr);
