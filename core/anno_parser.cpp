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
#include "anno_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Segment {
    std::string chrom;
    uint32_t start_1based = 0;
    uint32_t end_1based = 0;
    char strand = '+';
    uint32_t phase = 0;
};

struct TranscriptRecord {
    std::vector<Segment> exons;
    std::vector<Segment> cds;
    char strand = '+';
};

struct GeneRecord {
    std::string chrom;
    uint32_t start_1based = 0;
    uint32_t end_1based = 0;
    char strand = '+';
};

using TranscriptMap = std::unordered_map<std::string, TranscriptRecord>;

std::string to_lower_ascii(std::string value);
std::vector<std::string> split_tab_line(const std::string& line);
std::string get_attribute_value(std::string_view attributes, std::string_view key);
std::string normalize_gene_from_parent(const std::string& parent);

bool parse_annotation_records(const std::string& annotation_file,
                              std::unordered_map<std::string, TranscriptMap>& data,
                              std::unordered_map<std::string, GeneRecord>& gene_records,
                              std::string* error_message) {
    std::ifstream in(annotation_file);
    if (!in.is_open()) {
        if (error_message) *error_message = "Could not open annotation file: " + annotation_file;
        return false;
    }

    data.clear();
    gene_records.clear();
    data.reserve(1024);
    gene_records.reserve(1024);

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> cols = split_tab_line(line);
        if (cols.size() < 9) continue;

        const std::string feature = to_lower_ascii(cols[2]);
        const std::string chrom = cols[0];
        const uint32_t start = static_cast<uint32_t>(std::stoul(cols[3]));
        const uint32_t end = static_cast<uint32_t>(std::stoul(cols[4]));
        const char strand = cols[6].empty() ? '+' : cols[6][0];
        const uint32_t phase = (cols[7].empty() || cols[7] == ".") ? 0u : static_cast<uint32_t>(std::stoul(cols[7]));
        const std::string_view attrs(cols[8]);

        if (feature == "gene") {
            std::string gid = get_attribute_value(attrs, "ID");
            if (gid.empty()) gid = get_attribute_value(attrs, "gene_id");
            if (gid.empty()) gid = get_attribute_value(attrs, "gene_name");
            if (gid.empty()) gid = get_attribute_value(attrs, "Name");
            if (gid.empty()) continue;
            gene_records[gid] = GeneRecord{chrom, start, end, strand};
            data.try_emplace(gid);
            continue;
        }

        if (feature == "mrna" || feature == "transcript") {
            std::string tid = get_attribute_value(attrs, "ID");
            if (tid.empty()) tid = get_attribute_value(attrs, "transcript_id");
            if (tid.empty()) continue;

            std::string gid = get_attribute_value(attrs, "Parent");
            if (gid.empty()) gid = get_attribute_value(attrs, "gene_id");
            if (gid.empty()) gid = normalize_gene_from_parent(tid);
            gid = normalize_gene_from_parent(gid);

            data[gid].try_emplace(tid);
            data[gid][tid].strand = strand;
            continue;
        }

        if (feature != "exon" && feature != "cds") continue;

        std::string parent = get_attribute_value(attrs, "Parent");
        if (parent.empty()) parent = get_attribute_value(attrs, "transcript_id");
        if (parent.empty()) continue;

        std::string gid = get_attribute_value(attrs, "gene_id");
        if (gid.empty()) gid = normalize_gene_from_parent(parent);

        TranscriptRecord& record = data[gid][parent];
        if (record.strand == '+') record.strand = strand;

        Segment segment{chrom, start, end, strand, phase};
        if (feature == "exon") {
            record.exons.push_back(segment);
        } else {
            record.cds.push_back(segment);
        }
    }

    return true;
}

std::string_view safe_slice_view(std::string_view chrom_sequence,
                                 uint32_t start_1based,
                                 uint32_t end_1based) {
    if (start_1based == 0 || end_1based < start_1based) return "";
    const size_t start = static_cast<size_t>(start_1based - 1);
    const size_t len = static_cast<size_t>(end_1based - start_1based + 1);
    if (start + len > chrom_sequence.size()) return "";
    return chrom_sequence.substr(start, len);
}

void append_segment_sequence_view(std::string& out,
                                  std::string_view chrom_sequence,
                                  const AnnotationPlanSegment& segment,
                                  uint32_t trim_from_left,
                                  bool& ok) {
    std::string_view seq = safe_slice_view(chrom_sequence, segment.start_1based, segment.end_1based);
    if (seq.empty()) {
        ok = false;
        return;
    }
    if (trim_from_left >= seq.size()) return;
    out.append(seq.substr(trim_from_left));
}

std::string to_lower_ascii(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

char complement_base(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

std::string reverse_complement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc.push_back(complement_base(*it));
    }
    return rc;
}

std::string trim_copy(std::string_view text) {
    size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) ++start;
    size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
    return std::string(text.substr(start, end - start));
}

std::vector<std::string> split_tab_line(const std::string& line) {
    std::vector<std::string> cols;
    cols.reserve(9);
    size_t start = 0;
    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == '\t') {
            cols.emplace_back(line.substr(start, i - start));
            start = i + 1;
        }
    }
    return cols;
}

std::string get_attribute_value(std::string_view attributes, std::string_view key) {
    size_t pos = 0;
    while (pos < attributes.size()) {
        while (pos < attributes.size() &&
               (attributes[pos] == ';' || std::isspace(static_cast<unsigned char>(attributes[pos])))) {
            ++pos;
        }
        if (pos >= attributes.size()) break;

        size_t key_start = pos;
        while (pos < attributes.size() &&
               attributes[pos] != '=' &&
               !std::isspace(static_cast<unsigned char>(attributes[pos]))) {
            ++pos;
        }
        const std::string_view current_key = attributes.substr(key_start, pos - key_start);

        while (pos < attributes.size() && std::isspace(static_cast<unsigned char>(attributes[pos]))) ++pos;

        bool has_equals = false;
        if (pos < attributes.size() && attributes[pos] == '=') {
            has_equals = true;
            ++pos;
        }

        while (pos < attributes.size() && std::isspace(static_cast<unsigned char>(attributes[pos]))) ++pos;

        std::string value;
        if (pos < attributes.size() && attributes[pos] == '"') {
            ++pos;
            size_t end = attributes.find('"', pos);
            if (end == std::string_view::npos) end = attributes.size();
            value = std::string(attributes.substr(pos, end - pos));
            pos = (end < attributes.size()) ? end + 1 : end;
        } else {
            size_t value_start = pos;
            while (pos < attributes.size() &&
                   attributes[pos] != ';' &&
                   !(std::isspace(static_cast<unsigned char>(attributes[pos])) && !has_equals)) {
                ++pos;
            }
            value = trim_copy(attributes.substr(value_start, pos - value_start));
        }

        if (current_key == key) return value;

        while (pos < attributes.size() && attributes[pos] != ';') ++pos;
        if (pos < attributes.size() && attributes[pos] == ';') ++pos;
    }
    return "";
}

std::string normalize_gene_from_parent(const std::string& parent) {
    const std::string suffix = "-mRNA";
    const size_t pos = parent.find(suffix);
    if (pos != std::string::npos) {
        return parent.substr(0, pos);
    }
    return parent;
}

std::string_view safe_get_sequence(const GenomeReference& genome,
                                   const std::string& chrom,
                                   uint32_t start_1based,
                                   uint32_t end_1based) {
    const size_t chrom_idx = genome.get_chrom_idx(chrom);
    if (chrom_idx == static_cast<size_t>(-1)) return "";
    if (start_1based == 0 || end_1based < start_1based) return "";
    return genome.get_sequence(chrom_idx,
                               static_cast<size_t>(start_1based - 1),
                               static_cast<size_t>(end_1based - start_1based + 1));
}

void append_segment_sequence(std::string& out,
                             const GenomeReference& genome,
                             const Segment& segment,
                             uint32_t trim_from_left,
                             bool& ok) {
    std::string_view seq = safe_get_sequence(genome, segment.chrom, segment.start_1based, segment.end_1based);
    if (seq.empty()) {
        ok = false;
        return;
    }
    if (trim_from_left >= seq.size()) return;
    out.append(seq.substr(trim_from_left));
}

}  // namespace

bool build_annotation_query_plan(const std::string& gene_id,
                                 const std::string& annotation_file,
                                 AnnotationQueryPlan& out_plan,
                                 std::string* error_message) {
    out_plan = AnnotationQueryPlan{};

    std::unordered_map<std::string, TranscriptMap> data;
    std::unordered_map<std::string, GeneRecord> gene_records;
    if (!parse_annotation_records(annotation_file, data, gene_records, error_message)) {
        return false;
    }

    auto gene_it = data.find(gene_id);
    if (gene_it == data.end()) {
        if (error_message) *error_message = "gene_id '" + gene_id + "' not found in annotation";
        return false;
    }

    auto coord_it = gene_records.find(gene_id);
    if (coord_it == gene_records.end()) {
        if (error_message) *error_message = "gene feature for gene_id '" + gene_id + "' not found";
        return false;
    }

    out_plan.gene_id = gene_id;
    out_plan.chrom = coord_it->second.chrom;
    out_plan.gene_start_1based = coord_it->second.start_1based;
    out_plan.gene_end_1based = coord_it->second.end_1based;
    out_plan.strand = coord_it->second.strand;

    std::vector<std::pair<std::string, TranscriptRecord*>> ordered_transcripts;
    ordered_transcripts.reserve(gene_it->second.size());
    for (auto& entry : gene_it->second) {
        ordered_transcripts.push_back({entry.first, &entry.second});
    }
    std::sort(ordered_transcripts.begin(), ordered_transcripts.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    out_plan.transcripts.reserve(ordered_transcripts.size());
    for (const auto& entry : ordered_transcripts) {
        const std::string& tx_id = entry.first;
        const TranscriptRecord& record = *entry.second;
        AnnotationTranscriptPlan tx_plan;
        tx_plan.id = tx_id;
        tx_plan.strand = record.strand;

        tx_plan.exons.reserve(record.exons.size());
        for (const Segment& exon : record.exons) {
            if (exon.chrom != out_plan.chrom) {
                if (error_message) *error_message = "Transcript exon chromosome mismatch for '" + tx_id + "'";
                return false;
            }
            tx_plan.exons.push_back({exon.start_1based, exon.end_1based, 0});
        }

        tx_plan.cds.reserve(record.cds.size());
        for (const Segment& cds : record.cds) {
            if (cds.chrom != out_plan.chrom) {
                if (error_message) *error_message = "Transcript CDS chromosome mismatch for '" + tx_id + "'";
                return false;
            }
            tx_plan.cds.push_back({cds.start_1based, cds.end_1based, cds.phase});
        }

        out_plan.transcripts.push_back(std::move(tx_plan));
    }

    return true;
}

bool materialize_annotation_sequences(const AnnotationQueryPlan& plan,
                                      std::string_view chrom_sequence,
                                      const std::string& seq_type,
                                      SequenceList& out_sequences,
                                      std::string* error_message) {
    out_sequences.headers.clear();
    out_sequences.sequences.clear();

    const std::string seq_type_lower = to_lower_ascii(seq_type);
    if (seq_type_lower != "gene" && seq_type_lower != "mrna" && seq_type_lower != "cds") {
        if (error_message) {
            *error_message = "Unsupported --seq-type: " + seq_type + " (expected gene/mrna/cds)";
        }
        return false;
    }

    if (seq_type_lower == "gene") {
        std::string_view seq = safe_slice_view(chrom_sequence, plan.gene_start_1based, plan.gene_end_1based);
        if (seq.empty()) {
            if (error_message) *error_message = "Could not extract gene sequence for '" + plan.gene_id + "'";
            return false;
        }
        std::string final_seq = (plan.strand == '-') ? reverse_complement(seq) : std::string(seq);
        out_sequences.push_back(plan.gene_id, final_seq);
        return true;
    }

    for (const auto& tx_plan : plan.transcripts) {
        std::string concat;
        bool ok = true;

        if (seq_type_lower == "mrna") {
            if (tx_plan.exons.empty()) continue;
            std::vector<AnnotationPlanSegment> exons = tx_plan.exons;
            std::sort(exons.begin(), exons.end(),
                      [](const AnnotationPlanSegment& lhs, const AnnotationPlanSegment& rhs) {
                          return lhs.start_1based < rhs.start_1based;
                      });
            for (const auto& exon : exons) {
                append_segment_sequence_view(concat, chrom_sequence, exon, 0, ok);
                if (!ok) break;
            }
        } else {
            if (tx_plan.cds.empty()) continue;
            std::vector<AnnotationPlanSegment> cds = tx_plan.cds;
            std::sort(cds.begin(), cds.end(),
                      [](const AnnotationPlanSegment& lhs, const AnnotationPlanSegment& rhs) {
                          return lhs.start_1based < rhs.start_1based;
                      });
            for (const auto& cds_part : cds) {
                append_segment_sequence_view(concat, chrom_sequence, cds_part, cds_part.phase, ok);
                if (!ok) break;
            }
        }

        if (!ok) {
            if (error_message) *error_message = "Failed to extract transcript sequence for '" + tx_plan.id + "'";
            return false;
        }
        if (concat.empty()) continue;
        if (tx_plan.strand == '-') concat = reverse_complement(concat);
        out_sequences.push_back(tx_plan.id, concat);
    }

    if (out_sequences.size() == 0) {
        if (error_message) {
            *error_message = "No sequences extracted for gene_id '" + plan.gene_id + "' and seq-type '" + seq_type + "'";
        }
        return false;
    }
    return true;
}

bool extract_annotation_sequences(const std::string& gene_id,
                                  const std::string& annotation_file,
                                  const GenomeReference& genome,
                                  const std::string& seq_type,
                                  SequenceList& out_sequences,
                                  std::string* error_message) {
    out_sequences.headers.clear();
    out_sequences.sequences.clear();
    AnnotationQueryPlan plan;
    if (!build_annotation_query_plan(gene_id, annotation_file, plan, error_message)) {
        return false;
    }
    const size_t chrom_idx = genome.get_chrom_idx(plan.chrom);
    if (chrom_idx == static_cast<size_t>(-1) || chrom_idx >= genome.chrom_sequences.size()) {
        if (error_message) *error_message = "Reference chromosome '" + plan.chrom + "' not found";
        return false;
    }
    return materialize_annotation_sequences(plan, genome.chrom_sequences[chrom_idx], seq_type, out_sequences, error_message);
}
