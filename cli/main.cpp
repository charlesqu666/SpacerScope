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
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <thread>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <map>
#include <iomanip>
#include <cstdint>
#include "types.hpp"
#include "fasta_parser.hpp"
#include "spacer_extraction.hpp"
#include "dust_filter.hpp"
#include "count_filter.hpp"
#include "anno_parser.hpp"
#include "mismatch_search.hpp"
#include "indel_search.hpp"
#include "on_target_filter.hpp"
#include "postprocess.hpp"
#include "progress.hpp"
#include "report_html.hpp"
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#ifndef SPACERSCOPE_VERSION
#define SPACERSCOPE_VERSION "dev"
#endif

// --- Helper Functions ---

// OPTIMIZATION: Parse FASTA directly into GenomeReference (single load, zero-copy)
GenomeReference parse_fasta_to_genome(const std::string& filename) {
    GenomeReference genome;
    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string sequence_str(seq->seq.s);
        std::transform(sequence_str.begin(), sequence_str.end(), sequence_str.begin(), ::toupper);
        genome.add_chromosome(seq->name.s, std::move(sequence_str));
    }

    kseq_destroy(seq);
    gzclose(fp);
    return genome;
}

template <typename Fn>
void for_each_fasta_record(const std::string& filename, Fn&& fn) {
    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Could not open file " + filename);
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string sequence_str(seq->seq.s);
        std::transform(sequence_str.begin(), sequence_str.end(), sequence_str.begin(), ::toupper);
        if (!fn(std::string(seq->name.s), std::move(sequence_str))) {
            break;
        }
    }

    kseq_destroy(seq);
    gzclose(fp);
}

struct FastaSizeStats {
    std::uint64_t total_bases = 0;
    std::uint64_t max_chrom_bases = 0;
    std::uint64_t chrom_count = 0;
};

FastaSizeStats scan_fasta_size_stats(const std::string& filename) {
    FastaSizeStats stats;
    for_each_fasta_record(filename, [&](const std::string&, std::string&& chrom_seq) {
        const std::uint64_t len = static_cast<std::uint64_t>(chrom_seq.size());
        stats.total_bases += len;
        stats.max_chrom_bases = std::max(stats.max_chrom_bases, len);
        ++stats.chrom_count;
        return true;
    });
    return stats;
}

void load_annotation_queries(const std::string& gene_id,
                             const std::string& annotation_file,
                             const std::string& seq_type,
                             const GenomeReference& ref_genome,
                             GenomeReference& query_genome) {
    SequenceList extracted_queries;
    std::string error_message;
    if (!extract_annotation_sequences(gene_id, annotation_file, ref_genome, seq_type, extracted_queries, &error_message)) {
        throw std::runtime_error(error_message);
    }

    for (size_t i = 0; i < extracted_queries.size(); ++i) {
        query_genome.add_chromosome(extracted_queries.headers[i], std::move(extracted_queries.sequences[i]));
    }
}

void build_query_spacers(const GenomeReference& query_genome,
                         const std::string& pam,
                         int spacer_len,
                         bool upstream,
                         int num_threads,
                         int dust_threshold,
                         std::vector<SpacerLocation>& query_locs) {
    std::vector<SpacerLocation> unfiltered_locs;
    spacerExtractionLoc(query_genome, pam, spacer_len, upstream, unfiltered_locs, num_threads, dust_threshold);
    query_locs = std::move(unfiltered_locs);
}

void load_annotation_queries_cut(const AnnotationQueryPlan& plan,
                                 const std::string& ref_file,
                                 const std::string& seq_type,
                                 GenomeReference& query_genome) {
    bool found_chrom = false;
    bool extracted = false;
    std::string error_message;

    for_each_fasta_record(ref_file, [&](const std::string& chrom_name, std::string&& chrom_seq) {
        if (chrom_name != plan.chrom) {
            return true;
        }
        found_chrom = true;
        SequenceList extracted_queries;
        if (!materialize_annotation_sequences(plan, chrom_seq, seq_type, extracted_queries, &error_message)) {
            throw std::runtime_error(error_message);
        }
        for (size_t i = 0; i < extracted_queries.size(); ++i) {
            query_genome.add_chromosome(extracted_queries.headers[i], std::move(extracted_queries.sequences[i]));
        }
        extracted = true;
        return false;
    });

    if (!found_chrom) {
        throw std::runtime_error("Reference chromosome '" + plan.chrom + "' not found for anno-cut mode");
    }
    if (!extracted || query_genome.size() == 0) {
        throw std::runtime_error("No annotation-derived query sequence was extracted in anno-cut mode");
    }
}

void append_hits(std::vector<SearchHit>& dst, std::vector<SearchHit>& src) {
    if (src.empty()) return;
    dst.insert(dst.end(),
               std::make_move_iterator(src.begin()),
               std::make_move_iterator(src.end()));
}

size_t choose_query_processing_batch_size(size_t total_queries, bool cut_mode) {
    if (total_queries == 0) return 0;
    if (total_queries <= 128) return total_queries;
    if (total_queries <= 2048) return 256;
    if (total_queries <= 8192) return cut_mode ? 256 : 128;
    if (total_queries <= 32768) return 128;
    return 64;
}

void write_raw_results_header(std::ostream& out_stream) {
    out_stream << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\tSearchType\tGenomicFrequency\n";
}

void append_raw_results_tsv(std::ostream& out_stream,
                            const std::vector<CombinedResult>& combined_results) {
    for (const auto& data : combined_results) {
        out_stream << data.query_name << "\t" << data.query_seq << "\t"
                   << data.target_name << "\t" << data.target_seq << "\t"
                   << data.distance << "\t" << data.search_type << "\t"
                   << data.genomic_frequency << "\n";
    }
}

void write_scored_results_header(std::ostream& out_stream) {
    out_stream << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\tSearchType\tGenomicFrequency\tScore\tDetails\n";
}

void append_scored_results_tsv(std::ostream& out_stream,
                               const std::vector<ScoredResult>& results) {
    for (const auto& result : results) {
        out_stream << result.base.query_name << '\t'
                   << result.base.query_seq << '\t'
                   << result.base.target_name << '\t'
                   << result.base.target_seq << '\t'
                   << result.base.distance << '\t'
                   << result.base.search_type << '\t'
                   << result.base.genomic_frequency << '\t'
                   << format_score_cpp(result.score) << '\t'
                   << alignment_events_to_json(result.details) << '\n';
    }
}

void process_reference_chunk_collect(const GenomeReference& chunk_genome,
                                     const std::vector<SpacerLocation>& query_locs,
                                     const GenomeReference& query_genome,
                                     const std::string& pam,
                                     const std::string& alt_pam,
                                     int spacer_len,
                                     bool upstream,
                                     int max_mismatch,
                                     int max_indel,
                                     int num_threads,
                                     int dust_threshold,
                                     size_t& total_target_spacers,
                                     double& step2_seconds,
                                     double& step3_seconds,
                                     double& step4_seconds,
                                     bool count_target_spacers,
                                     std::vector<SearchHit>& mismatch_hits,
                                     std::vector<SearchHit>& indel_hits) {
    auto step_start_time = std::chrono::high_resolution_clock::now();
    std::vector<SpacerLocation> target_locs;
    spacerExtractionLoc(chunk_genome, pam, alt_pam, spacer_len, upstream, target_locs, num_threads, dust_threshold);
    step2_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
    if (count_target_spacers) total_target_spacers += target_locs.size();

    if (max_mismatch >= 0) {
        step_start_time = std::chrono::high_resolution_clock::now();
        std::vector<SearchHit> chunk_mismatch = process_mismatch_search_loc_collect(query_locs,
                                                                                     target_locs,
                                                                                     query_genome,
                                                                                     chunk_genome,
                                                                                     max_mismatch,
                                                                                     num_threads);
        step3_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
        append_hits(mismatch_hits, chunk_mismatch);
    }

    if (max_indel > 0) {
        step_start_time = std::chrono::high_resolution_clock::now();
        std::vector<SearchHit> chunk_indel = process_indel_search_loc_collect(query_locs,
                                                                               target_locs,
                                                                               query_genome,
                                                                               chunk_genome,
                                                                               max_indel,
                                                                               num_threads,
                                                                               spacer_len);
        step4_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
        append_hits(indel_hits, chunk_indel);
    }
}

// Combines mismatch/indel temporary hits into a deduplicated in-memory result set.
std::vector<CombinedResult> combine_results_direct(
        const std::vector<SearchHit>& mismatch_hits,
        const std::vector<SearchHit>& indel_hits,
        int max_mismatch,
        int max_indel) {

    std::map<std::pair<std::string, std::string>, CombinedResult> unique_results;

    auto store_results = [&](const std::vector<SearchHit>& hits, const std::string& type_label) {
        for (const auto& raw_hit : hits) {
            CombinedResult data = raw_hit;
            std::pair<std::string, std::string> key = {data.query_name, data.target_name};
            if (unique_results.count(key)) continue;

            if (type_label == "Mismatch") {
                data.search_type = (data.distance == 0) ? "Exact" : "Mismatch";
            } else {
                data.search_type = "Indel";
            }
            unique_results[key] = data;
        }
    };

    if (max_mismatch >= 0) store_results(mismatch_hits, "Mismatch");
    if (max_indel > 0) store_results(indel_hits, "Indel");

    std::map<std::pair<std::string, std::string>, int> new_frequency_map;
    for (const auto& pair : unique_results) {
        const auto& data = pair.second;
        new_frequency_map[{data.query_name, data.target_seq}]++;
    }

    std::vector<CombinedResult> combined_results;
    combined_results.reserve(unique_results.size());
    for (auto& pair : unique_results) {
        auto& data = pair.second;
        data.genomic_frequency = new_frequency_map.at({data.query_name, data.target_seq});
        combined_results.push_back(std::move(data));
    }
    return combined_results;
}

void write_raw_results_tsv(const std::string& output_path,
                           const std::vector<CombinedResult>& combined_results) {
    std::ofstream final_out_stream(output_path, std::ios::binary | std::ios::trunc);
    final_out_stream << "QueryName\tQuerySeq\tTargetName\tTargetSeq\tDistance\tSearchType\tGenomicFrequency\n";
    for (const auto& data : combined_results) {
        final_out_stream << data.query_name << "\t" << data.query_seq << "\t"
                         << data.target_name << "\t" << data.target_seq << "\t"
                         << data.distance << "\t" << data.search_type << "\t"
                         << data.genomic_frequency << "\n";
    }
    final_out_stream.flush();
}

enum class CliMode {
    Legacy,
    Fasta,
    FastaCut,
    Anno,
    AnnoCut,
};

bool is_mode_token(const std::string& arg) {
    return arg == "fasta" || arg == "fasta-cut" || arg == "anno" || arg == "anno-cut";
}

CliMode parse_mode_token(const std::string& arg) {
    if (arg == "fasta") return CliMode::Fasta;
    if (arg == "fasta-cut") return CliMode::FastaCut;
    if (arg == "anno") return CliMode::Anno;
    if (arg == "anno-cut") return CliMode::AnnoCut;
    return CliMode::Legacy;
}

bool mode_uses_query(CliMode mode) {
    return mode == CliMode::Fasta || mode == CliMode::FastaCut;
}

bool mode_uses_annotation(CliMode mode) {
    return mode == CliMode::Anno || mode == CliMode::AnnoCut;
}

bool is_cut_mode(CliMode mode) {
    return mode == CliMode::FastaCut || mode == CliMode::AnnoCut;
}

const char* mode_name(CliMode mode) {
    switch (mode) {
        case CliMode::Fasta: return "fasta";
        case CliMode::FastaCut: return "fasta-cut";
        case CliMode::Anno: return "anno";
        case CliMode::AnnoCut: return "anno-cut";
        default: return "legacy";
    }
}

double estimate_peak_memory_mib(CliMode mode, const FastaSizeStats& stats) {
    const double total_mbp = static_cast<double>(stats.total_bases) / 1'000'000.0;
    const double max_mbp = static_cast<double>(stats.max_chrom_bases) / 1'000'000.0;
    if (is_cut_mode(mode)) {
        return 19.676 + 0.123 * total_mbp + 9.268 * max_mbp;
    }
    return 19.320 + 8.640 * total_mbp - 0.132 * max_mbp;
}

std::string memory_formula_summary(CliMode mode) {
    if (is_cut_mode(mode)) {
        return "peak_mib ~= 19.7 + 9.27 * longest_chromosome_mbp";
    }
    return "peak_mib ~= 19.2 + 8.64 * total_genome_mbp";
}

FastaSizeStats emit_estimated_memory_usage(CliMode mode, const std::string& ref_file) {
    FastaSizeStats stats = scan_fasta_size_stats(ref_file);
    const double total_mbp = static_cast<double>(stats.total_bases) / 1'000'000.0;
    const double max_mbp = static_cast<double>(stats.max_chrom_bases) / 1'000'000.0;
    const double peak_mib = estimate_peak_memory_mib(mode, stats);
    const double peak_gib = peak_mib / 1024.0;
    progress::memory_estimate(stats.chrom_count,
                              total_mbp,
                              max_mbp,
                              peak_mib,
                              peak_gib,
                              memory_formula_summary(mode));
    return stats;
}

bool confirm_estimated_memory_usage(CliMode mode, const std::string& ref_file) {
    (void)emit_estimated_memory_usage(mode, ref_file);
    while (true) {
        std::ostream& prompt_stream = progress::is_jsonl() ? std::cerr : std::cout;
        prompt_stream << "Continue? Enter Y to proceed or N to abort: " << std::flush;
        std::string answer;
        if (!std::getline(std::cin, answer)) {
            std::cerr << "\nError: unable to read confirmation from stdin. "
                      << "Use --skip-memory-check to bypass this prompt." << std::endl;
            return false;
        }
        std::transform(answer.begin(), answer.end(), answer.begin(), ::toupper);
        if (answer == "Y" || answer == "YES") {
            return true;
        }
        if (answer == "N" || answer == "NO") {
            if (progress::is_jsonl()) {
                progress::log("Run canceled by user.");
            } else {
                std::cout << "Run canceled by user." << std::endl;
            }
            return false;
        }
        if (progress::is_jsonl()) {
            std::cerr << "Please enter Y or N." << std::endl;
        } else {
            std::cout << "Please enter Y or N." << std::endl;
        }
    }
}

void print_version(std::ostream& os = std::cout) {
    os << "spacerscope " << SPACERSCOPE_VERSION << "\n";
}

void print_general_usage(std::ostream& os = std::cout) {
    os << "Usage:\n"
       << "  spacerscope fasta <options>\n"
       << "  spacerscope fasta-cut <options>\n"
       << "  spacerscope anno <options>\n"
       << "  spacerscope anno-cut <options>\n"
       << "  spacerscope --report-tsv <file> --html-output <file>\n"
       << "  spacerscope <legacy options>\n\n"
       << "Modes:\n"
       << "  fasta       Query FASTA mode\n"
       << "  fasta-cut   Low-memory FASTA cut mode\n"
       << "  anno        Annotation-derived query mode\n"
       << "  anno-cut    Low-memory annotation cut mode\n\n"
       << "Common options:\n"
       << "  -R, --ref <file>         Reference FASTA file\n"
       << "  -o, --output <file>      Output TSV file\n"
       << "  --pam <string>           PAM sequence (default: NGG)\n"
       << "  --alt-pam <string>       Alternate PAM sequence (default: NAG)\n"
       << "  --spacer-len <int>       Spacer length (default: 20)\n"
       << "  --direction <up|down>    Search direction (default: up)\n"
       << "  -m, --mismatch <int>     Max mismatches (default: 4)\n"
       << "  -i, --indel <int>        Max indels (default: 2)\n"
       << "  --dust-threshold <int>   DUST complexity threshold (default: 12)\n"
       << "  -t, --threads <int>      Number of threads\n"
       << "  -y, --yes                Skip the startup memory confirmation prompt\n"
       << "  --skip-memory-check      Skip the startup memory confirmation prompt\n"
       << "  --raw-output             Write 7-column unscored merged TSV output\n"
       << "  --progress-format <fmt>  Progress output format: text | jsonl (default: text)\n"
       << "  --html <0|1>             Generate <output>.html when set to 1\n"
       << "  -ht <0|1>                Alias of --html\n"
       << "  --html-output <file>     Write HTML report to the given path\n"
       << "  --report-tsv <file>      Generate HTML from an existing scored TSV\n"
       << "  -h, --help               Show help\n"
       << "  -V, --version            Show version\n\n"
       << "Mode-specific options:\n"
       << "  fasta / fasta-cut:\n"
       << "    -I, --query <file>       Query FASTA file\n"
       << "  anno / anno-cut:\n"
       << "    -A, --annotation <file>  Annotation GFF/GFF3 file\n"
       << "    --geneID <id>            Exact gene ID from the annotation file\n"
       << "    --seq-type <type>        gene | mrna | cds (default: gene)\n\n"
       << "Notes:\n"
       << "  - anno/anno-cut are released with legacy annotation compatibility.\n"
       << "  - --geneID must match the annotation gene ID exactly; it is not a free gene symbol.\n"
       << "  - mrna/cds extraction follows the legacy GeneX-mRNA transcript naming convention.\n"
       << "  - Generic Parent=transcript:... GFF3 inputs are not guaranteed for mrna/cds mode.\n\n"
       << "Examples:\n"
       << "  spacerscope fasta --query q.fa --ref ref.fa -o out.tsv\n"
       << "  spacerscope fasta-cut --query q.fa --ref ref.fa -o out.tsv\n"
       << "  spacerscope anno --annotation genes.gff3 --geneID gene:Gene1 --seq-type cds --ref ref.fa -o out.tsv\n"
       << "  spacerscope --report-tsv out.tsv --html-output out.html\n";
}

void print_mode_usage(CliMode mode, std::ostream& os = std::cout) {
    os << "Mode: " << mode_name(mode) << "\n\n";
    switch (mode) {
        case CliMode::Fasta:
            os << "Usage:\n"
               << "  spacerscope fasta --query <file> --ref <file> -o <file> [common options]\n\n"
               << "Description:\n"
               << "  Search query FASTA sequences directly against the reference genome.\n";
            break;
        case CliMode::FastaCut:
            os << "Usage:\n"
               << "  spacerscope fasta-cut --query <file> --ref <file> -o <file> [common options]\n\n"
               << "Description:\n"
               << "  Low-memory mode. Streams the reference genome chromosome by chromosome.\n";
            break;
        case CliMode::Anno:
            os << "Usage:\n"
               << "  spacerscope anno --annotation <file> --geneID <id> --ref <file> -o <file> [common options]\n\n"
               << "Description:\n"
               << "  Extract query sequence(s) from annotation, then run the standard pipeline.\n"
               << "  --geneID must match the annotation gene ID exactly.\n"
               << "  mrna/cds extraction follows the legacy GeneX-mRNA transcript naming convention.\n"
               << "  Generic Parent=transcript:... GFF3 inputs are not guaranteed for mrna/cds mode.\n";
            break;
        case CliMode::AnnoCut:
            os << "Usage:\n"
               << "  spacerscope anno-cut --annotation <file> --geneID <id> --ref <file> -o <file> [common options]\n\n"
               << "Description:\n"
               << "  Extract annotation-derived query sequence(s), then run low-memory cut mode.\n"
               << "  --geneID must match the annotation gene ID exactly.\n"
               << "  mrna/cds extraction follows the legacy GeneX-mRNA transcript naming convention.\n"
               << "  Generic Parent=transcript:... GFF3 inputs are not guaranteed for mrna/cds mode.\n";
            break;
        default:
            print_general_usage(os);
            return;
    }
    os << "\n";
}

bool is_help_token(const std::string& arg) {
    return arg == "-h" || arg == "--help" || arg == "-help";
}

bool is_version_token(const std::string& arg) {
    return arg == "-V" || arg == "--version" || arg == "-version";
}

// --- Main Pipeline ---

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_general_usage(std::cerr);
        return 1;
    }

    {
        const std::string first = argv[1];
        if (is_help_token(first)) {
            print_general_usage();
            return 0;
        }
        if (is_version_token(first)) {
            print_version();
            return 0;
        }
        if (first == "help") {
            if (argc >= 3 && is_mode_token(argv[2])) {
                print_mode_usage(parse_mode_token(argv[2]));
            } else {
                print_general_usage();
            }
            return 0;
        }
    }

    CliMode cli_mode = CliMode::Legacy;
    int arg_start = 1;
    if (argc >= 2 && is_mode_token(argv[1])) {
        cli_mode = parse_mode_token(argv[1]);
        arg_start = 2;
    }

    std::string query_file, ref_file, annotation_file, gene_id, output_file;
    std::string html_output_file, report_tsv_file;
    std::string pam = "NGG", alt_pam = "NAG", direction_str = "up";
    std::string seq_type = "gene";
    int spacer_len = 20, max_mismatch = 4, max_indel = 2;
    int num_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (num_threads <= 0) num_threads = 1;
    int dust_threshold = 12;
    bool raw_output = false;
    bool skip_memory_check = false;
    int html_flag = -1;
    std::string progress_format_str = "text";

    // Command line parsing
    for (int i = arg_start; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const std::string& opt) -> bool {
            if (i + 1 >= argc) {
                std::cerr << "Error: option '" << opt << "' requires a value." << std::endl;
                return false;
            }
            return true;
        };
        if (is_help_token(arg)) {
            print_mode_usage(cli_mode);
            return 0;
        }
        if (is_version_token(arg)) {
            print_version();
            return 0;
        }
        if (arg == "-I" || arg == "--query") {
            if (!require_value(arg)) return 1;
            query_file = argv[++i];
        }
        else if (arg == "-R" || arg == "--ref") {
            if (!require_value(arg)) return 1;
            ref_file = argv[++i];
        }
        else if (arg == "-A" || arg == "--annotation") {
            if (!require_value(arg)) return 1;
            annotation_file = argv[++i];
        }
        else if (arg == "--geneID") {
            if (!require_value(arg)) return 1;
            gene_id = argv[++i];
        }
        else if (arg == "--seq-type") {
            if (!require_value(arg)) return 1;
            seq_type = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (!require_value(arg)) return 1;
            output_file = argv[++i];
        }
        else if (arg == "--pam") {
            if (!require_value(arg)) return 1;
            pam = argv[++i];
        }
        else if (arg == "--alt-pam") {
            if (!require_value(arg)) return 1;
            alt_pam = argv[++i];
        }
        else if (arg == "--spacer-len") {
            if (!require_value(arg)) return 1;
            spacer_len = std::stoi(argv[++i]);
        }
        else if (arg == "--direction") {
            if (!require_value(arg)) return 1;
            direction_str = argv[++i];
        }
        else if (arg == "-m" || arg == "--mismatch") {
            if (!require_value(arg)) return 1;
            max_mismatch = std::stoi(argv[++i]);
        }
        else if (arg == "-i" || arg == "--indel") {
            if (!require_value(arg)) return 1;
            max_indel = std::stoi(argv[++i]);
        }
        else if (arg == "--dust-threshold") {
            if (!require_value(arg)) return 1;
            dust_threshold = std::stoi(argv[++i]);
        }
        else if (arg == "-t" || arg == "--threads") {
            if (!require_value(arg)) return 1;
            num_threads = std::stoi(argv[++i]);
        }
        else if (arg == "-y" || arg == "--yes" || arg == "--skip-memory-check") skip_memory_check = true;
        else if (arg == "--html" || arg == "-ht") {
            if (!require_value(arg)) return 1;
            html_flag = std::stoi(argv[++i]);
        }
        else if (arg == "--raw-output") raw_output = true;
        else if (arg == "--progress-format") {
            if (!require_value(arg)) return 1;
            progress_format_str = argv[++i];
        }
        else if (arg == "--html-output") {
            if (!require_value(arg)) return 1;
            html_output_file = argv[++i];
        }
        else if (arg == "--report-tsv") {
            if (!require_value(arg)) return 1;
            report_tsv_file = argv[++i];
        }
        else {
            std::cerr << "Error: unknown argument '" << arg << "'. Use --help for usage." << std::endl;
            return 1;
        }
    }

    std::transform(direction_str.begin(), direction_str.end(), direction_str.begin(), ::tolower);
    std::transform(seq_type.begin(), seq_type.end(), seq_type.begin(), ::tolower);
    std::transform(progress_format_str.begin(), progress_format_str.end(), progress_format_str.begin(), ::tolower);

    if (progress_format_str == "text") {
        progress::set_format(progress::Format::Text);
    } else if (progress_format_str == "jsonl") {
        progress::set_format(progress::Format::Jsonl);
    } else {
        std::cerr << "Error: --progress-format must be 'text' or 'jsonl'." << std::endl;
        return 1;
    }

    if (!report_tsv_file.empty()) {
        if (html_output_file.empty()) html_output_file = report_tsv_file + ".html";
        std::string report_error;
        progress::step_start("Report", "HTML report generation");
        if (!write_html_report_from_tsv(report_tsv_file, html_output_file, report_error)) {
            progress::error(report_error);
            return 1;
        }
        progress::log("HTML report generated: " + html_output_file);
        progress::step_done("Report", 0.0);
        progress::complete(0.0, report_tsv_file, html_output_file);
        return 0;
    }

    if (spacer_len <= 0) {
        std::cerr << "Error: --spacer-len must be > 0." << std::endl;
        return 1;
    }
    if (max_mismatch < 0) {
        std::cerr << "Error: --mismatch must be >= 0." << std::endl;
        return 1;
    }
    if (max_indel < 0) {
        std::cerr << "Error: --indel must be >= 0." << std::endl;
        return 1;
    }
    if (dust_threshold < 0) {
        std::cerr << "Error: --dust-threshold must be >= 0." << std::endl;
        return 1;
    }
    if (num_threads <= 0) {
        std::cerr << "Error: --threads must be > 0." << std::endl;
        return 1;
    }
    if (direction_str != "up" && direction_str != "down") {
        std::cerr << "Error: --direction must be 'up' or 'down'." << std::endl;
        return 1;
    }
    if (seq_type != "gene" && seq_type != "mrna" && seq_type != "cds") {
        std::cerr << "Error: --seq-type must be one of: gene, mrna, cds." << std::endl;
        return 1;
    }
    if (html_flag != -1 && html_flag != 0 && html_flag != 1) {
        std::cerr << "Error: --html / -ht must be 0 or 1." << std::endl;
        return 1;
    }

    if (html_flag == 1 && html_output_file.empty() && !output_file.empty()) {
        html_output_file = output_file + ".html";
    } else if (html_flag == 0) {
        html_output_file.clear();
    }

    if (ref_file.empty() || output_file.empty()) {
        std::cerr << "Error: Reference file and output file are required." << std::endl;
        return 1;
    }

    if (cli_mode == CliMode::Legacy) {
        if (query_file.empty() && (annotation_file.empty() || gene_id.empty())) {
            std::cerr << "Error: Provide either --query or both --annotation and --geneID." << std::endl;
            return 1;
        }
    } else if (mode_uses_query(cli_mode)) {
        if (!annotation_file.empty() || !gene_id.empty()) {
            std::cerr << "Error: fasta/fasta-cut mode does not accept --annotation or --geneID." << std::endl;
            return 1;
        }
        if (query_file.empty()) {
            std::cerr << "Error: mode requires --query." << std::endl;
            return 1;
        }
    } else if (mode_uses_annotation(cli_mode)) {
        if (!query_file.empty()) {
            std::cerr << "Error: anno/anno-cut mode does not accept --query." << std::endl;
            return 1;
        }
        if (annotation_file.empty() || gene_id.empty()) {
            std::cerr << "Error: anno/anno-cut mode requires both --annotation and --geneID." << std::endl;
            return 1;
        }
    }

    if (is_cut_mode(cli_mode)) {
        progress::log(std::string("[mode] ") +
                      (cli_mode == CliMode::FastaCut ? "fasta-cut" : "anno-cut") +
                      " will stream the reference genome chromosome by chromosome to reduce peak memory.");
    }

    if (!skip_memory_check) {
        CliMode estimate_mode = cli_mode;
        if (estimate_mode == CliMode::Legacy) {
            estimate_mode = query_file.empty() ? CliMode::Anno : CliMode::Fasta;
        }
        if (!confirm_estimated_memory_usage(estimate_mode, ref_file)) {
            return 0;
        }
    } else if (progress::is_jsonl()) {
        CliMode estimate_mode = cli_mode;
        if (estimate_mode == CliMode::Legacy) {
            estimate_mode = query_file.empty() ? CliMode::Anno : CliMode::Fasta;
        }
        (void)emit_estimated_memory_usage(estimate_mode, ref_file);
    }

    bool upstream = (direction_str == "up");
    auto total_start_time = std::chrono::high_resolution_clock::now();

    try {
        if (is_cut_mode(cli_mode)) {
            auto step_start_time = std::chrono::high_resolution_clock::now();
            progress::log("Cut mode: streaming reference genome chromosome by chromosome...");

            progress::step_start("Step 1", "Preparing initial on-target spacers");
            GenomeReference query_genome;
            std::vector<SpacerLocation> query_locs;

            if (mode_uses_query(cli_mode)) {
                query_genome = parse_fasta_to_genome(query_file);
            } else {
                AnnotationQueryPlan plan;
                std::string error_message;
                if (!build_annotation_query_plan(gene_id, annotation_file, plan, &error_message)) {
                    throw std::runtime_error(error_message);
                }
                load_annotation_queries_cut(plan, ref_file, seq_type, query_genome);
                progress::log("  Extracted " + std::to_string(query_genome.size()) +
                              " annotation-derived query sequence(s) for gene " + gene_id +
                              " (seq-type: " + seq_type + ")");
            }

            build_query_spacers(query_genome, pam, spacer_len, upstream, num_threads, dust_threshold, query_locs);
            progress::log("  Found " + std::to_string(query_locs.size()) + " query spacers");
            progress::step_done("Step 1",
                                std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count());

            const size_t query_batch_size = std::max<size_t>(1, choose_query_processing_batch_size(query_locs.size(), true));
            if (query_batch_size < query_locs.size()) {
                progress::log("  Low-memory query batching enabled (" +
                              std::to_string(query_batch_size) + " query spacers per batch)");
            }

            std::ofstream final_out_stream(output_file, std::ios::binary | std::ios::trunc);
            if (!final_out_stream.is_open()) {
                throw std::runtime_error("Could not open output file: " + output_file);
            }
            if (raw_output) {
                if (!html_output_file.empty()) {
                    std::cerr << "Error: --html-output requires scored output, not --raw-output." << std::endl;
                    return 1;
                }
                write_raw_results_header(final_out_stream);
            } else {
                write_scored_results_header(final_out_stream);
            }

            double step2_seconds = 0.0;
            double step3_seconds = 0.0;
            double step4_seconds = 0.0;
            double step5_seconds = 0.0;
            double step6_seconds = 0.0;
            size_t total_target_spacers = 0;
            size_t total_discarded = 0;

            for (size_t batch_start = 0; batch_start < query_locs.size(); batch_start += query_batch_size) {
                const size_t batch_count = std::min(query_batch_size, query_locs.size() - batch_start);
                std::vector<SpacerLocation> batch_query_locs(query_locs.begin() + static_cast<std::ptrdiff_t>(batch_start),
                                                             query_locs.begin() + static_cast<std::ptrdiff_t>(batch_start + batch_count));
                std::vector<SearchHit> mismatch_hits;
                std::vector<SearchHit> indel_hits;
                size_t chrom_count = 0;
                const bool count_target_spacers = (batch_start == 0);

                for_each_fasta_record(ref_file, [&](const std::string& chrom_name, std::string&& chrom_seq) {
                    GenomeReference chunk_genome;
                    chunk_genome.add_chromosome(chrom_name, std::move(chrom_seq));
                    ++chrom_count;
                    progress::log("  [cut] processing chromosome " + chrom_name + " (" + std::to_string(chrom_count) +
                                  "), query batch " + std::to_string(batch_start + batch_count) + "/" +
                                  std::to_string(query_locs.size()));
                    process_reference_chunk_collect(chunk_genome,
                                                   batch_query_locs,
                                                   query_genome,
                                                   pam,
                                                   alt_pam,
                                                   spacer_len,
                                                   upstream,
                                                   max_mismatch,
                                                   max_indel,
                                                   num_threads,
                                                   dust_threshold,
                                                   total_target_spacers,
                                                   step2_seconds,
                                                   step3_seconds,
                                                   step4_seconds,
                                                   count_target_spacers,
                                                   mismatch_hits,
                                                   indel_hits);
                    return true;
                });

                step_start_time = std::chrono::high_resolution_clock::now();
                std::vector<CombinedResult> combined_results =
                    combine_results_direct(mismatch_hits, indel_hits, max_mismatch, max_indel);
                step5_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();

                if (raw_output) {
                    append_raw_results_tsv(final_out_stream, combined_results);
                } else {
                    step_start_time = std::chrono::high_resolution_clock::now();
                    std::unordered_map<std::string, std::vector<CombinedResult>> results_by_chrom;
                    results_by_chrom.reserve(128);
                    for (const auto& result : combined_results) {
                        RegionHeader region;
                        if (!parse_region_header(result.target_name, region)) {
                            continue;
                        }
                        results_by_chrom[region.chrom].push_back(result);
                    }

                    size_t batch_discarded = 0;
                    std::vector<ScoredResult> scored_results;
                    for_each_fasta_record(ref_file, [&](const std::string& chrom_name, std::string&& chrom_seq) {
                        auto it = results_by_chrom.find(chrom_name);
                        if (it == results_by_chrom.end() || it->second.empty()) {
                            return true;
                        }
                        GenomeReference chunk_genome;
                        chunk_genome.add_chromosome(chrom_name, std::move(chrom_seq));
                        size_t chunk_discarded = 0;
                        std::vector<ScoredResult> chunk_scored = score_results_cpp(it->second,
                                                                                   chunk_genome,
                                                                                   max_indel,
                                                                                   max_mismatch,
                                                                                   upstream ? -1 : 1,
                                                                                   num_threads,
                                                                                   chunk_discarded);
                        batch_discarded += chunk_discarded;
                        scored_results.insert(scored_results.end(),
                                              std::make_move_iterator(chunk_scored.begin()),
                                              std::make_move_iterator(chunk_scored.end()));
                        return true;
                    });
                    total_discarded += batch_discarded;
                    append_scored_results_tsv(final_out_stream, scored_results);
                    step6_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
                }
            }

            final_out_stream.flush();

            progress::step_start("Step 2", "Extracting reference spacer locations");
            progress::log("  Found " + std::to_string(total_target_spacers) + " target spacers");
            progress::step_done("Step 2", step2_seconds);

            if (max_mismatch >= 0) {
                progress::step_start("Step 3", "Mismatch search (Location-based, Parallel Indexing)");
                progress::step_done("Step 3", step3_seconds);
            }

            if (max_indel > 0) {
                progress::step_start("Step 4", "Indel search (Location-based)");
                progress::step_done("Step 4", step4_seconds);
            }

            progress::step_start("Step 5", "Combining results");
            progress::step_done("Step 5", step5_seconds);

            if (!raw_output) {
                progress::step_start("Step 6", "Post-processing and scoring");
                progress::step_done("Step 6", step6_seconds);
                if (total_discarded > 0) {
                    progress::log("  Note: discarded " + std::to_string(total_discarded) +
                                  " records without a valid alignment path during post-processing.");
                }
                if (!html_output_file.empty()) {
                    step_start_time = std::chrono::high_resolution_clock::now();
                    progress::step_start("Step 7", "HTML report generation");
                    std::string report_error;
                    if (!write_html_report_from_tsv(output_file, html_output_file, report_error)) {
                        progress::error(report_error);
                        return 1;
                    }
                    progress::log("  HTML report written to " + html_output_file);
                    progress::step_done("Step 7",
                                        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count());
                }
            }

            progress::complete(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - total_start_time).count(),
                               output_file,
                               html_output_file);
            return 0;
        }

        // OPTIMIZATION: Load reference genome once into memory
        auto step_start_time = std::chrono::high_resolution_clock::now();
        progress::log("Loading reference genome...");
        GenomeReference genome = parse_fasta_to_genome(ref_file);
        progress::log("  Loaded " + std::to_string(genome.size()) + " sequences");

        // Step 1: On-target Preparation (extract query spacer locations)
        progress::step_start("Step 1", "Preparing initial on-target spacers");
        step_start_time = std::chrono::high_resolution_clock::now();
        std::vector<SpacerLocation> query_locs;
        GenomeReference query_genome;  // Keep query genome alive for search

        if (!query_file.empty()) {
            // Parse query file into genome structure
            query_genome = parse_fasta_to_genome(query_file);
            std::vector<SpacerLocation> unfiltered_locs;
            spacerExtractionLoc(query_genome, pam, spacer_len, upstream, unfiltered_locs, num_threads, dust_threshold);
            
            // Apply DUST filter by checking complexity on-demand
            for (const auto& loc : unfiltered_locs) {
                std::string_view seq = query_genome.get_sequence(loc);
                // Simple complexity check: just keep all for now (dust_filter would need adaptation)
                query_locs.push_back(loc);
            }
        } else {
            load_annotation_queries(gene_id, annotation_file, seq_type, genome, query_genome);
            progress::log("  Extracted " + std::to_string(query_genome.size()) +
                          " annotation-derived query sequence(s) for gene " + gene_id +
                          " (seq-type: " + seq_type + ")");
            std::vector<SpacerLocation> unfiltered_locs;
            spacerExtractionLoc(query_genome, pam, spacer_len, upstream, unfiltered_locs, num_threads, dust_threshold);
            
            for (const auto& loc : unfiltered_locs) {
                query_locs.push_back(loc);
            }
        }
        progress::log("  Found " + std::to_string(query_locs.size()) + " query spacers");
        progress::step_done("Step 1", std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count());

        // Step 2: Extract all reference spacer locations (streaming, zero-copy)
        step_start_time = std::chrono::high_resolution_clock::now();
        progress::step_start("Step 2", "Extracting reference spacer locations");
        std::vector<SpacerLocation> all_target_locs;
        spacerExtractionLoc(genome, pam, alt_pam, spacer_len, upstream, all_target_locs, num_threads, dust_threshold);
        progress::log("  Found " + std::to_string(all_target_locs.size()) + " target spacers");
        progress::step_done("Step 2", std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count());

        const size_t query_batch_size = std::max<size_t>(1, choose_query_processing_batch_size(query_locs.size(), false));
        if (query_batch_size < query_locs.size()) {
            progress::log("  Low-memory query batching enabled (" +
                          std::to_string(query_batch_size) + " query spacers per batch)");
        }

        std::ofstream final_out_stream(output_file, std::ios::binary | std::ios::trunc);
        if (!final_out_stream.is_open()) {
            throw std::runtime_error("Could not open output file: " + output_file);
        }
        if (raw_output) {
            if (!html_output_file.empty()) {
                std::cerr << "Error: --html-output requires scored output, not --raw-output." << std::endl;
                return 1;
            }
            write_raw_results_header(final_out_stream);
        } else {
            write_scored_results_header(final_out_stream);
        }

        double step3_seconds = 0.0;
        double step4_seconds = 0.0;
        double step5_seconds = 0.0;
        double step6_seconds = 0.0;
        size_t total_discarded = 0;

        MismatchPreencodedTargets mismatch_preencoded;
        if (max_mismatch >= 0) {
            progress::step_start("Step 3", "Mismatch search (Location-based, Parallel Indexing)");
            step_start_time = std::chrono::high_resolution_clock::now();
            mismatch_preencoded = preencode_targets_for_mismatch(all_target_locs, genome, num_threads);
            step3_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
        }

        PreencodedTargets indel_preencoded;
        if (max_indel > 0) {
            progress::step_start("Step 4", "Indel search (Location-based)");
            progress::log("Pre-encoding " + std::to_string(all_target_locs.size()) + " target sequences for indel search...");
            step_start_time = std::chrono::high_resolution_clock::now();
            indel_preencoded = preencode_targets_for_indel(all_target_locs, genome, num_threads);
            step4_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
            progress::log("Target pre-encoding complete.");
            progress::log("  Using " + std::to_string(num_threads) +
                          " threads for indel search (query batch size 16)");
        }

        for (size_t batch_start = 0; batch_start < query_locs.size(); batch_start += query_batch_size) {
            const size_t batch_count = std::min(query_batch_size, query_locs.size() - batch_start);
            std::vector<SpacerLocation> batch_query_locs(query_locs.begin() + static_cast<std::ptrdiff_t>(batch_start),
                                                         query_locs.begin() + static_cast<std::ptrdiff_t>(batch_start + batch_count));
            std::vector<SearchHit> batch_mismatch_hits;
            std::vector<SearchHit> batch_indel_hits;

            if (max_mismatch >= 0) {
                step_start_time = std::chrono::high_resolution_clock::now();
                batch_mismatch_hits = process_mismatch_search_loc_collect_preencoded(batch_query_locs,
                                                                                     all_target_locs,
                                                                                     mismatch_preencoded,
                                                                                     query_genome,
                                                                                     genome,
                                                                                     max_mismatch,
                                                                                     num_threads);
                step3_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
            }

            if (max_indel > 0) {
                step_start_time = std::chrono::high_resolution_clock::now();
                batch_indel_hits = process_indel_search_loc_collect_preencoded(batch_query_locs,
                                                                               all_target_locs,
                                                                               query_genome,
                                                                               genome,
                                                                               indel_preencoded,
                                                                               max_indel,
                                                                               num_threads,
                                                                               spacer_len,
                                                                               batch_start,
                                                                               query_locs.size());
                step4_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
            }

            step_start_time = std::chrono::high_resolution_clock::now();
            std::vector<CombinedResult> combined_results =
                combine_results_direct(batch_mismatch_hits, batch_indel_hits, max_mismatch, max_indel);
            step5_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();

            if (raw_output) {
                append_raw_results_tsv(final_out_stream, combined_results);
            } else {
                step_start_time = std::chrono::high_resolution_clock::now();
                size_t discarded_count = 0;
                std::vector<ScoredResult> scored_results = score_results_cpp(combined_results,
                                                                             genome,
                                                                             max_indel,
                                                                             max_mismatch,
                                                                             upstream ? -1 : 1,
                                                                             num_threads,
                                                                             discarded_count);
                total_discarded += discarded_count;
                append_scored_results_tsv(final_out_stream, scored_results);
                step6_seconds += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count();
            }
        }

        final_out_stream.flush();

        if (max_mismatch >= 0) {
            progress::step_done("Step 3", step3_seconds);
        }
        if (max_indel > 0) {
            progress::step_done("Step 4", step4_seconds);
        }

        progress::step_start("Step 5", "Combining results");
        progress::step_done("Step 5", step5_seconds);

        if (!raw_output) {
            progress::step_start("Step 6", "Post-processing and scoring");
            progress::step_done("Step 6", step6_seconds);
            if (total_discarded > 0) {
                progress::log("  Note: discarded " + std::to_string(total_discarded) +
                              " records without a valid alignment path during post-processing.");
            }

            if (!html_output_file.empty()) {
                step_start_time = std::chrono::high_resolution_clock::now();
                progress::step_start("Step 7", "HTML report generation");
                std::string report_error;
                if (!write_html_report_from_tsv(output_file, html_output_file, report_error)) {
                    progress::error(report_error);
                    return 1;
                }
                progress::log("  HTML report written to " + html_output_file);
                progress::step_done("Step 7",
                                    std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - step_start_time).count());
            }
        }

    } catch (const std::exception& e) {
        progress::error(e.what());
        return 1;
    }

    progress::complete(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - total_start_time).count(),
                       output_file,
                       html_output_file);
    return 0;
}
