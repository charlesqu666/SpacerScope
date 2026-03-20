#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace progress {

enum class Format {
    Text,
    Jsonl,
};

void set_format(Format format);
Format format();
bool is_text();
bool is_jsonl();

void log(const std::string& message);
void memory_estimate(std::uint64_t chrom_count,
                     double total_mbp,
                     double max_mbp,
                     double peak_mib,
                     double peak_gib,
                     const std::string& model);
void step_start(const std::string& step_id, const std::string& title);
void step_done(const std::string& step_id, double seconds);
void progress_update(const std::string& step_id,
                     std::size_t done,
                     std::size_t total,
                     const std::string& unit);
void complete(double total_seconds,
              const std::string& output_path = std::string(),
              const std::string& html_path = std::string());
void error(const std::string& message);

}  // namespace progress
