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
#include "progress.hpp"

#include "third_party/nlohmann/json.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>

namespace progress {
namespace {

std::mutex g_progress_mutex;
Format g_format = Format::Text;

std::string now_iso8601_utc() {
    using clock = std::chrono::system_clock;
    const auto now = clock::now();
    const std::time_t now_time = clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &now_time);
#else
    gmtime_r(&now_time, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

void emit_json_locked(nlohmann::json payload) {
    payload["ts"] = now_iso8601_utc();
    std::cout << payload.dump() << '\n';
    std::cout.flush();
}

void emit_text_locked(const std::string& message) {
    std::cout << message << std::endl;
}

void emit_event(const nlohmann::json& payload, const std::string& text_fallback) {
    std::lock_guard<std::mutex> lock(g_progress_mutex);
    if (g_format == Format::Jsonl) {
        emit_json_locked(payload);
    } else {
        emit_text_locked(text_fallback);
    }
}

}  // namespace

void set_format(Format format) {
    std::lock_guard<std::mutex> lock(g_progress_mutex);
    g_format = format;
}

Format format() {
    std::lock_guard<std::mutex> lock(g_progress_mutex);
    return g_format;
}

bool is_text() {
    return format() == Format::Text;
}

bool is_jsonl() {
    return format() == Format::Jsonl;
}

void log(const std::string& message) {
    emit_event({{"event", "log"}, {"message", message}}, message);
}

void memory_estimate(std::uint64_t chrom_count,
                     double total_mbp,
                     double max_mbp,
                     double peak_mib,
                     double peak_gib,
                     const std::string& model) {
    if (is_jsonl()) {
        emit_event({
                       {"event", "memory_estimate"},
                       {"chromosomes", chrom_count},
                       {"total_genome_mbp", total_mbp},
                       {"longest_chromosome_mbp", max_mbp},
                       {"peak_mib", peak_mib},
                       {"peak_gib", peak_gib},
                       {"model", model},
                   },
                   "");
        return;
    }

    std::ostringstream oss;
    oss << "\nMemory estimate:\n"
        << "  Reference chromosomes: " << chrom_count << "\n"
        << std::fixed << std::setprecision(2)
        << "  Total genome size: " << total_mbp << " Mb\n"
        << "  Longest chromosome: " << max_mbp << " Mb\n"
        << "  Estimated peak memory: " << peak_mib << " MiB (" << peak_gib << " GiB)\n"
        << "  Estimate model: " << model;
    std::lock_guard<std::mutex> lock(g_progress_mutex);
    std::cout << oss.str() << std::endl;
}

void step_start(const std::string& step_id, const std::string& title) {
    emit_event({{"event", "step_start"}, {"step", step_id}, {"title", title}},
               "\n" + step_id + ": " + title + "...");
}

void step_done(const std::string& step_id, double seconds) {
    std::ostringstream oss;
    oss << "  " << step_id << " took " << seconds << "s";
    emit_event({{"event", "step_done"}, {"step", step_id}, {"seconds", seconds}}, oss.str());
}

void progress_update(const std::string& step_id,
                     std::size_t done,
                     std::size_t total,
                     const std::string& unit) {
    const std::size_t percent = total == 0 ? 0 : (done * 100 / total);
    std::ostringstream oss;
    oss << "  " << step_id << " progress: " << done << "/" << total << " " << unit
        << " processed (" << percent << "%)";
    emit_event({
                   {"event", "progress"},
                   {"step", step_id},
                   {"done", done},
                   {"total", total},
                   {"unit", unit},
                   {"percent", percent},
               },
               oss.str());
}

void complete(double total_seconds, const std::string& output_path, const std::string& html_path) {
    if (is_jsonl()) {
        nlohmann::json payload{
            {"event", "complete"},
            {"total_seconds", total_seconds},
        };
        if (!output_path.empty()) payload["output"] = output_path;
        if (!html_path.empty()) payload["html"] = html_path;
        emit_event(payload, "");
        return;
    }

    std::ostringstream oss;
    oss << "\nComplete. Total Time: " << total_seconds << "s";
    emit_event({{"event", "complete"}, {"total_seconds", total_seconds}}, oss.str());
}

void error(const std::string& message) {
    std::lock_guard<std::mutex> lock(g_progress_mutex);
    if (g_format == Format::Jsonl) {
        emit_json_locked({{"event", "error"}, {"message", message}});
    } else {
        std::cerr << "Error: " << message << std::endl;
    }
}

}  // namespace progress
