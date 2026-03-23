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

#include "third_party/nlohmann/json.hpp"

#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <string>
#include <vector>

struct WebRunRequest {
    std::string mode;
    std::filesystem::path ref_path;
    std::filesystem::path query_path;
    std::filesystem::path annotation_path;
    std::string gene_id;
    std::string seq_type = "gene";
    std::string pam = "NGG";
    std::string alt_pam = "NAG";
    std::string direction = "up";
    int spacer_len = 20;
    int mismatch = 4;
    int indel = 2;
    int dust_threshold = 12;
    int threads = 1;
    bool raw_output = false;
    bool generate_html = true;
};

struct WebServerInfo {
    std::string bind_host;
    int port = 8787;
    std::string version;
    std::filesystem::path workspace_root;
    std::filesystem::path spacerscope_bin;
};

class JobManager {
public:
    enum class State {
        Idle,
        Starting,
        Running,
        Completed,
        Failed,
        Canceled,
    };

    explicit JobManager(WebServerInfo info);
    ~JobManager();

    nlohmann::json info_json() const;
    nlohmann::json status_json() const;

    bool start(const WebRunRequest& request, std::string& error_message);
    bool cancel(std::string& error_message);

    std::vector<std::string> events_since(std::size_t& cursor) const;
    bool wait_for_events(std::size_t cursor, std::chrono::milliseconds timeout) const;
    std::size_t event_count() const;

    bool artifact_path(const std::string& kind,
                       std::filesystem::path& path,
                       std::string& content_type) const;

private:
    void publish_event(nlohmann::json payload);
    void publish_event_locked(nlohmann::json payload);
    void publish_status_locked();
    nlohmann::json build_status_locked() const;
    void set_state_locked(State new_state, const std::string& error_message = std::string());
    void clear_previous_workspace_locked();
    std::filesystem::path create_run_workspace_locked();
    bool launch_process_locked(const WebRunRequest& request, std::string& error_message);
    void finalize_process(int exit_code);

    WebServerInfo info_;

    mutable std::mutex mutex_;
    mutable std::condition_variable cv_;

    State state_ = State::Idle;
    std::uint64_t run_id_ = 0;
    std::string current_mode_;
    std::string last_error_;
    nlohmann::json last_memory_estimate_;
    std::vector<std::string> events_;

    std::filesystem::path current_workspace_;
    std::filesystem::path current_tsv_;
    std::filesystem::path current_html_;
    std::filesystem::path current_stdout_;
    std::filesystem::path current_stderr_;

    struct ProcessHandle;
    std::unique_ptr<ProcessHandle> process_;
};
