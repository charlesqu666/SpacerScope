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
#include "web_job.hpp"

#include "third_party/nlohmann/json.hpp"

#include <atomic>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#include <fcntl.h>
#include <io.h>
#else
#include <csignal>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace {

std::string state_to_string(JobManager::State state) {
    switch (state) {
        case JobManager::State::Idle: return "idle";
        case JobManager::State::Starting: return "starting";
        case JobManager::State::Running: return "running";
        case JobManager::State::Completed: return "completed";
        case JobManager::State::Failed: return "failed";
        case JobManager::State::Canceled: return "canceled";
    }
    return "unknown";
}

bool truthy(std::string value) {
    for (char& ch : value) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return value == "1" || value == "true" || value == "yes" || value == "on";
}

std::string sanitize_filename(std::string name, const std::string& fallback) {
    if (name.empty()) return fallback;
    for (char& ch : name) {
        const bool ok = std::isalnum(static_cast<unsigned char>(ch)) || ch == '.' || ch == '_' || ch == '-';
        if (!ok) ch = '_';
    }
    return name.empty() ? fallback : name;
}

std::string current_timestamp_label() {
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    return std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(now).count());
}

std::string join_command_for_log(const std::vector<std::string>& args) {
    std::ostringstream oss;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) oss << ' ';
        oss << args[i];
    }
    return oss.str();
}

#ifdef _WIN32
std::wstring widen_best_effort(const std::string& text) {
    if (text.empty()) return std::wstring();
    auto convert = [&](UINT code_page) -> std::wstring {
        const int required = MultiByteToWideChar(code_page, 0, text.c_str(), -1, nullptr, 0);
        if (required <= 0) return std::wstring();
        std::wstring wide(static_cast<size_t>(required - 1), L'\0');
        const int written = MultiByteToWideChar(code_page, 0, text.c_str(), -1, wide.data(), required);
        if (written <= 0) return std::wstring();
        return wide;
    };
    std::wstring wide = convert(CP_UTF8);
    if (!wide.empty()) return wide;
    return convert(CP_ACP);
}

std::string path_arg_string(const fs::path& path) {
    return path.u8string();
}

std::string quote_windows_arg(const std::string& arg) {
    if (arg.find_first_of(" \t\"") == std::string::npos) return arg;
    std::string quoted = "\"";
    size_t backslashes = 0;
    for (char ch : arg) {
        if (ch == '\\') {
            ++backslashes;
            continue;
        }
        if (ch == '"') {
            quoted.append(backslashes * 2 + 1, '\\');
            quoted.push_back('"');
            backslashes = 0;
            continue;
        }
        quoted.append(backslashes, '\\');
        backslashes = 0;
        quoted.push_back(ch);
    }
    quoted.append(backslashes * 2, '\\');
    quoted.push_back('"');
    return quoted;
}

std::wstring quote_windows_arg(const std::wstring& arg) {
    if (arg.find_first_of(L" \t\"") == std::wstring::npos) return arg;
    std::wstring quoted = L"\"";
    size_t backslashes = 0;
    for (wchar_t ch : arg) {
        if (ch == L'\\') {
            ++backslashes;
            continue;
        }
        if (ch == L'"') {
            quoted.append(backslashes * 2 + 1, L'\\');
            quoted.push_back(L'"');
            backslashes = 0;
            continue;
        }
        quoted.append(backslashes, L'\\');
        backslashes = 0;
        quoted.push_back(ch);
    }
    quoted.append(backslashes * 2, L'\\');
    quoted.push_back(L'"');
    return quoted;
}
#else
std::string path_arg_string(const fs::path& path) {
    return path.string();
}
#endif

}  // namespace

struct JobManager::ProcessHandle {
    std::atomic<bool> cancel_requested{false};
    std::thread stdout_thread;
    std::thread stderr_thread;
    std::thread wait_thread;

#ifdef _WIN32
    HANDLE process = nullptr;
    HANDLE job = nullptr;
    HANDLE stdout_read = nullptr;
    HANDLE stderr_read = nullptr;
#else
    pid_t pid = -1;
    int stdout_fd = -1;
    int stderr_fd = -1;
#endif

    ~ProcessHandle() {
        if (stdout_thread.joinable()) stdout_thread.join();
        if (stderr_thread.joinable()) stderr_thread.join();
        if (wait_thread.joinable()) wait_thread.join();
#ifdef _WIN32
        if (stdout_read) CloseHandle(stdout_read);
        if (stderr_read) CloseHandle(stderr_read);
        if (process) CloseHandle(process);
        if (job) CloseHandle(job);
#else
        if (stdout_fd >= 0) close(stdout_fd);
        if (stderr_fd >= 0) close(stderr_fd);
#endif
    }
};

JobManager::JobManager(WebServerInfo info) : info_(std::move(info)) {
    fs::create_directories(info_.workspace_root);
}

JobManager::~JobManager() {
    std::string ignored;
    cancel(ignored);
    std::error_code ec;
    fs::remove_all(info_.workspace_root, ec);
}

nlohmann::json JobManager::info_json() const {
    return {
        {"product", "spacerscope-web"},
        {"version", info_.version},
        {"bind", info_.bind_host},
        {"port", info_.port},
        {"workspace_root", path_arg_string(info_.workspace_root)},
        {"spacerscope_bin", path_arg_string(info_.spacerscope_bin)},
        {"defaults",
            {
                {"pam", "NGG"},
                {"alt_pam", "NAG"},
                {"spacer_len", 20},
                {"direction", "up"},
                {"mismatch", 4},
                {"indel", 2},
                {"dust_threshold", 12},
                {"threads", (std::max)(1u, std::thread::hardware_concurrency())},
            }},
        {"modes", {"fasta", "fasta-cut", "anno", "anno-cut"}},
        {"lan_warning", "LAN mode in v1 has no authentication; only use it on a trusted network."},
    };
}

nlohmann::json JobManager::build_status_locked() const {
    nlohmann::json artifacts = {
        {"tsv", !current_tsv_.empty() && fs::exists(current_tsv_)},
        {"html", !current_html_.empty() && fs::exists(current_html_)},
        {"stdout", !current_stdout_.empty() && fs::exists(current_stdout_)},
        {"stderr", !current_stderr_.empty() && fs::exists(current_stderr_)},
    };

    nlohmann::json payload = {
        {"run_id", run_id_},
        {"state", state_to_string(state_)},
        {"mode", current_mode_},
        {"workspace", current_workspace_.empty() ? std::string() : path_arg_string(current_workspace_)},
        {"last_error", last_error_},
        {"can_cancel", state_ == State::Starting || state_ == State::Running},
        {"artifacts", artifacts},
    };
    if (!last_memory_estimate_.is_null()) {
        payload["memory_estimate"] = last_memory_estimate_;
    }
    return payload;
}

nlohmann::json JobManager::status_json() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return build_status_locked();
}

void JobManager::publish_event(nlohmann::json payload) {
    std::lock_guard<std::mutex> lock(mutex_);
    publish_event_locked(std::move(payload));
}

void JobManager::publish_event_locked(nlohmann::json payload) {
    payload["run_id"] = run_id_;
    if (payload.value("event", "") == "memory_estimate") {
        last_memory_estimate_ = payload;
    } else if (payload.value("event", "") == "error") {
        last_error_ = payload.value("message", "");
    }
    events_.push_back(payload.dump());
    cv_.notify_all();
}

void JobManager::publish_status_locked() {
    nlohmann::json payload = build_status_locked();
    payload["event"] = "status";
    payload["run_id"] = run_id_;
    events_.push_back(payload.dump());
    cv_.notify_all();
}

void JobManager::set_state_locked(State new_state, const std::string& error_message) {
    state_ = new_state;
    if (!error_message.empty()) last_error_ = error_message;
    publish_status_locked();
}

void JobManager::clear_previous_workspace_locked() {
    std::error_code ec;
    if (!current_workspace_.empty()) {
        fs::remove_all(current_workspace_, ec);
    }
    current_workspace_.clear();
    current_tsv_.clear();
    current_html_.clear();
    current_stdout_.clear();
    current_stderr_.clear();
    last_error_.clear();
    last_memory_estimate_ = nlohmann::json();
}

fs::path JobManager::create_run_workspace_locked() {
    clear_previous_workspace_locked();
    const fs::path dir = info_.workspace_root / ("run-" + current_timestamp_label());
    fs::create_directories(dir);
    current_workspace_ = dir;
    current_tsv_ = dir / "result.tsv";
    current_html_ = dir / "result.html";
    current_stdout_ = dir / "stdout.log";
    current_stderr_ = dir / "stderr.log";
    return dir;
}

bool JobManager::artifact_path(const std::string& kind,
                               fs::path& path,
                               std::string& content_type) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (kind == "tsv" && !current_tsv_.empty() && fs::exists(current_tsv_)) {
        path = current_tsv_;
        content_type = "text/tab-separated-values; charset=utf-8";
        return true;
    }
    if (kind == "html" && !current_html_.empty() && fs::exists(current_html_)) {
        path = current_html_;
        content_type = "text/html; charset=utf-8";
        return true;
    }
    if (kind == "stdout" && !current_stdout_.empty() && fs::exists(current_stdout_)) {
        path = current_stdout_;
        content_type = "text/plain; charset=utf-8";
        return true;
    }
    if (kind == "stderr" && !current_stderr_.empty() && fs::exists(current_stderr_)) {
        path = current_stderr_;
        content_type = "text/plain; charset=utf-8";
        return true;
    }
    return false;
}

std::vector<std::string> JobManager::events_since(std::size_t& cursor) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (cursor > events_.size()) cursor = 0;
    std::vector<std::string> out;
    out.reserve(events_.size() - cursor);
    for (std::size_t i = cursor; i < events_.size(); ++i) out.push_back(events_[i]);
    cursor = events_.size();
    return out;
}

bool JobManager::wait_for_events(std::size_t cursor, std::chrono::milliseconds timeout) const {
    std::unique_lock<std::mutex> lock(mutex_);
    return cv_.wait_for(lock, timeout, [&] { return cursor > events_.size() ? true : cursor < events_.size(); });
}

std::size_t JobManager::event_count() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return events_.size();
}

static std::vector<std::string> build_spacerscope_args(const fs::path& bin,
                                                       const WebRunRequest& request,
                                                       const fs::path& output_tsv,
                                                       const fs::path& output_html) {
    std::vector<std::string> args;
    args.push_back(path_arg_string(bin));
    args.push_back(request.mode);
    args.push_back("--ref");
    args.push_back(path_arg_string(request.ref_path));
    args.push_back("-o");
    args.push_back(path_arg_string(output_tsv));
    args.push_back("--pam");
    args.push_back(request.pam);
    args.push_back("--alt-pam");
    args.push_back(request.alt_pam);
    args.push_back("--spacer-len");
    args.push_back(std::to_string(request.spacer_len));
    args.push_back("--direction");
    args.push_back(request.direction);
    args.push_back("--mismatch");
    args.push_back(std::to_string(request.mismatch));
    args.push_back("--indel");
    args.push_back(std::to_string(request.indel));
    args.push_back("--dust-threshold");
    args.push_back(std::to_string(request.dust_threshold));
    args.push_back("--threads");
    args.push_back(std::to_string(request.threads));
    args.push_back("--progress-format");
    args.push_back("jsonl");
    args.push_back("--yes");

    if (request.raw_output) {
        args.push_back("--raw-output");
    } else if (request.generate_html) {
        args.push_back("--html-output");
        args.push_back(path_arg_string(output_html));
    }

    if (request.mode == "fasta" || request.mode == "fasta-cut") {
        args.push_back("--query");
        args.push_back(path_arg_string(request.query_path));
    } else {
        args.push_back("--annotation");
        args.push_back(path_arg_string(request.annotation_path));
        args.push_back("--geneID");
        args.push_back(request.gene_id);
        args.push_back("--seq-type");
        args.push_back(request.seq_type);
    }
    return args;
}

static bool save_bytes(const fs::path& path, const std::string& data, std::string& error_message) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        error_message = "Unable to write uploaded file: " + path.string();
        return false;
    }
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    return true;
}

bool JobManager::start(const WebRunRequest& request, std::string& error_message) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (state_ == State::Starting || state_ == State::Running) {
        error_message = "A task is already running.";
        return false;
    }
    if (request.ref_path.empty()) {
        error_message = "Reference FASTA is required.";
        return false;
    }
    if ((request.mode == "fasta" || request.mode == "fasta-cut") && request.query_path.empty()) {
        error_message = "Query FASTA is required for fasta modes.";
        return false;
    }
    if ((request.mode == "anno" || request.mode == "anno-cut") &&
        (request.annotation_path.empty() || request.gene_id.empty())) {
        error_message = "Annotation file and gene ID are required for anno modes.";
        return false;
    }
    if (request.raw_output && request.generate_html) {
        error_message = "Raw output and HTML generation cannot be enabled together.";
        return false;
    }

    process_.reset();
    ++run_id_;
    events_.clear();
    current_mode_ = request.mode;
    create_run_workspace_locked();
    set_state_locked(State::Starting);
    publish_event_locked({
        {"event", "log"},
        {"message", "[web] launching spacerscope subprocess"},
    });
    publish_event_locked({
        {"event", "log"},
        {"message", "[web] command: " + join_command_for_log(
            build_spacerscope_args(info_.spacerscope_bin, request, current_tsv_, current_html_))},
    });
    if (!launch_process_locked(request, error_message)) {
        set_state_locked(State::Failed, error_message);
        publish_event_locked({{"event", "error"}, {"message", error_message}});
        return false;
    }
    return true;
}

bool JobManager::cancel(std::string& error_message) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!process_ || (state_ != State::Starting && state_ != State::Running)) {
        error_message = "No running task to cancel.";
        return false;
    }
    process_->cancel_requested = true;
    auto* process = process_.get();
    lock.unlock();

#ifdef _WIN32
    if (process->job) {
        TerminateJobObject(process->job, 1);
    }
#else
    if (process->pid > 0) {
        killpg(process->pid, SIGTERM);
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        killpg(process->pid, SIGKILL);
    }
#endif
    publish_event({{"event", "log"}, {"message", "[web] cancel requested"}});
    return true;
}

bool JobManager::launch_process_locked(const WebRunRequest& request, std::string& error_message) {
    auto process = std::make_unique<ProcessHandle>();
    const std::vector<std::string> args = build_spacerscope_args(info_.spacerscope_bin, request, current_tsv_, current_html_);

    auto handle_stdout_line = [this](const std::string& line) {
        {
            std::ofstream out(current_stdout_, std::ios::app | std::ios::binary);
            out << line << '\n';
        }
        auto parsed = nlohmann::json::parse(line, nullptr, false);
        if (!parsed.is_discarded() && parsed.is_object() && parsed.contains("event")) {
            publish_event(parsed);
        } else {
            publish_event({{"event", "log"}, {"message", line}});
        }
    };
    auto handle_stderr_line = [this](const std::string& line) {
        {
            std::ofstream out(current_stderr_, std::ios::app | std::ios::binary);
            out << line << '\n';
        }
        publish_event({{"event", "stderr"}, {"message", line}});
    };

#ifdef _WIN32
    SECURITY_ATTRIBUTES sa{};
    sa.nLength = sizeof(sa);
    sa.bInheritHandle = TRUE;

    HANDLE stdout_read = nullptr, stdout_write = nullptr;
    HANDLE stderr_read = nullptr, stderr_write = nullptr;
    if (!CreatePipe(&stdout_read, &stdout_write, &sa, 0) ||
        !CreatePipe(&stderr_read, &stderr_write, &sa, 0)) {
        error_message = "CreatePipe failed.";
        return false;
    }
    SetHandleInformation(stdout_read, HANDLE_FLAG_INHERIT, 0);
    SetHandleInformation(stderr_read, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOW si{};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = stdout_write;
    si.hStdError = stderr_write;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);

    PROCESS_INFORMATION pi{};
    std::wstring command_line;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) command_line.push_back(L' ');
        command_line += quote_windows_arg(widen_best_effort(args[i]));
    }

    HANDLE job = CreateJobObjectA(nullptr, nullptr);
    if (!job) {
        error_message = "CreateJobObject failed.";
        CloseHandle(stdout_read);
        CloseHandle(stdout_write);
        CloseHandle(stderr_read);
        CloseHandle(stderr_write);
        return false;
    }
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION limit_info{};
    limit_info.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    SetInformationJobObject(job, JobObjectExtendedLimitInformation, &limit_info, sizeof(limit_info));

    std::wstring workspace_w = current_workspace_.wstring();
    std::vector<wchar_t> mutable_command(command_line.begin(), command_line.end());
    mutable_command.push_back(L'\0');

    const BOOL created = CreateProcessW(
        nullptr,
        mutable_command.data(),
        nullptr,
        nullptr,
        TRUE,
        CREATE_NO_WINDOW,
        nullptr,
        workspace_w.c_str(),
        &si,
        &pi);

    CloseHandle(stdout_write);
    CloseHandle(stderr_write);

    if (!created) {
        error_message = "CreateProcess failed.";
        CloseHandle(stdout_read);
        CloseHandle(stderr_read);
        CloseHandle(job);
        return false;
    }

    AssignProcessToJobObject(job, pi.hProcess);
    CloseHandle(pi.hThread);

    process->process = pi.hProcess;
    process->job = job;
    process->stdout_read = stdout_read;
    process->stderr_read = stderr_read;

    auto start_reader = [](HANDLE handle, auto on_line) {
        int fd = _open_osfhandle(reinterpret_cast<intptr_t>(handle), _O_RDONLY);
        FILE* file = _fdopen(fd, "rb");
        std::string line;
        char buffer[4096];
        while (file && std::fgets(buffer, static_cast<int>(sizeof(buffer)), file)) {
            line.assign(buffer);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            if (!line.empty()) on_line(line);
        }
        if (file) std::fclose(file);
    };

    process->stdout_thread = std::thread(start_reader, stdout_read, handle_stdout_line);
    process->stderr_thread = std::thread(start_reader, stderr_read, handle_stderr_line);
    process->stdout_read = nullptr;
    process->stderr_read = nullptr;
    process->wait_thread = std::thread([this, raw = process.get()] {
        WaitForSingleObject(raw->process, INFINITE);
        DWORD exit_code = 1;
        GetExitCodeProcess(raw->process, &exit_code);
        finalize_process(static_cast<int>(exit_code));
    });
#else
    int stdout_pipe[2];
    int stderr_pipe[2];
    if (pipe(stdout_pipe) != 0 || pipe(stderr_pipe) != 0) {
        error_message = "pipe() failed.";
        return false;
    }

    const pid_t pid = fork();
    if (pid < 0) {
        error_message = "fork() failed.";
        close(stdout_pipe[0]); close(stdout_pipe[1]);
        close(stderr_pipe[0]); close(stderr_pipe[1]);
        return false;
    }

    if (pid == 0) {
        setpgid(0, 0);
        dup2(stdout_pipe[1], STDOUT_FILENO);
        dup2(stderr_pipe[1], STDERR_FILENO);
        close(stdout_pipe[0]); close(stdout_pipe[1]);
        close(stderr_pipe[0]); close(stderr_pipe[1]);
        if (chdir(current_workspace_.string().c_str()) != 0) {
            _exit(127);
        }

        std::vector<char*> argv;
        argv.reserve(args.size() + 1);
        for (const auto& arg : args) argv.push_back(const_cast<char*>(arg.c_str()));
        argv.push_back(nullptr);
        execv(args[0].c_str(), argv.data());
        _exit(127);
    }

    setpgid(pid, pid);
    close(stdout_pipe[1]);
    close(stderr_pipe[1]);

    process->pid = pid;
    process->stdout_fd = stdout_pipe[0];
    process->stderr_fd = stderr_pipe[0];

    auto start_reader = [](int fd, auto on_line) {
        FILE* file = fdopen(fd, "r");
        char* line_ptr = nullptr;
        size_t cap = 0;
        while (file && getline(&line_ptr, &cap, file) != -1) {
            std::string line(line_ptr ? line_ptr : "");
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            if (!line.empty()) on_line(line);
        }
        if (line_ptr) std::free(line_ptr);
        if (file) fclose(file);
    };

    process->stdout_thread = std::thread(start_reader, stdout_pipe[0], handle_stdout_line);
    process->stderr_thread = std::thread(start_reader, stderr_pipe[0], handle_stderr_line);
    process->stdout_fd = -1;
    process->stderr_fd = -1;
    process->wait_thread = std::thread([this, raw = process.get()] {
        int status = 0;
        waitpid(raw->pid, &status, 0);
        int exit_code = 1;
        if (WIFEXITED(status)) exit_code = WEXITSTATUS(status);
        else if (WIFSIGNALED(status)) exit_code = 128 + WTERMSIG(status);
        finalize_process(exit_code);
    });
#endif

    process_ = std::move(process);
    set_state_locked(State::Running);
    return true;
}

void JobManager::finalize_process(int exit_code) {
    std::lock_guard<std::mutex> lock(mutex_);
    const bool canceled = process_ && process_->cancel_requested.load();
    if (canceled) {
        set_state_locked(State::Canceled);
        publish_event_locked({{"event", "log"}, {"message", "[web] task canceled"}});
        return;
    }
    if (exit_code == 0) {
        set_state_locked(State::Completed);
    } else {
        set_state_locked(State::Failed, "spacerscope exited with code " + std::to_string(exit_code));
        publish_event_locked({{"event", "error"}, {"message", "spacerscope exited with code " + std::to_string(exit_code)}});
    }
}
