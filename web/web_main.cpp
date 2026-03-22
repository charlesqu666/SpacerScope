#include "web_assets.hpp"
#include "web_job.hpp"

#include "third_party/httplib/httplib.h"
#include "third_party/nlohmann/json.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <cctype>
#include <string_view>

#ifdef _WIN32
#include <Windows.h>
#include <shellapi.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <sys/socket.h>
#endif

#ifndef SPACERSCOPE_VERSION
#define SPACERSCOPE_VERSION "dev"
#endif

namespace fs = std::filesystem;

namespace {

void set_json(httplib::Response& res, const nlohmann::json& payload, int status = 200) {
    res.status = status;
    res.set_content(payload.dump(2), "application/json; charset=utf-8");
}

std::string get_field(const httplib::Request& req, const std::string& key) {
    if (req.has_param(key)) {
        return req.get_param_value(key);
    }
    if (req.has_file(key)) {
        return req.get_file_value(key).content;
    }
    const std::string content_type = req.get_header_value("Content-Type");
    if (content_type.find("multipart/form-data") != std::string::npos && !req.body.empty()) {
        const std::string marker = "name=\"" + key + "\"";
        std::size_t pos = req.body.find(marker);
        while (pos != std::string::npos) {
            const std::size_t header_start = req.body.rfind("\r\n", pos);
            const std::size_t line_end = req.body.find("\r\n", pos);
            if (line_end == std::string::npos) break;
            const std::size_t data_start = req.body.find("\r\n\r\n", line_end);
            if (data_start == std::string::npos) break;
            const std::string_view header_block(req.body.data() + (header_start == std::string::npos ? 0 : header_start + 2),
                                                data_start - (header_start == std::string::npos ? 0 : header_start + 2));
            if (header_block.find("filename=") == std::string_view::npos) {
                const std::size_t value_start = data_start + 4;
                const std::size_t value_end = req.body.find("\r\n", value_start);
                if (value_end == std::string::npos) break;
                return req.body.substr(value_start, value_end - value_start);
            }
            pos = req.body.find(marker, data_start + 4);
        }
    }
    return {};
}

int get_int_field(const httplib::Request& req, const std::string& key, int fallback) {
    const std::string value = get_field(req, key);
    if (value.empty()) return fallback;
    return std::stoi(value);
}

bool get_bool_field(const httplib::Request& req, const std::string& key, bool fallback) {
    const std::string value = get_field(req, key);
    if (value.empty()) return fallback;
    std::string normalized = value;
    for (char& ch : normalized) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on";
}

std::string sanitize_filename(const std::string& name, const std::string& fallback) {
    if (name.empty()) return fallback;
    std::string safe = name;
    for (char& ch : safe) {
        const bool ok = std::isalnum(static_cast<unsigned char>(ch)) || ch == '.' || ch == '_' || ch == '-';
        if (!ok) ch = '_';
    }
    return safe.empty() ? fallback : safe;
}

fs::path path_from_web_value(const std::string& text) {
#ifdef _WIN32
    if (text.empty()) return fs::path();
    return fs::u8path(text);
#else
    return fs::path(text);
#endif
}

std::string path_for_message(const fs::path& path) {
#ifdef _WIN32
    return path.u8string();
#else
    return path.string();
#endif
}

fs::path executable_dir() {
#ifdef _WIN32
    wchar_t buffer[MAX_PATH];
    DWORD len = GetModuleFileNameW(nullptr, buffer, MAX_PATH);
    if (len == 0 || len == MAX_PATH) return fs::current_path();
    return fs::path(std::wstring(buffer, len)).parent_path();
#elif defined(__APPLE__)
    char pathbuf[1024];
    uint32_t bufsize = sizeof(pathbuf);
    if (_NSGetExecutablePath(pathbuf, &bufsize) == 0) {
        std::error_code ec;
        fs::path p = fs::canonical(path, ec);
        if (!ec) return p.parent_path();
        return fs::path(pathbuf).parent_path();
    }
    return fs::current_path();
#else
    std::error_code ec;
    fs::path p = fs::read_symlink("/proc/self/exe", ec);
    if (ec) return fs::current_path();
    return p.parent_path();
#endif
}

#ifdef _WIN32
std::wstring widen_best_effort_local(const std::string& text) {
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

std::wstring browser_url_for_host(const std::string& bind_host, int port) {
    std::string host = bind_host;
    if (host == "0.0.0.0" || host == "::" || host == "[::]") {
        host = "127.0.0.1";
    }
    std::ostringstream oss;
    oss << "http://" << host << ":" << port;
    return widen_best_effort_local(oss.str());
}

void open_browser_async(const std::string& bind_host, int port) {
    const std::wstring url = browser_url_for_host(bind_host, port);
    std::thread([url]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        ShellExecuteW(nullptr, L"open", url.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
    }).detach();
}

#elif defined(__APPLE__)
void open_browser_async(const std::string& bind_host, int port) {
    std::string host = bind_host;
    if (host == "0.0.0.0" || host == "::" || host == "[::]") {
        host = "127.0.0.1";
    }
    std::string url = "http://" + host + ":" + std::to_string(port);
    std::thread([url]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        std::string cmd = "open " + url;
        system(cmd.c_str());
    }).detach();
}
#endif

void print_usage() {
    std::cout
        << "Usage:\n"
        << "  spacerscope-web [options]\n\n"
        << "Options:\n"
        << "  --bind <host>            Bind host (default: 127.0.0.1)\n"
        << "  --port <int>             Bind port (default: 8787)\n"
        << "  --workspace-root <dir>   Workspace root for uploads and artifacts\n"
        << "  --spacerscope-bin <path> Path to spacerscope executable\n"
        << "  -h, --help               Show help\n"
        << "  -V, --version            Show version\n\n"
        << "Notes:\n"
        << "  - Default bind is localhost only.\n"
        << "  - Binding to 0.0.0.0 enables LAN access with no authentication in v1.\n";
}

std::string file_to_string(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

bool save_uploaded_file(const httplib::MultipartFormData& file,
                        const fs::path& dir,
                        const std::string& fallback_name,
                        fs::path& saved_path,
                        std::string& error_message) {
    if (file.content.empty()) {
        error_message = "Uploaded file is empty: " + fallback_name;
        return false;
    }
    fs::create_directories(dir);
    saved_path = dir / sanitize_filename(file.filename, fallback_name);
    std::ofstream out(saved_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        error_message = "Unable to write uploaded file: " + saved_path.string();
        return false;
    }
    out.write(file.content.data(), static_cast<std::streamsize>(file.content.size()));
    return true;
}

bool parse_run_request(const httplib::Request& req,
                       const fs::path& workspace_root,
                       WebRunRequest& out,
                       std::string& error_message) {
    std::string input_mode = get_field(req, "input_mode");
    for (char& ch : input_mode) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (input_mode.empty()) input_mode = "path";
    if (input_mode != "path" && input_mode != "upload") {
        error_message = "Invalid input mode.";
        return false;
    }

    out.mode = get_field(req, "mode");
    out.seq_type = get_field(req, "seq_type");
    if (out.seq_type.empty()) out.seq_type = "gene";
    out.pam = get_field(req, "pam");
    if (out.pam.empty()) out.pam = "NGG";
    out.alt_pam = get_field(req, "alt_pam");
    if (out.alt_pam.empty()) out.alt_pam = "NAG";
    out.direction = get_field(req, "direction");
    if (out.direction.empty()) out.direction = "up";
    out.gene_id = get_field(req, "gene_id");
    out.spacer_len = get_int_field(req, "spacer_len", 20);
    out.mismatch = get_int_field(req, "mismatch", 4);
    out.indel = get_int_field(req, "indel", 2);
    out.dust_threshold = get_int_field(req, "dust_threshold", 12);
    out.threads = get_int_field(req, "threads", 1);
    out.raw_output = get_bool_field(req, "raw_output", false);
    out.generate_html = get_bool_field(req, "generate_html", true);

    if (out.mode != "fasta" && out.mode != "fasta-cut" &&
        out.mode != "anno" && out.mode != "anno-cut") {
        error_message = "Invalid mode.";
        return false;
    }

    const fs::path upload_dir = workspace_root / "_uploads" / std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count());

    if (input_mode == "upload") {
        if (!req.has_file("ref_upload")) {
            error_message = "Reference FASTA upload is required in upload mode.";
            return false;
        }
        if (!save_uploaded_file(req.get_file_value("ref_upload"), upload_dir, "reference.fa", out.ref_path, error_message)) {
            return false;
        }
    } else {
        out.ref_path = path_from_web_value(get_field(req, "ref_path"));
    }

    if (input_mode == "upload") {
        if (out.mode == "fasta" || out.mode == "fasta-cut") {
            if (!req.has_file("query_upload")) {
                error_message = "Query FASTA upload is required for fasta modes in upload mode.";
                return false;
            }
            if (!save_uploaded_file(req.get_file_value("query_upload"), upload_dir, "query.fa", out.query_path, error_message)) {
                return false;
            }
        }
    } else {
        out.query_path = path_from_web_value(get_field(req, "query_path"));
    }

    if (input_mode == "upload") {
        if (out.mode == "anno" || out.mode == "anno-cut") {
            if (!req.has_file("annotation_upload")) {
                error_message = "Annotation upload is required for anno modes in upload mode.";
                return false;
            }
            if (!save_uploaded_file(req.get_file_value("annotation_upload"), upload_dir, "annotation.gff3", out.annotation_path, error_message)) {
                return false;
            }
        }
    } else {
        out.annotation_path = path_from_web_value(get_field(req, "annotation_path"));
    }

    if (out.ref_path.empty()) {
        error_message = "Reference FASTA is required.";
        return false;
    }
    if (!fs::exists(out.ref_path)) {
        error_message = "Reference FASTA not found: " + path_for_message(out.ref_path);
        return false;
    }

    if (out.mode == "fasta" || out.mode == "fasta-cut") {
        if (out.query_path.empty()) {
            error_message = "Query FASTA is required for fasta modes.";
            return false;
        }
        if (!fs::exists(out.query_path)) {
            error_message = "Query FASTA not found: " + path_for_message(out.query_path);
            return false;
        }
    } else {
        if (out.annotation_path.empty()) {
            error_message = "Annotation file is required for anno modes.";
            return false;
        }
        if (!fs::exists(out.annotation_path)) {
            error_message = "Annotation file not found: " + path_for_message(out.annotation_path);
            return false;
        }
        if (out.gene_id.empty()) {
            error_message = "Gene ID is required for anno modes.";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string bind_host = "127.0.0.1";
    int port = 8787;
    fs::path workspace_root = fs::temp_directory_path() / "spacerscope-web";
    fs::path spacerscope_bin = executable_dir() /
#ifdef _WIN32
        "spacerscope.exe";
#else
        "spacerscope";
#endif

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& opt) -> bool {
            if (i + 1 >= argc) {
                std::cerr << "Error: option '" << opt << "' requires a value." << std::endl;
                return false;
            }
            return true;
        };
        if (arg == "-h" || arg == "--help") {
            print_usage();
            return 0;
        }
        if (arg == "-V" || arg == "--version") {
            std::cout << "spacerscope-web " << SPACERSCOPE_VERSION << std::endl;
            return 0;
        }
        if (arg == "--bind") {
            if (!require_value(arg)) return 1;
            bind_host = argv[++i];
        } else if (arg == "--port") {
            if (!require_value(arg)) return 1;
            port = std::stoi(argv[++i]);
        } else if (arg == "--workspace-root") {
            if (!require_value(arg)) return 1;
            workspace_root = argv[++i];
        } else if (arg == "--spacerscope-bin") {
            if (!require_value(arg)) return 1;
            spacerscope_bin = argv[++i];
        } else {
            std::cerr << "Error: unknown argument '" << arg << "'." << std::endl;
            return 1;
        }
    }

    if (!fs::exists(spacerscope_bin)) {
        std::cerr << "Error: spacerscope binary not found: " << spacerscope_bin << std::endl;
        return 1;
    }

    if (port <= 0 || port > 65535) {
        std::cerr << "Error: --port must be between 1 and 65535." << std::endl;
        return 1;
    }

    fs::create_directories(workspace_root);

    WebServerInfo info;
    info.bind_host = bind_host;
    info.port = port;
    info.version = SPACERSCOPE_VERSION;
    info.workspace_root = workspace_root;
    info.spacerscope_bin = spacerscope_bin;
    JobManager manager(info);

    httplib::Server server;
#ifndef _WIN32
    server.set_socket_options([](int sock) {
        int yes = 1;
        setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));
    });
#endif

    server.Get("/", [&](const httplib::Request&, httplib::Response& res) {
        res.set_content(web_assets::index_html(), "text/html; charset=utf-8");
    });

    server.Get("/api/info", [&](const httplib::Request&, httplib::Response& res) {
        set_json(res, manager.info_json());
    });

    server.Get("/api/status", [&](const httplib::Request&, httplib::Response& res) {
        set_json(res, manager.status_json());
    });

    server.Post("/api/run", [&](const httplib::Request& req, httplib::Response& res) {
        WebRunRequest run_request;
        std::string error_message;
        if (!parse_run_request(req, workspace_root, run_request, error_message)) {
            set_json(res, {{"error", error_message}}, 400);
            return;
        }
        if (!manager.start(run_request, error_message)) {
            set_json(res, {{"error", error_message}}, 409);
            return;
        }
        nlohmann::json payload = manager.status_json();
        payload["accepted"] = true;
        set_json(res, payload, 202);
    });

    server.Post("/api/cancel", [&](const httplib::Request&, httplib::Response& res) {
        std::string error_message;
        if (!manager.cancel(error_message)) {
            set_json(res, {{"error", error_message}}, 409);
            return;
        }
        set_json(res, {{"accepted", true}});
    });

    server.Get("/api/events", [&](const httplib::Request&, httplib::Response& res) {
        res.set_header("Cache-Control", "no-cache");
        res.set_header("Connection", "keep-alive");
        res.set_header("X-Accel-Buffering", "no");

        const auto initial_status = manager.status_json();
        const auto initial_cursor = manager.event_count();

        res.set_chunked_content_provider(
            "text/event-stream; charset=utf-8",
            [initial_status, initial_cursor, &manager, sent_initial = false, cursor = initial_cursor]
            (size_t, httplib::DataSink& sink) mutable {
                auto send_payload = [&](const std::string& payload) -> bool {
                    const std::string frame = "data: " + payload + "\n\n";
                    return sink.write(frame.data(), frame.size());
                };

                if (!sent_initial) {
                    nlohmann::json payload = initial_status;
                    payload["event"] = "status";
                    sent_initial = true;
                    return send_payload(payload.dump());
                }

                if (!manager.wait_for_events(cursor, std::chrono::milliseconds(750))) {
                    static const std::string heartbeat = ": keep-alive\n\n";
                    return sink.write(heartbeat.data(), heartbeat.size());
                }

                auto events = manager.events_since(cursor);
                for (const auto& event_payload : events) {
                    if (!send_payload(event_payload)) return false;
                }
                return true;
            });
    });

    server.Get(R"(/api/artifacts/(tsv|html|stdout|stderr))", [&](const httplib::Request& req, httplib::Response& res) {
        const std::string kind = req.matches[1];
        fs::path path;
        std::string content_type;
        if (!manager.artifact_path(kind, path, content_type)) {
            set_json(res, {{"error", "Artifact not available."}}, 404);
            return;
        }
        res.set_content(file_to_string(path), content_type);
    });

    if (!server.bind_to_port(bind_host, port)) {
        std::cerr << "Error: failed to bind to " << bind_host << ":" << port
                  << ". The port may already be in use. Try --port <other>." << std::endl;
        return 1;
    }
    std::cout << "spacerscope-web " << SPACERSCOPE_VERSION
              << " listening on http://" << bind_host << ":" << port << std::endl;
    if (bind_host != "127.0.0.1" && bind_host != "::1" && bind_host != "localhost") {
        std::cout << "Warning: LAN mode has no authentication in v1. Use only on a trusted network." << std::endl;
    }
#if defined(_WIN32) || defined(__APPLE__)
    open_browser_async(bind_host, port);
#endif
    if (!server.listen_after_bind()) {
        std::cerr << "Error: server stopped unexpectedly after binding "
                  << bind_host << ":" << port << std::endl;
        return 1;
    }
    return 0;
}
