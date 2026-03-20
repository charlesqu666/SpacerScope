#include "report_html.hpp"

#include "report_logo.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace {

constexpr int kHotThreshold = 10;

struct QuerySummary {
    std::string query_name;
    std::string safe_id;
    size_t total = 0;
    size_t exact = 0;
    size_t mismatch = 0;
    size_t indel = 0;
    size_t uniq = 0;
    double success_rate = 1.0;
    std::vector<size_t> row_indices;
};

std::string_view trim_ascii_whitespace(std::string_view text) {
    size_t start = 0;
    while (start < text.size()) {
        char ch = text[start];
        if (ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t') break;
        ++start;
    }
    size_t end = text.size();
    while (end > start) {
        char ch = text[end - 1];
        if (ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t') break;
        --end;
    }
    return text.substr(start, end - start);
}

std::string html_escape(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    for (char ch : text) {
        switch (ch) {
            case '&': out.append("&amp;"); break;
            case '<': out.append("&lt;"); break;
            case '>': out.append("&gt;"); break;
            case '"': out.append("&quot;"); break;
            case '\'': out.append("&#39;"); break;
            default: out.push_back(ch); break;
        }
    }
    return out;
}

std::string json_escape(std::string_view text) {
    std::string out;
    out.reserve(text.size() + 8);
    for (unsigned char ch : text) {
        switch (ch) {
            case '\\': out.append("\\\\"); break;
            case '"': out.append("\\\""); break;
            case '\b': out.append("\\b"); break;
            case '\f': out.append("\\f"); break;
            case '\n': out.append("\\n"); break;
            case '\r': out.append("\\r"); break;
            case '\t': out.append("\\t"); break;
            default:
                if (ch < 0x20) {
                    std::ostringstream oss;
                    oss << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(ch);
                    out.append(oss.str());
                } else {
                    out.push_back(static_cast<char>(ch));
                }
                break;
        }
    }
    return out;
}

uint32_t fnv1a_32(std::string_view text) {
    uint32_t hash = 2166136261u;
    for (unsigned char ch : text) {
        hash ^= ch;
        hash *= 16777619u;
    }
    return hash;
}

std::string safe_class_name(std::string_view text) {
    std::string safe;
    safe.reserve(text.size() + 12);
    safe.append("q_");
    for (unsigned char ch : text) {
        if (std::isalnum(ch)) {
            safe.push_back(static_cast<char>(ch));
        } else {
            safe.push_back('_');
        }
    }
    std::ostringstream oss;
    oss << std::hex << std::setw(8) << std::setfill('0') << fnv1a_32(text);
    safe.push_back('_');
    safe.append(oss.str());
    return safe;
}

double parse_double_or_zero(const std::string& value) {
    try {
        return std::stod(value);
    } catch (...) {
        return 0.0;
    }
}

std::vector<std::string> split_tsv_line(const std::string& line) {
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

std::string row_to_json_array(const ReportRow& row) {
    std::string json = "[";
    for (size_t i = 0; i < row.cols.size(); ++i) {
        if (i > 0) json.push_back(',');
        json.push_back('"');
        json.append(json_escape(row.cols[i]));
        json.push_back('"');
    }
    json.push_back(']');
    return json;
}

std::vector<QuerySummary> build_summaries(const std::vector<ReportRow>& rows) {
    std::map<std::string, QuerySummary> summary_map;
    std::map<std::string, std::set<std::string>> uniq_targets;

    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        const std::string& qname = row.cols[0];
        const std::string& search_type = row.cols[5];

        auto& summary = summary_map[qname];
        if (summary.query_name.empty()) {
            summary.query_name = qname;
            summary.safe_id = safe_class_name(qname);
        }
        ++summary.total;
        summary.row_indices.push_back(i);

        if (search_type == "Exact") {
            ++summary.exact;
        } else if (search_type == "Mismatch") {
            ++summary.mismatch;
        } else if (search_type == "Indel") {
            ++summary.indel;
        }

        if (search_type != "Exact") {
            const double score = parse_double_or_zero(row.cols[7]);
            summary.success_rate *= (1.0 - (score / 2.0));
        }
        uniq_targets[qname].insert(row.cols[2]);
    }

    std::vector<QuerySummary> ordered;
    ordered.reserve(summary_map.size());
    for (auto& entry : summary_map) {
        entry.second.uniq = uniq_targets[entry.first].size();
        ordered.push_back(std::move(entry.second));
    }

    std::sort(ordered.begin(), ordered.end(), [](const QuerySummary& lhs, const QuerySummary& rhs) {
        if (lhs.uniq != rhs.uniq) return lhs.uniq < rhs.uniq;
        return lhs.query_name < rhs.query_name;
    });
    return ordered;
}

void write_html_header(std::ostream& out) {
    const std::string_view logo_b64 = trim_ascii_whitespace(kReportLogoBase64);
    out << "<!DOCTYPE html>\n"
        << "<html lang=\"en\">\n"
        << "<head>\n"
        << "<meta charset=\"utf-8\">\n"
        << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        << "<title>SpacerScope Report</title>\n"
        << "<style>\n"
        << "*{box-sizing:border-box;}\n"
        << "body{margin:0;font-family:Segoe UI,Arial,sans-serif;background:#f4f7fb;color:#1f2937;}\n"
        << ".container{width:min(1400px,92vw);margin:24px auto;background:#fff;border-radius:16px;box-shadow:0 8px 24px rgba(15,23,42,.08);overflow:hidden;}\n"
        << "header{display:flex;align-items:center;gap:18px;padding:22px 28px;background:linear-gradient(135deg,#1d4ed8,#0f766e);color:#fff;}\n"
        << ".logo{height:96px;width:auto;flex:0 0 auto;}\n"
        << "h1{margin:0;font-size:28px;}\n"
        << ".sub{margin-top:6px;font-size:14px;opacity:.9;}\n"
        << ".toolbar{display:flex;align-items:center;gap:10px;flex-wrap:wrap;padding:16px 28px;background:#eef4ff;border-bottom:1px solid #dbe5f3;}\n"
        << ".btn{border:none;border-radius:8px;background:#1d4ed8;color:#fff;padding:8px 14px;cursor:pointer;font-size:14px;}\n"
        << ".btn:hover{background:#1e40af;}\n"
        << ".note{margin-left:auto;font-size:13px;color:#334155;}\n"
        << ".section{padding:24px 28px;}\n"
        << "table{width:100%;border-collapse:collapse;}\n"
        << "thead{background:#eff6ff;}\n"
        << "th,td{padding:10px 12px;border-bottom:1px solid #e5e7eb;text-align:center;vertical-align:top;}\n"
        << "td.left{text-align:left;}\n"
        << "tbody tr:nth-child(even){background:#fafcff;}\n"
        << "tbody tr:hover{background:#eef6ff;}\n"
        << "code{font-family:Consolas,Monaco,monospace;font-size:12px;}\n"
        << ".badge{display:inline-block;padding:3px 8px;border-radius:999px;font-size:12px;font-weight:600;color:#fff;}\n"
        << ".hot{background:#dc2626;}\n"
        << ".ok{background:#059669;}\n"
        << "details{border-top:1px solid #e5e7eb;padding:0 28px;}\n"
        << "summary{cursor:pointer;padding:18px 0;font-weight:600;}\n"
        << ".scroll-box{max-height:420px;overflow:auto;border:1px solid #dbe5f3;border-radius:12px;background:#fafcff;margin-bottom:24px;}\n"
        << ".detail-table{min-width:960px;}\n"
        << ".detail-table th{position:sticky;top:0;background:#fff;z-index:1;}\n"
        << "</style>\n"
        << "</head>\n"
        << "<body>\n"
        << "<div class=\"container\">\n"
        << "<header>\n"
        << "<img src=\"data:image/png;base64," << logo_b64 << "\" alt=\"SpacerScope\" class=\"logo\">\n"
        << "<div><h1>SpacerScope Match Report</h1>"
        << "<div class=\"sub\">Interactive summary, detailed rows, and export of selected matches.</div></div>\n"
        << "</header>\n";
}

void write_html_footer(std::ostream& out) {
    out << "<script>\n"
        << "function selectAllRows(flag){document.querySelectorAll('input[type=checkbox]').forEach(function(el){el.checked=flag;});}\n"
        << "function toggleMaster(master){document.querySelectorAll('.rowCk').forEach(function(el){el.checked=master.checked;});}\n"
        << "function toggleQuery(master){var safe=master.dataset.safe;document.querySelectorAll('.'+safe).forEach(function(el){el.checked=master.checked;});}\n"
        << "function toggleDetail(master,safe){document.querySelectorAll('.'+safe).forEach(function(el){el.checked=master.checked;});}\n"
        << "function exportSelected(){"
        << "const header=['QueryName','QuerySeq','TargetName','TargetSeq','Distance','SearchType','GenomicFrequency','Score','Details'];"
        << "const rows=[];"
        << "document.querySelectorAll('.detailCk:checked').forEach(function(el){"
        << "const raw=el.closest('tr').dataset.raw;"
        << "if(raw){rows.push(JSON.parse(raw));}"
        << "});"
        << "if(!rows.length){alert('Select rows before exporting.');return;}"
        << "const tsv=[header.join('\\t')].concat(rows.map(function(r){return r.join('\\t');})).join('\\n');"
        << "const blob=new Blob([tsv],{type:'text/tab-separated-values;charset=utf-8'});"
        << "const url=URL.createObjectURL(blob);"
        << "const a=document.createElement('a');a.href=url;a.download='selected_rows.tsv';a.click();URL.revokeObjectURL(url);"
        << "}\n"
        << "</script>\n"
        << "</body>\n"
        << "</html>\n";
}

}  // namespace

std::vector<ReportRow> build_report_rows(const std::vector<ScoredResult>& results) {
    std::vector<ReportRow> rows;
    rows.reserve(results.size());
    for (const auto& result : results) {
        ReportRow row;
        row.cols[0] = result.base.query_name;
        row.cols[1] = result.base.query_seq;
        row.cols[2] = result.base.target_name;
        row.cols[3] = result.base.target_seq;
        row.cols[4] = std::to_string(result.base.distance);
        row.cols[5] = result.base.search_type;
        row.cols[6] = std::to_string(result.base.genomic_frequency);
        row.cols[7] = format_score_cpp(result.score);
        row.cols[8] = alignment_events_to_json(result.details);
        rows.push_back(std::move(row));
    }
    return rows;
}

bool load_report_rows_from_tsv(const std::string& tsv_path,
                               std::vector<ReportRow>& rows,
                               std::string& error_message) {
    rows.clear();
    std::ifstream in(tsv_path, std::ios::binary);
    if (!in) {
        error_message = "Could not open scored TSV: " + tsv_path;
        return false;
    }

    std::string line;
    if (!std::getline(in, line)) {
        error_message = "Scored TSV is empty: " + tsv_path;
        return false;
    }

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> cols = split_tsv_line(line);
        if (cols.size() < 9) continue;
        ReportRow row;
        for (size_t i = 0; i < 9; ++i) row.cols[i] = std::move(cols[i]);
        rows.push_back(std::move(row));
    }
    return true;
}

bool write_html_report(const std::vector<ReportRow>& rows,
                       const std::string& html_path,
                       std::string& error_message) {
    std::ofstream out(html_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        error_message = "Could not open HTML output: " + html_path;
        return false;
    }

    const std::vector<QuerySummary> summaries = build_summaries(rows);

    write_html_header(out);

    out << "<div class=\"toolbar\">\n"
        << "<button class=\"btn\" onclick=\"selectAllRows(true)\">Select all</button>\n"
        << "<button class=\"btn\" onclick=\"selectAllRows(false)\">Clear all</button>\n"
        << "<button class=\"btn\" onclick=\"exportSelected()\">Export selected</button>\n"
        << "<div class=\"note\">Queries with uniq <= " << kHotThreshold << " are highlighted as hot.</div>\n"
        << "</div>\n";

    out << "<div class=\"section\">\n"
        << "<table>\n"
        << "<thead><tr>"
        << "<th><input type=\"checkbox\" id=\"masterCk\" onchange=\"toggleMaster(this)\"></th>"
        << "<th>QueryName</th><th>Total</th><th>Exact</th><th>Mismatch</th><th>Indel</th>"
        << "<th>uniq</th><th>Unique Editing Rate</th><th>Exact Sites</th>"
        << "</tr></thead>\n"
        << "<tbody>\n";

    for (const auto& summary : summaries) {
        const bool hot = static_cast<int>(summary.uniq) <= kHotThreshold;
        out << "<tr>"
            << "<td><input type=\"checkbox\" class=\"rowCk\" data-safe=\"" << html_escape(summary.safe_id)
            << "\" onchange=\"toggleQuery(this)\"></td>"
            << "<td class=\"left\"><a href=\"#q_" << html_escape(summary.safe_id)
            << "\" style=\"text-decoration:none;color:#1f2937;font-weight:600;\">"
            << html_escape(summary.query_name) << "</a></td>"
            << "<td>" << summary.total << "</td>"
            << "<td>" << summary.exact << "</td>"
            << "<td>" << summary.mismatch << "</td>"
            << "<td>" << summary.indel << "</td>"
            << "<td>" << summary.uniq << " <span class=\"badge " << (hot ? "hot" : "ok")
            << "\">" << (hot ? "hot" : "ok") << "</span></td>"
            << "<td>" << std::fixed << std::setprecision(2) << (summary.success_rate * 100.0) << "%</td>"
            << "<td>" << summary.exact << "</td>"
            << "</tr>\n";
    }

    out << "</tbody>\n</table>\n</div>\n";

    out << "<section class=\"section\" style=\"padding-top:0;\">\n"
        << "<h2 style=\"margin:0 0 16px;\">Details</h2>\n";

    for (const auto& summary : summaries) {
        out << "<details id=\"q_" << html_escape(summary.safe_id) << "\">\n"
            << "<summary>" << html_escape(summary.query_name) << " | total=" << summary.total
            << " | uniq=" << summary.uniq << "</summary>\n"
            << "<div class=\"scroll-box\">\n"
            << "<table class=\"detail-table\">\n"
            << "<thead><tr>"
            << "<th><input type=\"checkbox\" onchange=\"toggleDetail(this,'" << html_escape(summary.safe_id) << "')\"></th>"
            << "<th>QueryName</th><th>QuerySeq</th><th>TargetName</th><th>TargetSeq</th>"
            << "<th>Distance</th><th>SearchType</th><th>GenomicFrequency</th><th>Score</th><th>Details</th>"
            << "</tr></thead>\n<tbody>\n";

        for (size_t idx : summary.row_indices) {
            const ReportRow& row = rows[idx];
            const std::string row_json = row_to_json_array(row);
            out << "<tr data-raw=\"" << html_escape(row_json) << "\">"
                << "<td><input type=\"checkbox\" class=\"detailCk " << html_escape(summary.safe_id) << "\"></td>"
                << "<td>" << html_escape(row.cols[0]) << "</td>"
                << "<td><code>" << html_escape(row.cols[1]) << "</code></td>"
                << "<td>" << html_escape(row.cols[2]) << "</td>"
                << "<td><code>" << html_escape(row.cols[3]) << "</code></td>"
                << "<td>" << html_escape(row.cols[4]) << "</td>"
                << "<td>" << html_escape(row.cols[5]) << "</td>"
                << "<td>" << html_escape(row.cols[6]) << "</td>"
                << "<td>" << html_escape(row.cols[7]) << "</td>"
                << "<td><code>" << html_escape(row.cols[8]) << "</code></td>"
                << "</tr>\n";
        }

        out << "</tbody>\n</table>\n</div>\n</details>\n";
    }

    out << "</section>\n</div>\n";

    write_html_footer(out);
    out.flush();
    return true;
}

bool write_html_report_from_tsv(const std::string& tsv_path,
                                const std::string& html_path,
                                std::string& error_message) {
    std::vector<ReportRow> rows;
    if (!load_report_rows_from_tsv(tsv_path, rows, error_message)) return false;
    return write_html_report(rows, html_path, error_message);
}
