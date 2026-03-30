#pragma once

#include "postprocess.hpp"

#include <array>
#include <string>
#include <vector>

struct ReportRow {
    std::array<std::string, 9> cols;
};

std::vector<ReportRow> build_report_rows(const std::vector<ScoredResult>& results);

bool load_report_rows_from_tsv(const std::string& tsv_path,
                               std::vector<ReportRow>& rows,
                               std::string& error_message);

bool write_html_report(const std::vector<ReportRow>& rows,
                       const std::string& html_path,
                       std::string& error_message);

bool write_html_report_from_tsv(const std::string& tsv_path,
                                const std::string& html_path,
                                std::string& error_message);

