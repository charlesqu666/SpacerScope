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

