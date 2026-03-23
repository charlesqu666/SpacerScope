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
﻿#include "fasta_parser.hpp"
#include "types.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

SequenceList parse_fasta(const std::string& filename) {
    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    SequenceList sequences;
    kseq_t* seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        std::string sequence_str(seq->seq.s);
        std::transform(sequence_str.begin(),
                       sequence_str.end(),
                       sequence_str.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
        sequences.push_back(seq->name.s, sequence_str);
    }

    kseq_destroy(seq);
    gzclose(fp);
    return sequences;
}
