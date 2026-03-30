#include "fasta_parser.hpp"
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
