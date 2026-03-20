#pragma once

#include "types.hpp"

#include <string_view>

// Return true when the sequence passes the low-complexity filter.
bool is_complex(std::string_view seq, int threshold);

// Apply the DUST-based complexity filter in parallel.
void dustFilter(const SequenceList& input,
                SequenceList& output,
                int threshold,
                int num_threads = 1);
