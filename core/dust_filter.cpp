#include "dust_filter.hpp"

#include <algorithm>
#include <cctype>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

namespace {

std::mutex g_dust_mutex;

int base_to_idx(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

} // namespace

bool is_complex(std::string_view seq, int threshold) {
    if (seq.length() < 3) return true;

    uint16_t counts[64] = {0};
    for (size_t j = 0; j + 2 < seq.length(); ++j) {
        const int b1 = base_to_idx(seq[j]);
        const int b2 = base_to_idx(seq[j + 1]);
        const int b3 = base_to_idx(seq[j + 2]);
        if (b1 != -1 && b2 != -1 && b3 != -1) {
            counts[(b1 << 4) | (b2 << 2) | b3]++;
        }
    }

    int score = 0;
    for (int c : counts) {
        if (c > 1) {
            score += c * (c - 1) / 2;
        }
    }
    return score <= threshold;
}

void dustFilter(const SequenceList& input, SequenceList& output, int threshold, int num_threads) {
    const size_t total_size = input.size();
    if (total_size == 0) return;

    if (num_threads <= 1) {
        for (size_t i = 0; i < total_size; ++i) {
            if (is_complex(input.sequences[i], threshold)) {
                output.push_back(input.headers[i], input.sequences[i]);
            }
        }
        return;
    }

    const size_t chunk_size = (total_size + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);
    std::vector<std::future<void>> futures;

    for (int t = 0; t < num_threads; ++t) {
        const size_t start = static_cast<size_t>(t) * chunk_size;
        const size_t end = (std::min)(start + chunk_size, total_size);
        if (start >= end) break;

        futures.push_back(std::async(std::launch::async, [&input, &output, threshold, start, end]() {
            SequenceList local_list;
            for (size_t i = start; i < end; ++i) {
                if (is_complex(input.sequences[i], threshold)) {
                    local_list.push_back(input.headers[i], input.sequences[i]);
                }
            }

            if (local_list.size() > 0) {
                std::lock_guard<std::mutex> lock(g_dust_mutex);
                for (size_t i = 0; i < local_list.size(); ++i) {
                    output.push_back(local_list.headers[i], local_list.sequences[i]);
                }
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
}
