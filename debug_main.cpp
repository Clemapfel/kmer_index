//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <functional>
#include <set>
#include <thread_pool.hpp>
#include <choose_best_k.hpp>
#include <thread_pool.hpp>

std::atomic<size_t> seed = 200;
size_t text_length = 1000000;

int main()
{
    using namespace seqan3;

    constexpr size_t n = 3000;

    std::vector<size_t> _all_ks = {3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19};
    std::sort(_all_ks.begin(), _all_ks.end(), [](size_t a, size_t b) -> bool { return a > b; });

    std::vector<size_t> high_ks;
    for (size_t i : _all_ks)
        if (i >= 9)
            high_ks.push_back(i);

    std::array<std::vector<size_t>, n> _optimal_k;
    std::fill(_optimal_k.begin(), _optimal_k.end(), std::vector<size_t>());

    std::array<bool, n> _found_in_first_pass;
    std::fill(_found_in_first_pass.begin(), _found_in_first_pass.end(), false);

    std::vector<size_t> not_found{};

    // first pass: dynamic programming
    for (size_t k : high_ks)
    {
        _optimal_k[k] = {k};
        _found_in_first_pass[k] = true;
    }

    for (size_t q = _all_ks.front()+1; q < n; ++q)
    {
        bool found = false;
        for (size_t k : high_ks)
        {
            if (not _optimal_k[q - k].empty())
            {
                _optimal_k[q] = _optimal_k[q - k];
                _optimal_k[q].push_back(k);
                found = true;
                break;
            }
        }

        if (found)
            _found_in_first_pass[q] = true;
    }

    // second pass: fill rest that can only be searched with subk
    for (size_t q = 1; q < n; ++q)
    {
        if (not _optimal_k[q].empty())
            continue;

        if (q < _all_ks.front())
        {
            size_t optimal_k = _all_ks.front();
            for (size_t k : _all_ks)
            {
                if (q <= k and (k - q < optimal_k - q))
                {
                    optimal_k = k;
                    continue;
                }
            }
            _optimal_k[q] = {optimal_k};
        }
        else
        {
            size_t optimal_k = _all_ks.front();
            for (size_t k : _all_ks)
            {
                if ((ceil(q / float(k)) * k - q) <
                    (ceil(q / float(optimal_k)) * optimal_k - q))
                    optimal_k = k;
            }

            _optimal_k[q] = {optimal_k};
        }
    }

    // debug print
    for (size_t i = 0; i < 300; ++i)
        seqan3::debug_stream << i << " | " << (_found_in_first_pass[i] ? _optimal_k[i] : std::vector<size_t>{}) << "\n";
}
