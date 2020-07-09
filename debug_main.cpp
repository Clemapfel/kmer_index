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
#include <thread_pool.hpp>
#include <choose_best_k.hpp>
#include <thread_pool.hpp>

std::atomic<size_t> seed = 200;
size_t text_length = 1000000;
int main()
{
    using namespace seqan3;

    constexpr size_t n = 107;

    std::vector<size_t> _all_ks = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
    std::sort(_all_ks.begin(), _all_ks.end(), [](size_t a, size_t b) -> bool {return a > b;});

    std::vector<size_t> not_found{};

    // first pass
    std::array<std::pair<size_t, uint64_t>, n> _optimal_k;
    for (size_t query_size = 0; query_size < _optimal_k.size(); ++query_size)
    {
        // pick best k to search withthis is rea
        size_t optimal_k = _all_ks.front();
        bool found = false;

        for (size_t k : _all_ks)
        {
            // for actual kmers, prioritze absolute distance rather than divisibility
            if (query_size < 31)
            {
                if (query_size <= k and (k - query_size < optimal_k - query_size))
                {
                    optimal_k = k;
                    found = true;
                    continue;
                }
            }
            // for long queries divisibility is more important, also prefer higher k to lower k if mod is equal
            else
            {
                if (query_size % k == 0)
                {
                    optimal_k = k;
                    found = true;
                    break;
                }
            }
        }

        if (found) {
            _optimal_k[query_size] = std::make_pair(optimal_k, 0);
        }
        else {
            not_found.push_back(query_size);
        }
    }

    // second pass
    std::vector<size_t> not_found_2{};
    for (size_t query_size : not_found)
    {
        size_t optimal_k = _all_ks.front();
        bool found = false;

        for (size_t k : _all_ks)
        {
            if (query_size % k >= 10 and
                std::binary_search(_all_ks.begin(), _all_ks.end(), query_size % k, [](auto a, auto b) -> bool { return a > b; }))
            {
                optimal_k = k;
                found = true;
                break;
            }
        }

        if (found) {
            _optimal_k[query_size] = std::make_pair(optimal_k, query_size % optimal_k);
        }
        else {
            not_found_2.push_back(query_size);
        }
    }

    // third pass
    for (size_t query_size : not_found_2)
    {
        size_t optimal_k = _all_ks.front();
        bool found = false;

        for (size_t k : _all_ks)
        {
            for (size_t k : _all_ks)
            {
                if ((ceil(query_size / float(k)) * k - query_size) <
                    (ceil(query_size / float(optimal_k)) * optimal_k - query_size))
                    optimal_k = k;
            }
        }

        _optimal_k[query_size] = std::make_pair(optimal_k, -1 * (query_size % optimal_k));
    }

    seqan3::debug_stream << "m"

    for (size_t i = 1; i < n; ++i)
    {
        size_t k = _optimal_k[i].first;
        int rest = _optimal_k[i].second;
        seqan3::debug_stream << i << " | " << _optimal_k[i] << " \t\t(" << i/k << "*" << k << "= " << (i/k)*k << ")\n";
    }
}