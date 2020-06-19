// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <vector>
#include <array>
#include <set>
#include <iostream>
#include <algorithm>

// pick best multiple k (TODO: make constexpr)
template<seqan3::alphabet alphabet_t, std::ranges::range range_t>
std::vector<size_t> choose_best_k(
        range_t interval,             // [in] vector of query sizes
        size_t n_k = 4                // [in] number of k to choose
)
{
    std::vector<std::pair<size_t, size_t>> k_and_score;

    //if (alphabet_t <)

    // inversely ordered so loop prioritizes high k
    for (size_t k : {29, 27, 25, 23, 21, 19, 17, 13, 11, 10})
        k_and_score.emplace_back(std::make_pair(k, 0));

    for (size_t i : interval)
    {
        for (auto& p : k_and_score)
        {
            size_t k = p.first;

            // best case: 4 pts
            if (i % k == 0)
            {
                p.second += 3;
                break;
            }
            // good case: p = n*k - 1 or -2
            else if (k - (i % k) <= 3)
            {
                p.second += 4 - (k - (i % k));
                break;
            }
            // try to find a better one
            else
                continue;
        }
    }

    std::sort(k_and_score.begin(), k_and_score.end(),
            [](auto a, auto b) -> bool {return a.second > b.second;});

    seqan3::debug_stream << "(k, score): " << k_and_score << "\n";

    std::vector<size_t> output;
    for (size_t i = 0; i < n_k; ++i)
        output.push_back(k_and_score.at(i).first);

    return output;
}