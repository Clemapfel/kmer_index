// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <vector>
#include <array>
#include <set>
#include <iostream>
#include <algorithm>

// pick best multiple k (TODO: make constexpr)
std::vector<size_t> choose_best_k(
        std::vector<size_t> interval, // [in] vector of query sizes
        size_t n_k = 4                // [in] number of k to choose
)
{
    std::vector<std::pair<size_t, size_t>> k_and_score;

    // primes inversely ordered so loop prioritizes high k
    for (size_t k : {25, 23, 20, 19, 17, 14, 13, 11})
        k_and_score.emplace_back(std::make_pair(k, 0));

    size_t not_covered = 0;

    for (size_t i : interval)
    {
        bool found = false;
        for (auto& p : k_and_score)
        {
            size_t k = p.first;

            // best case: 4 pts
            if (i % k == 0)
            {
                p.second += 3;
                found = true;
                break;
            }
            // good case: p = n*k - 1 or -2
            else if (k - (i % k) < 3)
            {
                p.second += 3 - (k - (i % k));
                found = true;
                break;
            }
            // try to find a better one
            else
                continue;
        }

        if (not found)
            not_covered++;
        }
    }

    std::sort(k_and_score.begin(), k_and_score.end(),
            [](auto a, auto b) -> bool {return a.second > b.second;});

    std::cout << "(k, score): " << std::to_string(stdk_and_score) << "\n";
    std::cout << "not covered: " << std::to_string(not_covered) << "/" << std::to_string(interval.size()) << "\n";

    std::vector<size_t> output;
    for (size_t i = 0; i < n_k; ++i)
        output.push_back(k_and_score.at(i).first);

    return output;
}