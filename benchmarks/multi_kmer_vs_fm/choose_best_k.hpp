// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <vector>
#include <array>
#include <set>

// pick best multiple k (TODO: make constexpr)
std::vector<size_t> choose_best_k(std::vector<size_t> interval, size_t n_k = 4)
{
    std::vector<std::pair<size_t, size_t>> k_and_score;
    for (size_t k : {25, 23, 20, 19, 17, 14, 13, 11})
        k_and_score.emplace_back(std::make_pair(k, 0));

    std::set<size_t> uncovered;
    size_t min = std::numeric_limits<size_t>::max();
    size_t max = 0;
    for (size_t i : interval)
    {
        if (i < min)
            min = i;
        if (i > max)
            max = i;

        bool found = false;
        for (auto p : k_and_score)
        {
            size_t k = p.first;
            if (i == k-1 or i == k-2)
            {
                for (auto& j : k_and_score)
                {
                    if (j.first == k)
                    {
                        j.second += (3 - (k - i));
                        break;
                    }
                }
                found = true;
                break;
            }
            else if (i % k <= 3)
            {
               for (auto& j : k_and_score)
               {
                   if (j.first == k)
                   {
                       j.second += (4 - i % k);
                       break;
                   }
               }
               found = true;
               break;               // prioritize higher ks because of n*k search
            }
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