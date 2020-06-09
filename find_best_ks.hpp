// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <c++/7/cstddef>
#include <array>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <map>

namespace kmer
{
    // function that picks k based on query input
    constexpr std::array<size_t, 5> possible_primes = {5, 7, 11, 13, 17};

    template<size_t n, std::ranges::range range_t>
    std::array<size_t, n> pick_optimal_ks(const range_t query_sizes)
    {
        assert(n < possible_primes.size());

        std::array<size_t, n> hist;

        for (size_t size : query_size)
        {
            size_t min_i = 0;
            size_t min = query_size % possible_primes[0];

            size_t current_i = 0;
            for (size_t i = 1; i < possible_primes.size(); ++i)
            {
                auto mod = query_size % possible_primes.at(i);
                if (mod < min)
                {
                    min = mod;
                    min_i = i;
                }
                ++i;
            }

            hist[min_i]++;
        }


    }

}