// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <cstdint>

// copied from: https://gist.github.com/orlp/3551590 (slightly altered to be constexpr)

namespace kmer::detail
{
    // fastest pow according to pow benchmark c.f. pow_vs_pow/main.cpp
    constexpr size_t fast_pow(size_t base, uint8_t exp)
    {
        unsigned long result = 1;
        while (exp > 0)
        {
            if ((exp & 1ul) == 1)   // i % n = i & n-1
                result = result * base;

            exp = exp >> 1ul;
            base = base * base;
        }

        return result;
    }
}
