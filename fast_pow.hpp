// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <cstdint>

// copied from: https://gist.github.com/orlp/3551590 (slightly altered to be constexpr)

namespace kmer::detail
{
    constexpr uint8_t highest_bit_set[] =
    {
            0, 1, 2, 2, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 255, // anything past 63 is a guaranteed overflow with base > 1
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255,
    };

    constexpr int64_t fast_pow(int64_t base, uint8_t exp)
    {
        int64_t result = 1;

        switch (highest_bit_set[exp])
        {
            case 255: // we use 255 as an overflow marker and return 0 on overflow/underflow
                if (base == 1)
                {
                    return 1;
                }

                if (base == -1)
                {
                    return 1 - 2 * (exp & 1);
                }

                return 0;
            case 6:
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            case 5:
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            case 4:
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            case 3:
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            case 2:
                if (exp & 1) result *= base;
                exp >>= 1;
                base *= base;
            case 1:
                if (exp & 1) result *= base;
            default:
                return result;
        }
    }
}
