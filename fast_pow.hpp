#pragma once

#include <cstdint>
#include <cstdlib>

// source: https://gist.github.com/orlp/3551590

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
            6, 6, 6, 6, 6, 6, 6, 255,
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

    constexpr size_t fast_pow(size_t base, uint8_t exp)
    {
        size_t result = 1;

        // instead of loop result will fall through the switch starting at position according to highest_bit_set
        // will skip all computation if result would overflow
        switch (highest_bit_set[exp])
        {
            case 255: // overflow
                if (base == 1)
                {
                    return 1;
                }
                else
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
