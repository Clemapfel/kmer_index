//
// Created by clem on 6/1/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <random>

#include <benchmark/benchmark.h>

size_t seed = 1234;

// compare different pow implementation

namespace local
{
    constexpr size_t recursive_pow(size_t base, uint8_t exp)
    {
        if (exp == 0)
            return 1ul;

        size_t temp = recursive_pow(base, exp / 2ul);

        if ((exp & 1ul) == 0) // i % n = i & n-1
            return temp * temp;
        else
            return base * temp * temp;
    }

    constexpr size_t bit_pow(size_t base, uint8_t exp)
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

    constexpr int64_t switch_pow(int64_t base, uint8_t exp)
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
                if (exp & 1ul) result *= base;
                exp >>= 1ul;
                base *= base;
            case 5:
                if (exp & 1ul) result *= base;
                exp >>= 1ul;
                base *= base;
            case 4:
                if (exp & 1ul) result *= base;
                exp >>= 1ul;
                base *= base;
            case 3:
                if (exp & 1ul) result *= base;
                exp >>= 1ul;
                base *= base;
            case 2:
                if (exp & 1ul) result *= base;
                exp >>= 1;
                base *= base;
            case 1:
                if (exp & 1ul) result *= base;
            default:
                return result;
        }
    }
}

// std::pow
static void std_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist{};
    std::uniform_int_distribution<uint8_t> exp_dist{};
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::pow(base_dist(engine), exp_dist(engine)));
    }
}

// fast_pow
static void switch_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist{};
    std::uniform_int_distribution<uint8_t> exp_dist{};
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(local::switch_pow(base_dist(engine), exp_dist(engine)));
    }
}

// bit_pow
static void bit_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist{};
    std::uniform_int_distribution<uint8_t> exp_dist{};
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(local::bit_pow(base_dist(engine), exp_dist(engine)));
    }
}

// bit_pow
static void recursive_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist{};
    std::uniform_int_distribution<uint8_t> exp_dist{};
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(local::recursive_pow(base_dist(engine), exp_dist(engine)));
    }
}

bool check_correctness_first = true;

int main(int argc, char** argv)
{
    size_t seed_bp = seed;

    // check correctness
    for (size_t i = 0; i < 1e7; ++i)
    {
        std::uniform_int_distribution<size_t> base_dist{};
        std::uniform_int_distribution<uint8_t> exp_dist{};
        std::mt19937 engine(seed_bp++);

        auto base = base_dist(engine);
        auto exp = exp_dist(engine);

        auto correct = std::pow(base, exp);

        if (local::bit_pow(base, exp) != correct)
        {
            std::cout << "bit failed\n";
            exit(1);
        }
    }


    benchmark::RegisterBenchmark("std_pow", std_pow);
    benchmark::RegisterBenchmark("switch_pow", fast_pow);
    benchmark::RegisterBenchmark("recursive_pow", fast_pow);
    benchmark::RegisterBenchmark("bit_pow", fast_pow);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/pow_vs_pow/raw.csv");
}


