//
// Created by clem on 6/1/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <cmath>

#include <random>
#include <iostream>

#include <benchmark/benchmark.h>

size_t seed = 1234;

// compare different pow implementation

namespace local
{
    constexpr size_t trivial_pow(size_t base, size_t exp)
    {
        size_t result = 1;
        while (exp > 0)
        {
            result *= base;
            exp--;
        }

        return result;
    }

    constexpr size_t recursive_pow(size_t base, size_t exp)
    {
        if (exp == 0)
            return 1ul;

        size_t temp = recursive_pow(base, exp / 2ul);

        if ((exp & 1ul) == 0) // i % n = i & n-1
            return temp * temp;
        else
            return base * temp * temp;
    }

    constexpr size_t bit_pow(size_t base, size_t exp)
    {
        size_t result = 1;
        while (exp > 0)
        {
            if (exp & 1ul)   // i % n = i & n-1
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

    constexpr size_t switch_pow(size_t base, uint8_t exp)
    {
        size_t result = 1;

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

size_t base_min = 1, base_max = 1e6;
uint32_t exp_min = 0, exp_max = 128;

// std::pow
static void trivial_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist(base_min, base_max);
    std::uniform_int_distribution<uint8_t> exp_dist(exp_min, exp_max);
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(local::trivial_pow(base_dist(engine), exp_dist(engine)));
    }
}

// fast_pow
static void switch_pow(benchmark::State& state)
{
    state.counters["seed"] = seed;

    std::uniform_int_distribution<size_t> base_dist(base_min, base_max);
    std::uniform_int_distribution<uint8_t> exp_dist(exp_min, exp_max);
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

    std::uniform_int_distribution<size_t> base_dist(base_min, base_max);
    std::uniform_int_distribution<uint8_t> exp_dist(exp_min, exp_max);
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

    std::uniform_int_distribution<size_t> base_dist(base_min, base_max);
    std::uniform_int_distribution<uint8_t> exp_dist(exp_min, exp_max);
    std::mt19937 engine(seed++);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(local::recursive_pow(base_dist(engine), exp_dist(engine)));
    }
}

int main(int argc, char** argv)
{

    benchmark::RegisterBenchmark("x*x", trivial_pow);
    benchmark::RegisterBenchmark("switch_pow", switch_pow);
    benchmark::RegisterBenchmark("recursive_pow", recursive_pow);
    benchmark::RegisterBenchmark("bit_pow", bit_pow);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/pow_vs_pow/raw.csv");
}

/*
    if (check_correctness_first)
    {
        size_t seed_bp = seed;

        // check correctness
        for (size_t i = 0; i < 1e9; ++i)
        {
            std::cout << "starting correctness test\n";

            std::uniform_int_distribution<size_t> base_dist{1, 10000};
            std::uniform_int_distribution<uint8_t> exp_dist{0, 32};
            std::mt19937 engine(seed_bp++);

            size_t base = base_dist(engine);
            uint8_t exp = exp_dist(engine);

            auto trivial_res = local::trivial_pow(base, exp);
            auto bit_res = local::bit_pow(base, exp);
            auto switch_res = local::switch_pow(base, exp);
            auto recursive_res = local::recursive_pow(base, exp);
            auto std_res = powl(base, exp);

            std::cout
            << "x*x       : " << trivial_res << "\n"
            << "std       : " << std_res << "\n"
            << "bit       : " << bit_res << "\n"
            << "switch    : " << switch_res << "\n"
            << "recursive : " << recursive_res << "\n";
        }

        std::cout << "correctness test succeeded\n";
    }
 */


