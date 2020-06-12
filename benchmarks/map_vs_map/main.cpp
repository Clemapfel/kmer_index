#include <vector>
#include <benchmark/benchmark.h>
#include <random>
#include <benchmarks/cleanup_csv.hpp>
#include <unordered_map>
#include <robin_hood.h>
#include <fast_pow.hpp>

namespace kmer
{
    struct map_proxy
    {
        protected:
            const std::vector<uint32_t> _empty{};

        public:
            virtual void insert(size_t hash, uint32_t pos) = 0;
            virtual const std::vector<uint32_t>* at(size_t hash) = 0;
    };

    class direct_addressing_proxy : map_proxy
    {
        private:
            std::vector<std::vector<uint32_t>> _data;
            size_t _min_hash = std::numeric_limits<size_t>::max(), _max_hash = 0;

        public:
            direct_addressing_proxy()
            {
                _data.emplace_back();
            }

            void insert(size_t hash, uint32_t pos) override
            {
                // resize if necessary
                if (hash < _min_hash)
                {
                    _data.insert(_data.begin(), std::labs(int(_min_hash) - int(hash)), std::vector<uint32_t>{});
                    _min_hash = hash;
                }
                else if (hash > _max_hash)
                {
                    _data.insert(_data.end(), std::labs(int(hash) - int(_max_hash)), std::vector<uint32_t>{});
                    _max_hash = hash;
                }

                _data[hash - _min_hash].push_back(pos);
            }

            const std::vector<uint32_t>* at(size_t hash) override
            {
                if (hash < _min_hash or hash > _max_hash)
                    return &_empty;
                else
                    return &_data.at(hash - _min_hash);
            }
    };

    class std_proxy : map_proxy
    {
        private:
            std::unordered_map<size_t, std::vector<uint32_t>> _data;

        public:
            std_proxy() = default;

            void insert(size_t hash, uint32_t pos) override
            {
                if (_data.find(hash) == _data.end())
                    _data.emplace(std::make_pair(hash, std::vector<uint32_t>()));

                _data[hash].push_back(pos);
            }

            const std::vector<uint32_t>* at(size_t hash) override
            {
                if (_data.find(hash) == _data.end())
                    return &_empty;
                else
                    return &_data.at(hash);
            }
    };

    class robin_hood_proxy : map_proxy
    {
        private:
            robin_hood::unordered_map<size_t, std::vector<uint32_t>> _data;

        public:
            robin_hood_proxy() = default;

            void insert(size_t hash, uint32_t pos) override
            {
                if (_data.find(hash) == _data.end())
                    _data.emplace(hash, std::vector<uint32_t>());

                _data[hash].push_back(pos);
            }

            const std::vector<uint32_t>* at(size_t hash) override
            {
                if (_data.find(hash) == _data.end())
                    return &_empty;
                else
                    return &_data.at(hash);
            }
    };
}

// #################################################################

using hash_t = size_t;
using value_t = std::vector<uint32_t>;

constexpr size_t MIN_HASH = 0;
constexpr size_t MAX_HASH = kmer::detail::fast_pow(4, 7);    // dna4 for k=exp

constexpr size_t N_LEAFS = 1000000;

size_t seed = 1234;

// direct addressing
static void std_insertion(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::std_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        map.insert(hash_dist(engine), pos_dist(engine));
    }
}

static void std_at(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::std_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

// direct addressing
static void da_insertion(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::direct_addressing_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        map.insert(hash_dist(engine), pos_dist(engine));
    }
}

static void da_at(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::std_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

// robin hood
static void robin_hood_insertion(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::robin_hood_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        map.insert(hash_dist(engine), pos_dist(engine));
    }
}

static void robin_hood_at(benchmark::State& state)
{
    state.counters["seed"] = seed;

    kmer::std_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < N_LEAFS; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

//./MAP_BENCHMARK_k10 --benchmark_format=console --benchmark_out=/srv/public/clemenscords/map_vs_map/raw.csv --benchmark_out_format=csv --benchmark_repetitions=100 --benchmark_report_aggregates_only=false

// main
int main(int argc, char** argv)
{
    benchmark::RegisterBenchmark("da_insert", da_insertion);
    benchmark::RegisterBenchmark("robin_hood_insert", robin_hood_insertion);
    benchmark::RegisterBenchmark("std_insert", std_insertion);

    benchmark::RegisterBenchmark("da_at", da_at);
    benchmark::RegisterBenchmark("robin_hood_at", robin_hood_at);
    benchmark::RegisterBenchmark("std_at", std_at);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/map_vs_map/raw.csv");
}
