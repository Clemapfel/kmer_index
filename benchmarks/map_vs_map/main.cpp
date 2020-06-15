#include <vector>
#include <benchmark/benchmark.h>
#include <random>
#include <benchmarks/cleanup_csv.hpp>
#include <fast_pow.hpp>
#include <iostream>

#include <robin_hood.h>
#include <unordered_map>
#include <absl/container/node_hash_map.h>
#include <boost/unordered_map.hpp>

namespace kmer
{
    template<typename map_t, typename key_t=size_t, typename value_t=std::vector<uint32_t>>
    class map_proxy
    {
        using inner_key_t = uint32_t;

        private:
            map_t _data;

            const value_t _empty{};

        public:
            map_proxy() = default;

        void insert(key_t hash, inner_key_t pos)
        {
            if (_data.find(hash) == _data.end())
                _data.emplace(std::make_pair(hash, std::vector<uint32_t>()));

            _data[hash].push_back(pos);
        }

        const value_t* at(key_t hash)
        {
            if (_data.find(hash) == _data.end())
                return &_empty;
            else
                return &_data.at(hash);
        }
    };

    using hash_t = size_t;
    using vec_t = std::vector<uint32_t>;

    using std_proxy = map_proxy<std::unordered_map<hash_t, vec_t>>;
    using boost_proxy = map_proxy<boost::unordered_map<hash_t, vec_t>>;
    using robin_hood_proxy = map_proxy<robin_hood::unordered_map<hash_t, vec_t>>;
    using abseil_proxy = map_proxy<absl::node_hash_map<hash_t, vec_t>>;

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

    class boost_proxy : map_proxy
    {
        private:
            boost::unordered_map<size_t, std::vector<uint32_t>> _data;

        public:
            boost_proxy() = default;

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

    class abseil_proxy : map_proxy
    {
        private:
            absl::node_hash_map<size_t, std::vector<uint32_t>> _data;

        public:
            abseil_proxy() = default;

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

size_t seed = 1234;

static void std_at(benchmark::State& state, size_t size)
{
    state.counters["seed"] = seed;
    state.counters["size"] = size;
    state.counters["max_hash"] = MAX_HASH;

    kmer::std_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < size; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

static void abseil_at(benchmark::State& state, size_t size)
{
    state.counters["seed"] = seed;
    state.counters["size"] = size;
    state.counters["max_hash"] = MAX_HASH;

    kmer::abseil_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < size; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

static void boost_at(benchmark::State& state, size_t size)
{
    state.counters["seed"] = seed;
    state.counters["size"] = size;
    state.counters["max_hash"] = MAX_HASH;

    kmer::boost_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < size; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

static void robin_hood_at(benchmark::State& state, size_t size)
{
    state.counters["seed"] = seed;
    state.counters["size"] = size;
    state.counters["max_hash"] = MAX_HASH;

    kmer::robin_hood_proxy map;
    std::uniform_int_distribution<size_t> hash_dist(MIN_HASH, MAX_HASH);
    std::uniform_int_distribution<uint32_t> pos_dist(0, std::numeric_limits<uint32_t>::max() - 2);

    std::mt19937 engine(seed++);

    // fill first
    for (size_t i = 0; i < size; ++i)
        map.insert(hash_dist(engine), pos_dist(engine));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(map.at(hash_dist(engine)));
    }
}

//./MAP_BENCHMARK_K10 --benchmark_format=console --benchmark_out=/srv/public/clemenscords/map_vs_map/raw.csv --benchmark_out_format=csv --benchmark_repetitions=10 --benchmark_report_aggregates_only=false

// main
int main(int argc, char** argv)
{
    // at
    for (size_t n = 0; n <= 2000000; n += 50000)
    {
        benchmark::RegisterBenchmark("da_at", da_at, n);
        benchmark::RegisterBenchmark("robin_hood_at", robin_hood_at, n);
        benchmark::RegisterBenchmark("std_at", std_at, n);

        n_benchmarks += 3;
    }

    std::cout << std::to_string(n_benchmarks) << " benchmarks registered.\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    std::cout << "done.\n";

    cleanup_csv("/srv/public/clemenscords/map_vs_map/raw.csv");
}
