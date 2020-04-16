//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <input_generator.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/algorithm/search.hpp>

#include <thread>
#include <chrono>
#include <memory>

using namespace seqan3;

struct safe_print
{
    static void print(std::string str)
    {
        _mutex.lock();
        std::cout << str << "\n";
        _mutex.unlock();
    }

    private:
        static inline std::mutex _mutex{};
};



void thread_practice()
{
    auto text = input_generator<dna4>::generate_sequence(100000);


    // create vector for operations
    const size_t vec_size = 10000;
    std::vector<int> vec;
    vec.reserve(vec_size);

    for (size_t i = 0; i < vec_size; ++i)
    {
        vec.emplace_back(rand());
    }

    std::vector<std::thread> _threads;
    auto n_threads = std::thread::hardware_concurrency();
    std::cout << n_threads << " threads supported\n\n";

    for (size_t i = 0; i < n_threads; ++i)
    {
        _threads.emplace_back(&thread_print);
        _names[_threads.back().get_id()] = "thread " + std::to_string(i);
    }

    std::this_thread::sleep_for(std::chrono::seconds(10));
}

int main()
{
    thread_practice();
}
