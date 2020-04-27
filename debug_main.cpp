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

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <functional>
#include <thread_pool.hpp>
#include <thread_pool_test.hpp>

#include <syncstream.hpp>

using namespace seqan3;

int main()
{
    thread_pool pool{8};

    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 6000; ++i)
        futures.emplace_back(pool.execute(sync_print<std::string>, "success #" + std::to_string(i)));

    std::this_thread::sleep_for(std::chrono::microseconds(1));

    //pool.resize(16);

    for (auto& f : futures)
        f.get();
    //pool.wait_to_finish();
}
