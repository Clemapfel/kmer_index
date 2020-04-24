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
#include <future>
#include <functional>
#include <thread_pool.hpp>

using namespace seqan3;

std::map<std::thread::id, size_t> _thread_ids;

template<typename T>
void sync_print(T t)
{
    auto id = std::this_thread::get_id();

    if (_thread_ids.find(id) == _thread_ids.end())
        _thread_ids.insert(std::make_pair(id, _thread_ids.size()));

    debug_stream << t << "printed by (" << _thread_ids.at(id) << ")\n";
}

int main()
{
    thread_pool pool{8};

    for (size_t i = 0; i < 8; ++i)
        auto fut = pool.execute(sync_print<std::string>, "test");

    //pool.wait_to_finish();

}
