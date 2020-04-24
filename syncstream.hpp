
#pragma once

#include <seqan3/core/debug_stream.hpp>
#include <thread_safe_queue.hpp>

#include <map>
#include <thread>

inline std::map<std::thread::id, size_t> _thread_ids;

template<typename T>
void sync_print(T t)
{
    auto id = std::this_thread::get_id();

    if (_thread_ids.find(id) == _thread_ids.end())
        _thread_ids.insert(std::make_pair(id, _thread_ids.size()));

    auto out = "[" + std::to_string(_thread_ids.at(id)) + "] " + std::string(t) + "\n";

    std::cout << out;
}