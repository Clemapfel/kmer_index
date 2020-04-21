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

#include <thread_pool.hpp>

using namespace seqan3;

// single consumer single producer queue
template<typename T>
struct thread_safe_queue
{
    public:
        void push_back(T t)
        {
            {
                std::lock_guard<std::mutex> lock(_insert_mutex);    // unlocks at end of block scope
                _data.push_back(t);
            }
        }

        std::optional<T> pop_front()
        {
            if (_currently_aborting || _data.empty())
                return std::move(std::optional<T>(std::nullopt));

            std::lock_guard<std::mutex> lock(_insert_mutex); // not currently inserting
            auto output = std::move(_data.front());
            _data.pop_front();
            return std::move(output);
        }

        ~thread_safe_queue()
        {
            abort();
        };

    protected:
        void abort()
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            _currently_aborting = true;
            _data.clear();
        }

    private:
        std::deque<T> _data;
        std::mutex _insert_mutex;
        std::condition_variable _cv;

        bool _currently_aborting = false;
};

template<typename T>
void sync_print(T s)
{
    debug_stream << s << "\n";
}

template<typename function_t, typename... args_t>
auto execute(function_t* f, args_t... args )
{
    debug_stream << "01\n";
    auto thr = std::thread([&f, &args...](){f(args...);});
    thr.join();

    // create task that holds function and args
    // wrap into lambda void (int)
    // store lambda in queue
}

template<typename function_t>
auto execute(function_t&& f)
{
    debug_stream << "02\n";
    auto thr = std::thread(f);
    thr.join();
}




int main()
{

    //execute(sync_print<decltype(std::this_thread::get_id())>, std::this_thread::get_id());
    //execute([](){sync_print(std::this_thread::get_id());});

    thread_safe_queue<int> q;

    auto thr1 = std::thread([&q]
        {
            while(true)
            {
                auto to_push = rand();
                q.push_back(to_push);
                debug_stream << "pushed " << to_push << "\n";
                std::this_thread::sleep_for(std::chrono::milliseconds{10});
            }
        });

    auto thr2 = std::thread([&q]
        {
            while(true)
            {
                debug_stream << "popped " << q.pop_front() << "\n";
                std::this_thread::sleep_for(std::chrono::milliseconds{40});
            }
        });

    (thr1.join(), thr2.join());

    //std::async(std::launch::async, [&task]{while(true) {task("test02"); std::this_thread::sleep_for(std::chrono::milliseconds(10));});


    //std::cout << "task_bind:\t" << result.get() << '\n';

}
/*
template<typename F, typename... Rest>
auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))> {
    auto pck = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
            std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
    );

    auto _f = new std::function<void(int id)>([pck](int id) {
        (*pck)(id);
    });
    this->q.push(_f);

    std::unique_lock<std::mutex> lock(this->mutex);
    this->cv.notify_one();

    return pck->get_future();




std::vector<std::vector<dna4>> queries = input_generator<dna4, 1234>::generate_queries(100, 6);

// state not reset so new text everytime
auto text = input_generator<dna4, 1234>::generate_sequence(1e3);

auto da_index = make_kmer_index<true, 5, 6, 7>(text);
//auto map_index = make_kmer_index<false, 5, 6, 7>(text);
//auto fm = fm_index(text);

// search in paralell
for (size_t i = 0; i < 10; ++i)
{
auto paralell_results = da_index.search(queries);
}


// search in sequence
std::vector<std::vector<unsigned int>> sequence_results;
size_t in = 0;
for (size_t i = 0; i < 10; ++i)
{
debug_stream << "seq search " << std::to_string(i) << "\n";
for (auto &q : queries) {
da_index.search(q);
}
}

debug_stream << "seq done.\n";

//std::cout << (paralell_results == sequence_results ? "test passed" : "test failed");

bool results_equal;
try
{
    //auto fm_results = search(query, fm);
     auto size = da_index.search(query).size() + map_index.search(query).size() + fm_results.size();
     results_equal = (size - (fm_results.size() * 3)) == 0;
}
catch (std::out_of_range)
{
    results_equal = false;
    std::cerr << "out of range exception\n";
}

if (not results_equal)
{
    debug_stream << "results not equal for i = " << i << "\n";
    //debug_stream << text << "\n\n";
    debug_stream << "fm  : " << search(query, fm) << "\n";
    debug_stream << "da  : " << da_index.search(query) << "\n";
    debug_stream << "map : " << map_index.search(query) << "\n";

    return 1;
}

debug_stream << "test passed succesfully.";*/