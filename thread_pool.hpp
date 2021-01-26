//
// Created by Clemens Cords on 4/20/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <queue>
#include <functional>
#include <map>
#include <iostream>

namespace kmer::detail
{

// generic variable sized thread pool
struct thread_pool
{
    private:
        // callable wrapper for storing functions in thread pool queue
        class _task_wrapper_base
        {
            public:
                virtual ~_task_wrapper_base()
                {};

                virtual void operator()() = 0;
        };

        template<typename function_t>
        struct _task_wrapper : public _task_wrapper_base
        {
            public:
                _task_wrapper(function_t&& function)
                        : _function(std::forward<function_t>(function))
                {}

                virtual void operator()() override
                {
                    _function();
                }

            private:
                function_t _function;
        };

        // allocate wrapper
        template<typename function_t>
        static std::unique_ptr<_task_wrapper_base> wrap(function_t&& function)
        {
            return std::unique_ptr<_task_wrapper_base>(
                    new _task_wrapper<function_t>(std::forward<function_t>(function)));
        }

        // queue for storing tasks
        std::queue<std::unique_ptr<_task_wrapper_base>> _task_queue;

        // parallelization primitives
        std::condition_variable _task_cv;
        std::mutex _queue_mutex;

        // worker threads
        std::vector<std::thread> _threads;

        std::atomic<bool> _currently_aborting = false;
        std::atomic<bool> _shutdown_asap = false;

        // create worker threads
        void setup_threads(size_t);

    public:
        // DTOR
        ~thread_pool();

        // CTOR
        explicit thread_pool(size_t n_threads = std::thread::hardware_concurrency());

        // pause execution, reallocate threads and restart to work through leftover queue
        void resize(size_t n_threads);

        // abort all threads as soon as their current execute call is done and clear queue
        void abort();

        // add function call to task queue
        template<typename function_t, typename... args_t>
        auto execute(function_t&& f, args_t... args)
        {
            // bind f and args, then wrap in packaged_task
            std::packaged_task<std::invoke_result_t<function_t, args_t...>()> to_wrap(std::bind(f, args...));

            auto future = to_wrap.get_future();

            std::unique_lock<std::mutex> lock(_queue_mutex, std::defer_lock);
            lock.lock();

            // wrap packaged task operator() call in lambda and wrap lambda in task_wrapper
            // task moves so lambda (and thus _task_queue element) holds ownership of to_wrap
            _task_queue.emplace(wrap([task(std::move(to_wrap))]() mutable { task(); }));

            lock.unlock();

            _task_cv.notify_one();
            return std::move(future);
        }
};

} // end of namespace kmer::detail

// references used:
// https://codereview.stackexchange.com/questions/221626/c17-thread-pool
// https://github.com/vit-vit/ctpl
// https://livebook.manning.com/book/c-plus-plus-concurrency-in-action/chapter-9/17

// #####################################################################################################################

namespace debug
{
    // synchronized console print that preserves order of invokation
    inline std::map<std::thread::id, size_t> _thread_ids;
    inline std::mutex _stream_mutex;

    inline size_t n_messages = 0;

    template<typename T>
    void sync_print(T t)
    {
        std::lock_guard<std::mutex> lock(_stream_mutex);

        auto id = std::this_thread::get_id();

        if (_thread_ids.find(id) == _thread_ids.end())
            _thread_ids.insert(std::make_pair(id, _thread_ids.size()));

        auto out = "[" + std::to_string(_thread_ids.at(id)) + "] " + std::string(t) + "\n";

        // single string cout::operator<< are atomic on windows and unix systems
        std::cout << out;
    }

} // namespace debug

