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

#include "thread_safe_queue.hpp"
#include "syncstream.hpp"

struct thread_pool
{
    // reference used: https://codereview.stackexchange.com/questions/221626/c17-thread-pool

    private:
        // callable wrapper for storing functions in thread pool queue
        // abstract class so unique_ptr can hold wrapper with type erasure
        class _task_wrapper_base
        {
            public:
                virtual ~_task_wrapper_base() {};
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
            return std::unique_ptr<_task_wrapper_base>(new _task_wrapper<function_t>(std::forward<function_t>(function)));
        }

        // queue for storing tasks
        std::queue<std::unique_ptr<_task_wrapper_base>> _task_queue;

        std::condition_variable _task_cv;
        std::mutex _queue_mutex;

        std::condition_variable _pause_cv;

        std::vector<std::thread> _threads;

        bool _currently_aborting = false;
        bool _skip_to_return = false;

        void setup_threads(size_t);

    public:
        ~thread_pool();
        thread_pool(size_t n_threads);

        void resize(size_t n_threads);

        template<typename function_t, typename... args_t>
        auto execute(function_t&& f, args_t... args)// -> std::future<std::invoke_result_t<function_t, args_t...>>&&
        {
            std::packaged_task<std::invoke_result_t<function_t, args_t...>()> to_wrap(std::bind(f, args...));

            auto future = to_wrap.get_future();

            std::unique_lock<std::mutex> lock(_queue_mutex, std::defer_lock);
            lock.lock();

            _task_queue.emplace(wrap([task(std::move(to_wrap))]() mutable {task();}));
            // move into lambda scope, wrap and store in queue
            // lambda needs to be mutable to modify task even though it's captured by value (moved)
            lock.unlock();

            _task_cv.notify_one();
            return std::move(future);
        }

        void wait_to_finish();

};

// references:
// https://github.com/vit-vit/ctpl
// https://livebook.manning.com/book/c-plus-plus-concurrency-in-action/chapter-9/17
