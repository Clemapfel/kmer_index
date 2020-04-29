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

// thread pool of variable size, accepts tasks and puts them in a not lock-free queue,
// worker threads work through them as they free up.

struct thread_pool
{
    // reference used: https://codereview.stackexchange.com/questions/221626/c17-thread-pool

    private:
        // callable wrapper for storing functions in thread pool queue
        // abstract class so unique_ptr can hold wrapper which means function of all signatures can be held
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

        // parallelization primitives
        std::condition_variable _task_cv;
        std::mutex _queue_mutex;

        // worker threads
        std::vector<std::thread> _threads;

        bool _currently_aborting = false;
        bool _paused_for_resizing = false;

        // fill empty _threads with threads
        void setup_threads(size_t);

    public:
        ~thread_pool();
        thread_pool(size_t n_threads = std::thread::hardware_concurrency());

        // let all threads finisht their current task, then join.
        // Clear all threads then reallocated new n_threads that pickup leftover queue
        void resize(size_t n_threads);

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
            _task_queue.emplace(wrap([task(std::move(to_wrap))]() mutable {task();}));

            lock.unlock();

            _task_cv.notify_one();
            return std::move(future);
        }

        // TODO
        void wait_to_finish();

};

// references:
// https://github.com/vit-vit/ctpl
// https://livebook.manning.com/book/c-plus-plus-concurrency-in-action/chapter-9/17
