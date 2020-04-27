//
// Created by clem on 4/22/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "thread_pool.hpp"
#include <cassert>

void thread_pool::setup_threads(size_t n_threads)
{
    assert(_threads.empty());

    for (size_t i = 0; i < n_threads; ++i)
    {
        _threads.emplace_back([&]()
            {
                std::unique_lock<std::mutex> queue_lock{_queue_mutex, std::defer_lock};

                while (true)
                {
                    queue_lock.lock();

                    std::unique_lock<std::mutex> cv_lock(_task_mutex, std::defer_lock);
                    _task_cv.wait(cv_lock, [&]() -> bool {
                      return !_task_queue.empty() || _currently_aborting;
                    });

                    // shutdown by DTOR
                    sync_print("task queue: " + std::to_string(_task_queue.size()));
                    if (_currently_aborting and _task_queue.empty())
                    {
                        sync_print("shutting down");
                        return;
                    }

                    // grab task and execute
                    auto task = std::move(_task_queue.front()); // transfers ownership
                    _task_queue.pop();

                    if (_task_queue.empty())
                        _wait_to_finish_cv.notify_all();

                    queue_lock.unlock();

                    task->operator()();
                }
            });
    }
}

thread_pool::thread_pool(size_t n_threads)
{
    _threads.reserve(n_threads);
    setup_threads(n_threads);
}

// block new executes, work through queue
thread_pool::~thread_pool()
{
    sync_print("starting dtor");

    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    sync_print("finished dtor");
}

// halt execution, reinit threads, then resume with leftover queue
void thread_pool::resize(size_t n_threads)
{
    std::unique_lock<std::mutex> lock{_task_mutex, std::defer_lock};
    lock.lock();

    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    _threads.clear();

    setup_threads(n_threads);

    lock.unlock();
}

void thread_pool::wait_to_finish()
{
   std::unique_lock<std::mutex> cv_lock(_task_mutex, std::defer_lock);
   _wait_to_finish_cv.wait(cv_lock, [&]() -> bool {return _task_queue.empty() || _currently_aborting;});
}