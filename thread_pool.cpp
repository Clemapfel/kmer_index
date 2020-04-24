//
// Created by clem on 4/22/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "thread_pool.hpp"
#include <cassert>

void thread_pool::setup_threads(size_t n_threads)
{
    assert(_threads.empty());

    std::lock_guard<std::mutex> lock{_task_queue.get_mutex()};  // redundancy < safety

    for (size_t i = 0; i < n_threads; ++i)
    {
        _threads.emplace_back([&]()
            {
                while (true)
                {
                    std::unique_lock<std::mutex> lock{_pool_mutex};

                    _task_cv.wait(lock, [&]() -> bool {
                      return !_task_queue.empty() || _currently_aborting;
                    });

                    // shutdown by DTOR
                    if (_currently_aborting and _task_queue.empty())
                      return;

                    // grab task and execute
                    auto task = _task_queue.pop_front();
                    lock.unlock();
                    lock.release();

                    if (task)   // queue returns optional
                      task.value()->operator()();
                }
            });
    }

    sync_print("finished setting up.");

}

thread_pool::thread_pool(size_t n_threads)
{
    _threads.reserve(n_threads);
    setup_threads(n_threads);
}

// block new executes, work through queue
thread_pool::~thread_pool()
{
    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    _task_queue.abort();
}

// halt execution, reinit threads, then resume with leftover queue
void thread_pool::resize(size_t n_threads)
{
    std::unique_lock<std::mutex> lock{_task_queue.get_mutex()};
    lock.lock();

    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    _threads.clear();

    setup_threads(n_threads);

    lock.unlock();
    //lock.release();
}

void thread_pool::wait_to_finish()
{
    // wait until queue is empty
}