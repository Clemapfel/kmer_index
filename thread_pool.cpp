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

                    _task_cv.wait(queue_lock, [&]() -> bool {
                      return _paused_for_resizing || !_task_queue.empty() || _currently_aborting;
                    });

                    // shutdown by resize()
                    if (_paused_for_resizing)
                    {
                        //sync_print("shutting down");
                        return;
                    }

                    // shutdown by DTOR
                    if (_currently_aborting and _task_queue.empty())
                    {
                        //sync_print("shutting down");
                        return;
                    }

                    // grab task and execute
                    auto task = std::move(_task_queue.front());
                    _task_queue.pop();

                    queue_lock.unlock();

                    task->operator()();
                }
            });
    }
}

thread_pool::thread_pool(size_t n_threads)
{
    assert(n_threads > 0);

    _threads.reserve(n_threads);
    setup_threads(n_threads);
}

thread_pool::~thread_pool()
{
    //sync_print("starting dtor");

    // work through current queue until empty, then abort
    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    //sync_print("finished dtor");
}

// halt execution, reinit threads, then resume with leftover queue
void thread_pool::resize(size_t n_threads)
{
    //sync_print("starting resize");

    _paused_for_resizing = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    //sync_print("all_joined");

    _paused_for_resizing = false;

    _threads.clear();
    setup_threads(n_threads);
}

void thread_pool::wait_to_finish()
{

}