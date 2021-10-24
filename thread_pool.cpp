//
// Created by clem on 4/22/20.
// Copyright (c) 2020 Clemens Cords (mail@clemens-cords.com). All rights reserved.
//

#include <thread_pool.hpp>
#include <cassert>

namespace kmer::detail
{
// create threads
void thread_pool::setup_threads(size_t n_threads)
{
    assert(_threads.empty());

    for (size_t i = 0; i < n_threads; ++i)
    {
        _threads.emplace_back([&]() {
            std::unique_lock<std::mutex> queue_lock{_queue_mutex, std::defer_lock};

            while (true)
            {
                queue_lock.lock();

                _task_cv.wait(queue_lock, [&]() -> bool {
                    return _shutdown_asap || !_task_queue.empty() || _currently_aborting;
                });

                // shutdown : return asap
                if (_shutdown_asap)
                {
                    //debug::sync_print("shutting down");
                    return;
                }

                // shutdown by DTOR: finish task queue first
                if (_currently_aborting and _task_queue.empty())
                {
                    //debug::sync_print("shutting down");
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

// ctor
thread_pool::thread_pool(size_t n_threads)
{
    assert(n_threads > 0);

    _threads.reserve(n_threads);
    setup_threads(n_threads);
}

//dtor
thread_pool::~thread_pool()
{
    //debug::sync_print("starting dtor");

    // work through current queue until empty, then abort
    _currently_aborting = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    //debug::sync_print("finished dtor");
}

// halt execution, reinit threads, then resume with leftover queue
void thread_pool::resize(size_t n_threads)
{
    //debug::sync_print("starting resize");

    _shutdown_asap = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    //debug::sync_print("all_joined");

    _shutdown_asap = false;

    _threads.clear();
    setup_threads(n_threads);
}

// safely abort all threads
void thread_pool::abort()
{
    //debug::sync_print("starting abort");

    _shutdown_asap = true;
    _task_cv.notify_all();

    for (auto& thr : _threads)
        thr.join();

    while (not _task_queue.empty())
        _task_queue.pop();

    _shutdown_asap = false;

}
} // end of namespace kmer::detail
