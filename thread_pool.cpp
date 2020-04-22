//
// Created by clem on 4/22/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "thread_pool.hpp"

thread_pool::thread_pool(size_t n_threads)
{

}

// block new executes, free queue
thread_pool::~thread_pool()
{
    _currently_aborting = true;
    _task_cv.notify_all();

}