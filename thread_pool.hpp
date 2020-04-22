//
// Created by Clemens Cords on 4/20/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <functional>
#include <deque>
//#include <thread_pool.hpp>

// single consumer single producer queue
template<typename T>
struct thread_safe_queue
{
    public:
        thread_safe_queue() = default;

        ~thread_safe_queue()
        {
            abort();
        };

        void push_back(T t)
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            _data.push_back(t);
        }

        template<typename... args_t>
        void emplace_back(args_t... args)
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            _data.emplace_back(args...);
        }

        std::optional<T> pop_front()
        {
            if (_currently_aborting || _data.empty())
                return std::move(std::optional<T>(std::nullopt));

            std::lock_guard<std::mutex> lock(_insert_mutex);
            auto output = std::move(_data.front());
            _data.pop_front();
            return std::move(output);
        }

        template<typename function_t>
        void apply(function_t&& f)
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            std::apply(f, _data);
        }

    protected:
        void abort()
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            _currently_aborting = true;
            _data.clear();
        }

    private:
        std::deque<T> _data{};
        std::mutex _insert_mutex;

        bool _currently_aborting = false;
};

struct thread_pool
{
    private:
        // callable wrapper for storing functions in thread pool queue
        // abstract class so unique_ptr can hold wrapper with type erasure
        struct _task_wrapper_base
        {
            virtual void operator()() = 0;
        };

        template<typename function_t>
        struct _task_wrapper : protected _task_wrapper_base
        {
            public:
                _task_wrapper(function_t&& function)
                        : _function(std::forward<function_t>(function))
                {}

                void operator()()
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
            return std::make_unique<_task_wrapper_base>(new _task_wrapper<function_t>(std::forward<function_t>(function)));
        }

        // queue for storing tasks (tasks = wrapper std::packaged_task)
        thread_safe_queue<std::unique_ptr<_task_wrapper_base>> _queue;

        bool _currently_aborting = false;
        std::condition_variable _task_cv;

    public:
        ~thread_pool();
        thread_pool(size_t n_threads);

        void resize(size_t n_threads);

        template<typename function_t, typename... args_t>
        auto execute(function_t&& f, args_t... args) -> std::future<std::invoke_result_t<function_t, args_t...>>&&
        {
            std::packaged_task<std::invoke_result_t<function_t, args_t...>()> to_wrap(std::bind(f, args...));

            auto future = task.get_future();
            _queue.emplace_back(wrap([task{std::move(to_wrap)}]() mutable {task();}));
            // move into lambda scope, wrap and store in queue
            // lambda needs to be mutable to modify task even though it's captured by value

            _task_cv.notify_one();

            return std::move(future);
        }

};

// references:
// https://codereview.stackexchange.com/questions/221626/c17-thread-pool
// https://github.com/vit-vit/ctpl
