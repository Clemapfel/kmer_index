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

class thread_pool;

// single consumer single producer queue
template<typename T>
struct thread_safe_queue
{
    friend class thread_pool;

    public:
        thread_safe_queue() = default;

        ~thread_safe_queue()
        {
            abort();
        };

        void push_back(T t)
        {
            if (_currently_aborting)
                return;

            {
                std::lock_guard<std::mutex> lock(_insert_mutex);
                _data.push_back(t);
            }
        }

        template<typename... args_t>
        void emplace_back(args_t... args)
        {
            if (_currently_aborting)
                return;
            {
                std::lock_guard<std::mutex> lock(_insert_mutex);
                _data.emplace_back(args...);
            }
        }

        std::optional<T> pop_front()
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);

            if (_currently_aborting || _data.empty())
                return std::move(std::optional<T>(std::nullopt));

            auto output = std::move(_data.front());
            _data.pop_front();
            return std::move(output);
        }

        bool empty() const
        {
            return _data.empty();
        }


        template<typename function_t>
        void apply(function_t&& f)
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            for (auto& e : _data)
                std::apply(std::forward<function_t>(f), e);
        }

    protected:
        std::mutex& get_mutex()
        {
            return _insert_mutex;
        }

        void abort()
        {
            std::lock_guard<std::mutex> lock(_insert_mutex);
            _currently_aborting = true;
            _data.clear();
        }

    private:
        std::deque<T> _data{};
        std::mutex _insert_mutex;

        std::condition_variable _cv;

        bool _currently_aborting = false;
};

struct thread_pool
{
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

        // ###########################

        // queue for storing tasks (tasks = wrapper std::packaged_task)
        thread_safe_queue<std::unique_ptr<_task_wrapper_base>> _task_queue;

        std::condition_variable _task_cv;
        std::mutex _task_mutex;

        std::vector<std::thread> _threads;

        bool _currently_aborting = false;

        void setup_threads(size_t);

    public:
        ~thread_pool();
        thread_pool(size_t n_threads);

        void resize(size_t n_threads);

        template<typename function_t, typename... args_t>
        auto execute(function_t&& f, args_t... args) -> std::future<std::invoke_result_t<function_t, args_t...>>&&
        {
            std::packaged_task<std::invoke_result_t<function_t, args_t...>()> to_wrap(std::bind(f, args...));

            auto future = to_wrap.get_future();
            _task_queue.emplace_back(wrap([task(std::move(to_wrap))]() mutable {task();}));
            // move into lambda scope, wrap and store in queue
            // lambda needs to be mutable to modify task even though it's captured by value (moved)

            _task_cv.notify_one();

            return std::move(future);
        }

        void wait_to_finish();

};

// references:
// https://codereview.stackexchange.com/questions/221626/c17-thread-pool
// https://github.com/vit-vit/ctpl
// https://livebook.manning.com/book/c-plus-plus-concurrency-in-action/chapter-9/17
