#include <functional>
#include <deque>
//#include <thread_pool.hpp>

#pragma once

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
        void emplace(args_t... args)
        {
            if (_currently_aborting)
                return;
            {
                std::lock_guard<std::mutex> lock(_insert_mutex);
                _data.emplace_back(std::forward<args_t>(args)...);
            }
        }

        std::optional<T> pop()
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
            std::lock_guard<std::mutex> lock(_op_mutex);
            return _data.empty();
        }

        size_t size() const
        {
            std::lock_guard<std::mutex> lock(_op_mutex);
            return _data.size();
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
        mutable std::mutex _op_mutex;

        std::condition_variable _cv;

        bool _currently_aborting = false;
};