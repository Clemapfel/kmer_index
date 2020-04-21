//
// Created by clem on 4/20/20.
//

#pragma once

#include <thread>
#include <future>
#include <functional>
#include <type_traits>

template<size_t n_threads>
class thread_pool
{
    public:
        thread_pool()
        {

        }

        ~thread_pool()
        {
            // wait for all to finish
        }

        template<typename function_t, typename... args_t>
        auto push(function_t&& f, args_t&&... args)
        {
            //std::packaged_task<function_t(args_t...)> task(std::bind(std::forward<function_t>(f), std::forward<args_t>(args)...));
            //auto out = task.get_future();
            //task(args...);

            //return out;
        }

    protected:

    private:
        std::vector<std::thread> _threads;


        std::mutex _mutex;
        std::condition_variable _condition_variable;

};

/*
// ######################################################
// reference: https://github.com/vit-vit/CTPL/blob/master/ctpl.h

template<typename F, typename... Rest>
auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))> {
    auto pck = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
            std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
    );

    auto _f = new std::function<void(int id)>([pck](int id) {
        (*pck)(id);
    });
    this->q.push(_f);

    std::unique_lock<std::mutex> lock(this->mutex);
    this->cv.notify_one();

    return pck->get_future();
}

#endif //RAT_GAME_THREAD_POOL_HPP
*/