#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <future>
#include <thread>
#include <queue>
#include <functional>
#include <mutex>

class ThreadPool 
{
public:
    ThreadPool(size_t size)
        : _shutdown(false)
    {
        _threads.reserve(size);
        for (size_t i = 0; i < size; ++i)
            _threads.push_back(std::thread(&ThreadPool::thread_loop, this));
    }
    ~ThreadPool()
    {
        _shutdown = true;
        _cv.notify_all();
        for (std::thread& thread : _threads)
            thread.join();
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
    {
        typedef decltype(f(args...)) return_type;

        auto task = std::make_shared<std::packaged_task<return_type()>>
        (
            std::bind(std::forward<F>(f),
                      std::forward<Args>(args)...)
        );

        std::lock_guard<std::mutex> lock(_mutex); // lock to synchronize

        // place task in queue
        _tasks.emplace([task]()
        { 
            (*task)();
        });

        _cv.notify_one(); // notify a thread that a task is available

        return task->get_future();
    }

private:
    void thread_loop()
    {
        while (true)
        {
            std::unique_lock<std::mutex> lock(_mutex);

            _cv.wait(lock, [this] 
                { 
                    return !_tasks.empty() || _shutdown;
                });

            if (_shutdown)
                break;

            std::function task = _tasks.front();
            _tasks.pop();

            lock.unlock();

            task();
        }
    }

private:
    std::vector<std::thread> _threads;
    std::queue<std::function<void()>> _tasks;
    std::condition_variable _cv;
    std::mutex _mutex; // sync access to task queue
    bool _shutdown;
};

#endif
