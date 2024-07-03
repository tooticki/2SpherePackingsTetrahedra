#include <future>  // async
#include <iostream>  // cout
#include <semaphore>  // counting_semaphore
#include <vector>

static const size_t THREAD_POOL_SIZE_DEFAULT{ std::thread::hardware_concurrency() };
static const size_t THREAD_POOL_SIZE_MAX{ std::thread::hardware_concurrency() * 2 };
static const size_t NUM_TASKS_DEFAULT{ 20 };

template <typename F>
void run_tasks(
    F&& f,
    size_t thread_pool_size = THREAD_POOL_SIZE_DEFAULT,
    size_t num_tasks = NUM_TASKS_DEFAULT)
{
    thread_pool_size = std::min(thread_pool_size, THREAD_POOL_SIZE_MAX);

    std::counting_semaphore task_slots(thread_pool_size);
    
    auto futures{ std::vector<std::future<void>>(num_tasks) };
    auto task_results{ std::vector<int>(num_tasks) };

    // We can run thread_pool_size tasks in parallel
    // If all task slots are busy, we have to wait for a task to finish
    for (size_t i{ 0 }; i < num_tasks; ++i)
    {
        // Wait for a task slot to be free
        task_slots.acquire();

        futures[i] = std::async(
            std::launch::async,
            [i, &f, &task_result = task_results[i], &task_slots]() {
                // Execute task
                task_result = std::forward<F>(f)(i);

                // Release the task slot
                task_slots.release();
            }
        );
    }

    // Wait for all the tasks to finish
    for (auto& future : futures) { future.get(); };
    for (auto& result: task_results) { std::cout << result << " "; }
}
