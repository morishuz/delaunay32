// SPDX-License-Identifier: MIT

#pragma once

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace delaunay32::detail {

class ThreadBarrier {
public:
    explicit ThreadBarrier(std::size_t participants)
        : participants_(participants) {}

    void wait() {
        std::unique_lock<std::mutex> lock(mutex_);
        const std::size_t generation = generation_;
        if (++arrived_ == participants_) {
            arrived_ = 0;
            ++generation_;
            condition_.notify_all();
            return;
        }
        condition_.wait(lock, [&] { return generation_ != generation; });
    }

private:
    const std::size_t participants_;
    std::size_t arrived_ = 0;
    std::size_t generation_ = 0;
    std::mutex mutex_;
    std::condition_variable condition_;
};

// A triangulation uses the same workers for radix sorting, topology, and
// export. Each run publishes one operation and wakes only the requested
// prefix of worker indices; index zero always runs on the calling thread.
class WorkerTeam {
public:
    explicit WorkerTeam(std::size_t thread_count);
    ~WorkerTeam();

    WorkerTeam(const WorkerTeam&) = delete;
    WorkerTeam& operator=(const WorkerTeam&) = delete;

    void run(
        std::size_t active_threads,
        const std::function<void(std::size_t)>& operation);

private:
    void worker_loop(std::size_t index);

    const std::size_t thread_count_;
    std::vector<std::thread> workers_;
    std::function<void(std::size_t)> operation_;
    std::size_t active_threads_ = 1;
    std::size_t completed_workers_ = 0;
    std::size_t generation_ = 0;
    bool stopping_ = false;
    std::exception_ptr first_exception_;
    std::mutex mutex_;
    std::condition_variable start_;
    std::condition_variable complete_;
};

}  // namespace delaunay32::detail
