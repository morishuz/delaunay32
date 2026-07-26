// SPDX-License-Identifier: MIT

#include "internal.hpp"

#include <stdexcept>

namespace delaunay32::detail {

WorkerTeam::WorkerTeam(std::size_t thread_count)
    : thread_count_(thread_count) {
    if (thread_count_ == 0) {
        throw std::invalid_argument(
            "worker team needs at least one thread");
    }
    workers_.reserve(thread_count_ - 1);
    try {
        for (std::size_t index = 1; index < thread_count_; ++index) {
            workers_.emplace_back(
                [this, index] { worker_loop(index); });
        }
    } catch (...) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        start_.notify_all();
        for (std::thread& worker : workers_) {
            worker.join();
        }
        throw;
    }
}

WorkerTeam::~WorkerTeam() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stopping_ = true;
    }
    start_.notify_all();
    for (std::thread& worker : workers_) {
        worker.join();
    }
}

void WorkerTeam::run(
    std::size_t active_threads,
    const std::function<void(std::size_t)>& operation) {
    if (active_threads == 0 || active_threads > thread_count_) {
        throw std::invalid_argument(
            "active worker count exceeds worker team size");
    }
    {
        std::lock_guard<std::mutex> lock(mutex_);
        operation_ = operation;
        active_threads_ = active_threads;
        completed_workers_ = 0;
        first_exception_ = nullptr;
        ++generation_;
    }
    start_.notify_all();

    std::exception_ptr caller_exception;
    try {
        operation_(0);
    } catch (...) {
        caller_exception = std::current_exception();
    }

    std::unique_lock<std::mutex> lock(mutex_);
    complete_.wait(lock, [&] {
        return completed_workers_ == active_threads_ - 1;
    });
    const std::exception_ptr worker_exception = first_exception_;
    operation_ = {};
    lock.unlock();
    if (caller_exception != nullptr) {
        std::rethrow_exception(caller_exception);
    }
    if (worker_exception != nullptr) {
        std::rethrow_exception(worker_exception);
    }
}

void WorkerTeam::worker_loop(std::size_t index) {
    std::size_t observed_generation = 0;
    while (true) {
        std::unique_lock<std::mutex> lock(mutex_);
        start_.wait(lock, [&] {
            return stopping_ ||
                   generation_ != observed_generation;
        });
        if (stopping_) {
            return;
        }
        observed_generation = generation_;
        const bool active = index < active_threads_;
        lock.unlock();
        if (!active) {
            continue;
        }

        std::exception_ptr exception;
        try {
            operation_(index);
        } catch (...) {
            exception = std::current_exception();
        }

        lock.lock();
        if (exception != nullptr && first_exception_ == nullptr) {
            first_exception_ = exception;
        }
        if (++completed_workers_ == active_threads_ - 1) {
            complete_.notify_one();
        }
    }
}

}  // namespace delaunay32::detail
