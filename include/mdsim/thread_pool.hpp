// Deterministic thread pool.
//
// std::thread / OpenMP are unavailable on the MinGW.org (win32 thread model)
// toolchain used to validate this project, so workers are backed by native
// Win32 threads on Windows and std::thread elsewhere. Scheduling is static:
// callers submit independent blocks and combine results in a fixed order, so
// the number of threads never changes the numerical result.
#pragma once

#include <atomic>
#include <functional>
#include <vector>

#include <cstdint>

#ifndef _WIN32
#include <thread>
#endif

namespace mdsim {

int hardwareConcurrency();

class ThreadPool {
public:
    explicit ThreadPool(int nthreads);
    ~ThreadPool();

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    int size() const { return nworkers_; }

    // Runs fn(block) for block in [0, nblocks). Each block executes exactly
    // once (possibly on any worker). Determinism is the caller's
    // responsibility (independent blocks + fixed-order reduction).
    void runBlocks(int nblocks, const std::function<void(int)>& fn);

private:
    int nworkers_ = 1;

    std::function<void(int)> task_;
    int nblocks_ = 0;

    std::atomic<int> next_block_{0};
    std::atomic<int> done_count_{0};
    std::atomic<int> epoch_{0};
    std::atomic<bool> shutdown_{false};

#ifdef _WIN32
    std::vector<void*> handles_;
#else
    std::vector<std::thread> threads_;
#endif

    void workerLoop();

#ifdef _WIN32
    static unsigned long __stdcall trampoline(void* arg);
#else
    static void trampoline(ThreadPool* p);
#endif
};

} // namespace mdsim
