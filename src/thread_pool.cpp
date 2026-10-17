#include "mdsim/thread_pool.hpp"

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <sched.h>
#include <unistd.h>
#endif

namespace mdsim {

static void threadYield() {
#ifdef _WIN32
    SwitchToThread();
#else
    sched_yield();
#endif
}

int hardwareConcurrency() {
#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    int n = static_cast<int>(si.dwNumberOfProcessors);
    return n > 0 ? n : 1;
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? static_cast<int>(n) : 1;
#endif
}

ThreadPool::ThreadPool(int nthreads) : nworkers_(nthreads) {
    if (nworkers_ < 1) nworkers_ = 1;
    if (nworkers_ > 256) nworkers_ = 256;

#ifdef _WIN32
    handles_.resize(nworkers_ - 1);
    for (int i = 0; i < nworkers_ - 1; ++i) {
        handles_[i] = CreateThread(nullptr, 0, &ThreadPool::trampoline, this, 0, nullptr);
    }
#else
    threads_.resize(nworkers_ - 1);
    for (int i = 0; i < nworkers_ - 1; ++i) {
        threads_[i] = std::thread(&ThreadPool::trampoline, this);
    }
#endif
}

ThreadPool::~ThreadPool() {
    shutdown_.store(true, std::memory_order_release);
    epoch_.fetch_add(1, std::memory_order_release);

#ifdef _WIN32
    for (int i = 0; i < static_cast<int>(handles_.size()); ++i) {
        if (handles_[i]) {
            WaitForSingleObject(handles_[i], INFINITE);
            CloseHandle(handles_[i]);
        }
    }
#else
    for (auto& t : threads_) {
        if (t.joinable()) t.join();
    }
#endif
}

#ifdef _WIN32
unsigned long __stdcall ThreadPool::trampoline(void* arg) {
    static_cast<ThreadPool*>(arg)->workerLoop();
    return 0;
}
#else
void ThreadPool::trampoline(ThreadPool* p) { p->workerLoop(); }
#endif

void ThreadPool::workerLoop() {
    int seen = epoch_.load(std::memory_order_relaxed);
    for (;;) {
        while (epoch_.load(std::memory_order_acquire) == seen) {
            if (shutdown_.load(std::memory_order_relaxed)) return;
            threadYield();
        }
        seen = epoch_.load(std::memory_order_relaxed);
        if (shutdown_.load(std::memory_order_relaxed)) return;

        for (;;) {
            int b = next_block_.fetch_add(1, std::memory_order_relaxed);
            if (b >= nblocks_) break;
            task_(b);
        }
        // Last worker to finish resets the counters for the next round.
        if (done_count_.fetch_add(1, std::memory_order_acq_rel) == nworkers_ - 1) {
            done_count_.store(0, std::memory_order_relaxed);
        }
    }
}

void ThreadPool::runBlocks(int nblocks, const std::function<void(int)>& fn) {
    task_ = fn;
    nblocks_ = nblocks;
    next_block_.store(0, std::memory_order_relaxed);
    done_count_.store(0, std::memory_order_relaxed);
    epoch_.fetch_add(1, std::memory_order_release);

    // Caller thread participates.
    for (;;) {
        int b = next_block_.fetch_add(1, std::memory_order_relaxed);
        if (b >= nblocks) break;
        fn(b);
    }
    if (done_count_.fetch_add(1, std::memory_order_acq_rel) == nworkers_ - 1) {
        done_count_.store(0, std::memory_order_relaxed);
        return;
    }
    while (done_count_.load(std::memory_order_acquire) != 0) {
        threadYield();
    }
}

} // namespace mdsim
