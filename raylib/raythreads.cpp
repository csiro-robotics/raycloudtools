// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "raythreads.h"

#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <tbb/task_scheduler_init.h>
#endif  // RAYLIB_WITH_TBB

using namespace ray;

namespace
{
#if RAYLIB_WITH_TBB
std::unique_ptr<tbb::task_scheduler_init> scheduler;
#endif  // RAYLIB_WITH_TBB
}  // namespace

int Threads::availableThreads()
{
#if RAYLIB_WITH_TBB
  return tbb::task_scheduler_init::default_num_threads();
#else   // RAYLIB_WITH_TBB
  return 1;  // Single threaded.
#endif  // RAYLIB_WITH_TBB
}


int Threads::recommendedThreadCount()
{
#if RAYLIB_WITH_TBB
  // Thread performance seems to peek between 4-6. For optimal threads, we use at least 2 threads (if available) up to
  // 6 threads. We try to leave one thread free and unused for the system and other processes.
  const int target_thread_count = MaxRecommendedThreads;
  int thread_count = availableThreads();
  if (thread_count > 2)
  {
    thread_count = std::min(thread_count - 1, target_thread_count);
  }
  return thread_count;
#else   // RAYLIB_WITH_TBB
  return 1;  // Single threaded.
#endif  // RAYLIB_WITH_TBB
}


void Threads::init(int thread_count)
{
#if RAYLIB_WITH_TBB
  if (!scheduler)
  {
    // Translate to TBB thread count
    int init_thread_count = tbb::task_scheduler_init::automatic;
    if (thread_count == ThreadCountRecommended)
    {
      init_thread_count = recommendedThreadCount();
    }
    else if (thread_count > 0)
    {
      init_thread_count = thread_count;
    }
    scheduler = std::make_unique<tbb::task_scheduler_init>(init_thread_count);
  }
#else   // RAYLIB_WITH_TBB
  RAYLIB_UNUSED(thread_count);
#endif  // RAYLIB_WITH_TBB
}
