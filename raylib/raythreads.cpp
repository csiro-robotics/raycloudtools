// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "raythreads.h"

#include <memory>
#include <optional>

#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/info.h>
#endif  // RAYLIB_WITH_TBB

using namespace ray;

namespace
{
#if RAYLIB_WITH_TBB
/// TODO(TH)
std::optional<int> initial_thread_count{ std::nullopt };

/// TODO(TH)
std::unique_ptr<oneapi::tbb::global_control> thread_count_global_control{ nullptr };

/// Store the current value of
/// @c oneapi::tbb::global_control::max_allowed_parallelism in
/// @c initial_thread_count .
/// This must be run before any call to
/// @c oneapi::tbb::global_control(oneapi::tbb::global_control::max_allowed_parallelism)
/// in order to get the unconstrained value.
void InitialiseInitialThreadCount()
{
  if (!initial_thread_count.has_value() && (thread_count_global_control == nullptr))
  {
    initial_thread_count =
      oneapi::tbb::global_control::active_value(oneapi::tbb::global_control::max_allowed_parallelism);
  }
}

#endif  // RAYLIB_WITH_TBB

}  // namespace

int Threads::availableThreads()
{
#if RAYLIB_WITH_TBB
  InitialiseInitialThreadCount();
  return initial_thread_count.value();
#else   // RAYLIB_WITH_TBB
  return 1;  // Single threaded.
#endif  // RAYLIB_WITH_TBB
}

int Threads::recommendedThreadCount()
{
#if RAYLIB_WITH_TBB
  InitialiseInitialThreadCount();
  // Thread performance seems to peek between 4-6. For optimal threads, we use
  // at least 2 threads (if available) up to 6 threads. We try to leave one
  // thread free and unused for the system and other processes.
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

void Threads::init(const int thread_count)
{
#if RAYLIB_WITH_TBB
  InitialiseInitialThreadCount();
  // Translate to TBB thread count
  int init_thread_count;
  switch (thread_count)
  {
  case ThreadCountAll:
    init_thread_count = availableThreads();
    break;
  case ThreadCountRecommended:
    init_thread_count = recommendedThreadCount();
    break;
  default:
    init_thread_count = std::max(1, thread_count);
    break;
  }
  thread_count_global_control = std::make_unique<oneapi::tbb::global_control>(
    oneapi::tbb::global_control::max_allowed_parallelism, init_thread_count);
#else   // RAYLIB_WITH_TBB
  RAYLIB_UNUSED(thread_count);
#endif  // RAYLIB_WITH_TBB
}
