// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYTHREADS_H
#define RAYTHREADS_H

#include "raylib/raylibconfig.h"

#include <memory>

namespace ray
{
/// A utility class for initialising the thread pool size.
///
/// Typical usage is to call @c init() at the start of your program. This is optional and it not specified, all
/// available threads will be used.
///
/// Note: some parts of the program may use threads even when Intel TBB is not availabe (via OpenMP).
///
/// @todo Also cap OpenMP thread usage.
class RAYLIB_EXPORT Threads
{
public:
  /// Argument for use with @c init() indicating all @c availableThreads() should be used.
  static const int ThreadCountAll = -1;
  /// Argument for use with @c init() indicating the @c recommendedThreadCount() should be used.
  static const int ThreadCountRecommended = 0;

  /// The maximum number of threads to use for @c recommendedThreadCount() .
  static const int MaxRecommendedThreads = 8;

  /// Returns the number of available threads. When built with Intel TBB, this returns the number of available
  /// processors. Without TBB, this returns 1.
  static int availableThreads();

  /// Query the recommended thread count. This is set at least two threads if available, prefering one less than the
  /// @c availableThreads() up to @c MaxRecommendedThreads threads.
  static int recommendedThreadCount();

  /// Initialise the thread count.
  static void init(int thread_count = ThreadCountRecommended);
};
}  // namespace ray

#endif  // RAYTHREADS_H
