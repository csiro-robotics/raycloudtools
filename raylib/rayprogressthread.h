// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYPROGRESSTHREAD_H
#define RAYPROGRESSTHREAD_H

#include "raylib/raylibconfig.h"

#include "rayprogress.h"

#include <atomic>
#include <thread>

namespace ray
{
/// A utility class used to display @c Progess to @c std::cout .
///
/// This class is inteded for use from a @c main() entry point. Typical usage is to create a @c Progess object followed
/// by a @c ProgressThread in stack memory to start progress display. Before exit, @c join() should be called.
/// @c requestQuit() may be used to exit the thread early.
class RAYLIB_EXPORT ProgressThread
{
public:
  /// Initialise a new progress thread to the given @p progress tracker.
  /// The thread starts immediately.
  ProgressThread(Progress &progress);

  /// Destructor ensuring the thread is joined.
  ~ProgressThread();

  void requestQuit() { quit_flag_ = true; }
  void join();

private:
  void run();

  void showProgress(Progress &progress, bool finalise, const Progress::Clock::duration *duration) const;

  Progress &progress_;
  std::atomic_bool quit_flag_;
  std::atomic_bool running_;
  std::thread thread_;
};
}  // namespace ray

#endif  // RAYPROGRESSTHREAD_H
