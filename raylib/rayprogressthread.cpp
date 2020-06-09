// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "rayprogressthread.h"

#include "rayutils.h"

#include <iostream>
#include <locale>

using namespace ray;

ProgressThread::ProgressThread(Progress &progress)
  : progress_(progress)
  , quit_flag_(false)
  , running_(true)
  , thread_(&ProgressThread::run, this)
{}

ProgressThread::~ProgressThread()
{
  join();
}

void ProgressThread::join()
{
  if (running_)
  {
    quit_flag_ = true;
    thread_.join();
    running_ = false;
  }
}

void ProgressThread::run()
{
  ray::Progress last;
  ray::Progress current;
  progress_.read(&last);

  std::locale mylocale("");
  std::cout.imbue(mylocale);

  while (!quit_flag_)
  {
    progress_.read(&current);
    if (current.phase() != last.phase() || current.target() != last.target())
    {
      // Ensure we finalise the display.
      auto duration = current.lastDuration();
      showProgress(last, true, &duration);
    }

    if (current.progress() != last.progress() || current.target() != last.target())
    {
      showProgress(current, false, nullptr);
      current.read(&last);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }

  // Past update.
  progress_.read(&current);
  // Do not finalise in case we didn't go full progress.
  auto duration = current.lastDuration();
  showProgress(last, false, &duration);
  std::cout << std::endl;
}


void ProgressThread::showProgress(Progress &progress, bool finalise, const Progress::Clock::duration *duration) const
{
  if (finalise)
  {
    // Finalise progress == target.
    if (progress.target())
    {
      progress.setProgress(progress.target());
    }
  }

  if (progress.phase().length() || progress.target() || progress.progress())
  {
    std::cout << "\r                                    \r";
    std::cout << progress.phase() << ' ' << progress.progress();
    if (size_t target = progress.target())
    {
      std::cout << " / " << target;
    }

    if (duration)
    {
      // Log duration with streaming operator from ray
      std::cout << ' ';
      logDuration(std::cout, *duration);
    }

    if (finalise)
    {
      std::cout << std::endl;
    }
    else
    {
      std::cout << std::flush;
    }
  }
}
