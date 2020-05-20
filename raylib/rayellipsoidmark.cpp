// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "rayellipsoidmark.h"

#include "rayunused.h"

#include <algorithm>

using namespace ray;

void EllipsoidMark::reset(size_t id)
{
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  id_ = id;
  first_intersection_time_ = std::numeric_limits<double>::max();
  last_intersection_time_ = std::numeric_limits<double>::lowest();
  hits_ = 0u;
}

void EllipsoidMark::sortPassThroughIds()
{
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  std::sort(pass_through_ids_.begin(), pass_through_ids_.end());
}

void EllipsoidMark::hit(size_t ray_id, double time)
{
  RAYLIB_UNUSED(ray_id);
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  ++hits_;
  first_intersection_time_ = std::min(time, first_intersection_time_);
  last_intersection_time_ = std::max(time, last_intersection_time_);
}


void EllipsoidMark::passThrough(size_t ray_id)
{
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  pass_through_ids_.emplace_back(ray_id);
}
