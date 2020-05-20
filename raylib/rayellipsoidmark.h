// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYELLIPSOIDMARK_H
#define RAYELLIPSOIDMARK_H

#include "raylib/raylibconfig.h"

#if RAYLIB_WITH_TBB
#include <tbb/spin_mutex.h>
#endif  // RAYLIB_WITH_TBB

#include <limits>
#include <vector>

namespace ray
{
/// Helper structure for marking ellipsoids during @c TransientFilter opertaions.
///
/// Tracks the following meta data about an ellipsed:
/// - @c passThroughId() : IDs of the rays passing near the ellipse
/// - @c hits() : the number of rays intersecting the ellipse.
/// - @c firstIntersectionTime() : the time stamp of the first intersecting ray
/// - @c lastIntersectionTime() : the time stamp of the last intersecting ray
///
/// Supports threadsafe access when compiled with @c RAYLIB_WITH_TBB
class RAYLIB_EXPORT EllipsoidMark
{
public:
#if RAYLIB_WITH_TBB
  using Mutex = tbb::spin_mutex;
#endif  // RAYLIB_WITH_TBB

  inline explicit EllipsoidMark(size_t id = 0u)
    : id_(id)
    , first_intersection_time_(std::numeric_limits<double>::max())
    , last_intersection_time_(std::numeric_limits<double>::lowest())
    , hits_(0u)
#if RAYLIB_WITH_TBB
    , lock_(new Mutex)
#endif  // RAYLIB_WITH_TBB
  {}

  inline size_t id() const { return id_; }

  inline double firstIntersectionTime() const { return first_intersection_time_; }
  inline double lastIntersectionTime() const { return last_intersection_time_; }

  inline size_t hits() const { return hits_; }

  inline const std::vector<size_t> &passThroughIds() const { return pass_through_ids_; }

  void reset(size_t id = 0u);

  void sortPassThroughIds();

  void hit(size_t ray_id, double time);
  void passThrough(size_t ray_id);

private:
  std::vector<size_t> pass_through_ids_;
  size_t id_;
  double first_intersection_time_;
  double last_intersection_time_;
  size_t hits_;
#if RAYLIB_WITH_TBB
  std::shared_ptr<Mutex> lock_;
#endif  // RAYLIB_WITH_TBB
};

}  // namespace ray


#endif // RAYELLIPSOIDMARK_H
