// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPOSE_H
#define RAYLIB_RAYPOSE_H

#include "raylib/raylibconfig.h"

#include <Eigen/Geometry>
#include <ostream>
#include <string>
#include "rayutils.h"

namespace ray
{
/// A Euclidean transformation class (3d orientation and position).
class RAYLIB_EXPORT Pose
{
public:
  Eigen::Vector3d position;
  Eigen::Quaterniond rotation;

  Pose() {}
  // Constructors
  Pose(const Eigen::Vector3d &position, const Eigen::Quaterniond &rotation)
  {
    this->position = position;
    this->rotation = rotation;
  }
  void set(const Eigen::Vector3d &position, const Eigen::Quaterniond &rotation)
  {
    this->position = position;
    this->rotation = rotation;
  }

  /// Operators
  Pose operator*(const Eigen::Quaterniond &quat) const { return Pose(position, rotation * quat); }
  Eigen::Vector3d operator*(const Eigen::Vector3d &vec) const { return position + rotation * vec; }
  Pose operator*(const Pose &pose) const { return Pose(position + rotation * pose.position, rotation * pose.rotation); }
  void operator*=(const Pose &pose)
  {
    position += rotation * pose.position;
    rotation *= pose.rotation;
  }
  Pose operator~() const
  {
    Eigen::Quaterniond inv = rotation.conjugate();
    return Pose(inv * -position, inv);
  }
  double &operator[](int index)
  {
    ASSERT(index >= 0 && index < 7);
    if (index < 3)
      return position[index];
    if (index == 3)
      return rotation.w();
    else if (index == 4)
      return rotation.x();
    else if (index == 5)
      return rotation.y();
    else
      return rotation.z();
  }

  /// Normalise in-place
  void normalise() { rotation.normalize(); }
  /// Return normalised pose
  Pose normalised() const
  {
    Pose result = *this;
    result.normalise();
    return result;
  }

  static Pose identity() { return Pose(Eigen::Vector3d(0, 0, 0), Eigen::Quaterniond(1, 0, 0, 0)); }
  static Pose zero() { return Pose(Eigen::Vector3d(0, 0, 0), Eigen::Quaterniond(0, 0, 0, 0)); }
};

inline std::ostream &operator<<(std::ostream &os, const Pose &pose)
{
  os << "x,y,z: " << pose.position.transpose() << ", qw,qx,qy,qz: " << pose.rotation.w() << ", " << pose.rotation.x()
     << ", " << pose.rotation.y() << ", " << pose.rotation.z();
  return os;
}
}  // namespace ray

#endif  // RAYLIB_RAYPOSE_H
