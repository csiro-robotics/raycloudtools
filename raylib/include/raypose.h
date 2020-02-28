#pragma once
#include "rayutils.h"
#include <string>
namespace RAY
{
// 3d orientation and position
struct Pose
{
  Eigen::Vector3d position;
  Quaterniond rotation;

  Pose() {}
  // Constructors
  Pose(const Eigen::Vector3d &position, const Quat &rotation)
  {
    this->position = position;
    this->rotation = rotation;
  }
  void set(const Eigen::Vector3d &position, const Quat &rotation)
  {
    this->position = position;
    this->rotation = rotation;
  }
  
  /// Operators
  Pose operator *(const Quaterniond &quat) const
  {
    return Pose(position, rotation * quat);
  }
  Eigen::Vector3d operator *(const Eigen::Vector3d &vec) const
  {
    return position + rotation.rotateVector(vec);
  }
  Pose operator *(const Pose &pose) const
  {
    return Pose(position + rotation.rotateVector(pose.position), rotation * pose.rotation);
  }
  Pose operator *(double scale) const
  {
    return Pose(position*scale, rotation * scale);
  }
  Pose operator /(double scale) const
  {
    return Pose(position / scale, rotation / scale);
  }
  Pose operator +(const Pose &pose) const
  {
    return Pose(position + pose.position, rotation + pose.rotation);
  }
  Pose operator -(const Pose &pose) const
  {
    return Pose(position - pose.position, rotation - pose.rotation);
  }
  void operator *=(const Pose &pose) 
  {
    position += rotation.rotateVector(pose.position);
    rotation *= pose.rotation;
  }
  void operator +=(const Pose &pose) 
  {
    position += pose.position;
    rotation += pose.rotation;
  }
  void operator /=(double x) 
  {
    position /= x;
    rotation /= x;
  }
  Pose operator ~() const
  {
    Quat inv = ~rotation;
    return Pose( inv.rotateVector(-position), inv);
  }
  double &operator[](int index)
  {
    ASSERT(index >=0 && index < 7);
    if (index < 3)
      return position[index];
    return rotation[index - 3];
  }
  
  /// Normalise in-place
  void normalise()
  {
    rotation.normalize();
  }
  /// Return normalized pose
  Pose normalised() const
  {
    Pose result = *this;
    result.normalise();
    return result;
  }

  static Pose identity()
  {
    return Pose(Eigen::Vector3d(0,0,0), Quat(1,0,0,0));
  }
  static Pose zero()
  {
    return Pose(Eigen::Vector3d(0,0,0), Quat(0,0,0,0));
  }
};

}