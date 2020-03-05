#pragma once
#include "rayutils.h"
#include "raycloud.h"

namespace RAY
{

struct Mesh
{
  std::vector<Eigen::Vector3d> vertices;
  std::vector<Eigen::Vector3i> indexList;

  void splitCloud(const RayCloud &cloud, RayCloud &inside, RayCloud &outside);
};


}
