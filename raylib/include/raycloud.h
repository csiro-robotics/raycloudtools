#pragma once
#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"

namespace RAY
{

struct Cloud
{
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends; 
  std::vector<double> times;
  std::vector<double> intensities;

  void save(const std::string &fileName);
  bool load(const std::string &fileName);
  bool load(const std::string &pointCloud, const std::string &trajFile);

  void transform(const Pose &pose, double timeDelta);
  
protected:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &lazFile, const std::string &trajFile);
};


}
