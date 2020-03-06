#pragma once
#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"

namespace RAY
{
template<class T> 
struct Grid
{
  Grid(int dimensionX, int dimensionY, int dimensionZ): dims(dimensionX, dimensionY, dimensionZ)
  {
    cells.resize(dims[0]*dims[1]*dims[2]);
  }
  struct Cell
  {
    std::vector<T> data;
  };
  Cell &cell(int x, int y, int z)
  {
    if (x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2])
      return nullCell;
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }

  Eigen::Vector3i dims;
protected:
  std::vector<Cell> cells;
  Cell nullCell;
};

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
  void decimate(double voxelWidth);

  std::vector<Eigen::Vector3d> generateNormals(int searchSize = 16);
  void findTransients(Cloud &transient, Cloud &fixed, double timeDelta);
  
protected:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &lazFile, const std::string &trajFile);
};


}
