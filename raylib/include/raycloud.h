#pragma once
#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"

namespace RAY
{
struct Ellipsoid
{
  Eigen::Vector3d pos;
  Eigen::Vector3d vectors[3];
  double time;
  double size;
  bool transient;
};

template<class T> 
struct Grid
{
  Grid(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth)
  {
    this->boxMin = boxMin;
    this->boxMax = boxMax;
    this->voxelWidth = voxelWidth;
    Eigen::Vector3d diff = (boxMax - boxMin)/voxelWidth;
    dims = Eigen::Vector3i(ceil(diff[0]), ceil(diff[1]), ceil(diff[2]));   
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

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
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
  void combine(std::vector<Cloud> &clouds, Cloud &differences, bool maximal);
  void markIntersectedEllipsoids(Grid<Ellipsoid *> &grid, double timeDelta, bool maximal = false);
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids);

protected:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &lazFile, const std::string &trajFile);
};


}
