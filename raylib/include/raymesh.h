#pragma once
#include "rayutils.h"
#include "raycloud.h"

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

struct Mesh
{
  std::vector<Eigen::Vector3d> vertices;
  std::vector<Eigen::Vector3i> indexList;

  void splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside);
};


}
