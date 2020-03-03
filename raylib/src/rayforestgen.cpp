#include "rayforestgen.h"
#include "rayutils.h"
using namespace Ray;
using namespace std;
using namespace Eigen;

void ForestGen::make(double randomFactor)
{
  double fieldWidth = 10;
  double maxTreeRadius = 0.1;
  double dimension = 2.0; // #trees = radius^-d
  double adultTreeDensity = 0.04; // #trees per m^2

  double rad = maxTreeRadius;
  double numTrees = sqr(fieldWidth) * adultTreeDensity;
  for (int level = 0; level<2; level++)
  {
    for (int i = 0; i<(int)numTrees; i++)
    {
      double radius = rad * (1.0 + random(-0.25, 0.5)*randomFactor);
      Vector3d root;
      bool found = false;
      while (!found)
      {
        root = fieldWidth*0.5 * Vector3d(random(-1.0, 1.0), random(-1.0, 1.0), 0.0);
        found = true;
        for (auto &tree: trees)
        {
          double d = (root - tree.root).norm();
          if (d < (radius + tree.branches[0].radius) * 10.0)
          {
            found = false;
            break;
          }
        }
      }
      TreeGen tree;
      trees.push_back(tree);
      trees.back().make(root, radius, randomFactor);
    }
    rad /= 2.0;
    numTrees *= pow(2.0, dimension);
  }
}

void ForestGen::generateRays(double rayDensity)
{
  for (auto &tree: trees)
    tree.generateRays(rayDensity);
}

vector<Vector3d> ForestGen::getCanopy()
{
  vector<Vector3d> canopy;
  for (auto &tree: trees)
    canopy.insert(canopy.end(), tree.leaves.begin(), tree.leaves.end());

  return canopy;
}

vector<Vector3d> ForestGen::getPointCloud()
{
  vector<Vector3d> cloud;
  for (auto &tree: trees)
    cloud.insert(cloud.end(), tree.rayEnds.begin(), tree.rayEnds.end());

  return cloud;
}