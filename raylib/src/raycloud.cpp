#include "raycloud.h"
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
#include <set>

using namespace std;
using namespace Eigen;
using namespace RAY;

void Cloud::save(const std::string &fileName)
{
  string name = fileName;
  if (name.substr(name.length()-4) != ".ply")
    name += ".ply";
  writePly(name, starts, ends, times);
}

bool Cloud::load(const std::string &fileName)
{
  // look first for the raycloud PLY
  if (fileName.substr(fileName.size()-4)==".ply")
    return loadPLY(fileName);
  if(ifstream((fileName + ".ply").c_str(), ios::in))
    return loadPLY(fileName + ".ply");

  // otherwise, look for a .laz and _traj.txt file by that name
  if (ifstream((fileName + ".laz").c_str(), ios::in) && ifstream((fileName + "_traj.txt").c_str(), ios::in))
    return loadLazTraj(fileName + ".laz", fileName + "_traj.txt");

  return false;
}  

bool Cloud::load(const std::string &pointCloud, const std::string &trajFile)
{
  string nameEnd = pointCloud.substr(pointCloud.size()-4);
  if (nameEnd == ".ply")
    readPly(pointCloud, starts, ends, times);
  else if (nameEnd == ".laz" || nameEnd == ".las")
    readLas(pointCloud, ends, times, intensities, 1);    
  else
  {
    cout << "Error converting unknown type: " << pointCloud << endl;
    return false;
  }

  Trajectory trajectory;
  trajectory.load(trajFile);
  calculateStarts(trajectory);
  return true;  
}

bool Cloud::loadPLY(const string &file)
{
  return readPly(file, starts, ends, times);
}

bool Cloud::loadLazTraj(const string &lazFile, const string &trajFile)
{
  bool success = readLas(lazFile, ends, times, intensities, 1);
  if (!success)
    return false;
  Trajectory trajectory;
  trajectory.load(trajFile);
  calculateStarts(trajectory);
  return true;
}

void Cloud::calculateStarts(const Trajectory &trajectory)
{
  if (trajectory.nodes.size()>0)
  {
    int n = 1;
    starts.resize(ends.size());
    for (int i = 0; i<(int)ends.size(); i++)
    {
      while ((times[i] > trajectory.nodes[n].time) && n<(int)trajectory.nodes.size()-1)
        n++;
      double blend = (times[i] - trajectory.nodes[n-1].time)/(trajectory.nodes[n].time-trajectory.nodes[n-1].time);
      starts[i] = trajectory.nodes[n-1].pose.position + (trajectory.nodes[n].pose.position - trajectory.nodes[n-1].pose.position)*clamped(blend, 0.0, 1.0);
    }
  }
  else
    cout << "can only recalculate when a trajectory is available" << endl;
}

void Cloud::transform(const Pose &pose, double timeDelta)
{
  for (int i = 0; i<(int)starts.size(); i++)
  {
    starts[i] = pose * starts[i];
    ends[i] = pose * ends[i];
    times[i] += timeDelta;
  }  
}

struct Vector3iLess
{
  bool operator()(const Vector3i &a, const Vector3i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    return a[2] < b[2];
  }
};

vector<int> voxelRandomSubsample(const vector<Vector3d> &points, double voxelWidth)
{
  vector<int> indices;
  set<Vector3i, Vector3iLess> testSet;
  for (unsigned int i = 0; i<points.size(); i++)
  {
    Vector3i place = Vector3i(floor(points[i][0] / voxelWidth), floor(points[i][1] / voxelWidth), floor(points[i][2] / voxelWidth));
    if (testSet.find(place) == testSet.end())
    {
      testSet.insert(place);
      indices.push_back(i);
    }
  }
  return indices;
}

void Cloud::decimate(double voxelWidth)
{
  vector<int> subsample = voxelRandomSubsample(ends, voxelWidth);
  for (int i = 0; i<(int)subsample.size(); i++)
  {
    int id = subsample[i];
    starts[i] = starts[id];
    ends[i] = ends[id];
    if (intensities.size() > 0)
      intensities[i] = intensities[id];
    times[i] = times[id];
  }
  starts.resize(subsample.size());
  ends.resize(subsample.size());
  if (intensities.size() > 0)
    intensities.resize(subsample.size());
  times.resize(subsample.size());
}