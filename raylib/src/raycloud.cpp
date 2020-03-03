#include "raycloud.h"
#include "raytrajectory.h"
using namespace std;
using namespace RAY;

void Cloud::save(const std::string &fileName)
{
  writePLY(fileName + ".ply", starts, ends, times);
}

bool Cloud::load(const std::string &fileName)
{
  // look first for the raycloud PLY
  if (ifstream(fileName.c_str(), ios::in))
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
  if (pointCloud.substr(pointCloud.size()-4) == ".ply")
    readPly(pointCloud, starts, ends, times);
  else if (pointCloud.substr(pointCloud.size()-4) == ".laz")
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
  readPly(file, starts, ends, times);
}

bool Cloud::loadLazTraj(const string &lazFile, const string &trajFile)
{
  readLas(lazFile, ends, times, intensities, 1);
  Trajectory trajectory;
  trajectory.load(trajFile);
  calculateStarts(trajectory);
}

void Cloud::calculateStarts(const Trajectory &trajectory)
{
  if (trajectory.size()>0)
  {
    int n = 1;
    starts.resize(ends.size());
    for (int i = 0; i<(int)ends.size(); i++)
    {
      while (times[i] > trajectory.nodes[n].time) && n<(int)trajectory.size()-1)
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
  for (int i = 0; i<starts.size(); i++)
  {
    starts[i] = pose * starts[i];
    ends[i] = pose * ends[i];
    times[i] += timeDelta;
  }  
}