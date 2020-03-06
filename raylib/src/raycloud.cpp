#include "raycloud.h"
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
#include <nabo/nabo.h>
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

// starts are required to get the normal the right way around
vector<Vector3d> Cloud::generateNormals(int searchSize)
{
  // simplest scheme... find 3 nearest neighbours and do cross product
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  MatrixXd pointsP(3, ends.size());
  for (unsigned int i = 0; i<ends.size(); i++)
    pointsP.col(i) = ends[i];//col.transpose();
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);//, 0, params);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(searchSize, ends.size());
  dists2.resize(searchSize, ends.size());
  nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0);
  delete nns;

  vector<Vector3d> normals;
  for (unsigned int i = 0; i<ends.size(); i++)
  {
    Matrix3d scatter;
    scatter.setZero();
    Vector3d average(0,0,0);
    for (int j = 0; j<searchSize; j++)
      average += ends[indices(j, i)];
    average /= double(searchSize);
    for (int j = 0; j<searchSize; j++)
    {
      Vector3d offset = ends[indices(j,i)] - average;
      scatter += offset * offset.transpose();
    }

    SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
    ASSERT(eigenSolver.info() == Success); 
    Vector3d eigenValue = eigenSolver.eigenvalues();
    Matrix3d eigenVector = eigenSolver.eigenvectors();

    swap(eigenValue[0], eigenValue[2]);
    Vector3d temp = eigenVector.col(0);
    eigenVector.col(0) = eigenVector.col(2);
    eigenVector.col(2) = temp;

    Vector3d normal = eigenVector.col(2);
    normal.normalize();
    if ((ends[i] - starts[i]).dot(normal) > 0.0)
      normal = -normal;
    normals.push_back(normal);
  }
  return normals;
}


struct Event
{
  Vector3d pos;
  double time;
  bool transient;
};

void Cloud::findTransients(Cloud &transient, Cloud &fixed, double radius, double timeDelta)
{
  cout << "find transients" << endl;
  Vector3d boxMin(1e10,1e10,1e10);
  Vector3d boxMax(-1e10,-1e10,-1e10);
  for (auto &end: ends)
  {
    boxMin = minVector(boxMin, end);
    boxMax = maxVector(boxMax, end);
  }
  double voxelWidth = 2.0*max(radius, 0.1);
  Vector3d diff = (boxMax - boxMin)/voxelWidth;
  Grid<Event *> grid(ceil(diff[0]), ceil(diff[1]), ceil(diff[2]));
  cout << "created grid: " << grid.dims.transpose() << endl;
  
  // now get the point normals, in order to find the sphere centres...
  vector<Vector3d> normals = generateNormals();
  vector<Event> spheres(ends.size());
  for (int i = 0; i<(int)ends.size(); i++)
  {
    spheres[i].pos = ends[i] - normals[i]*radius;
    spheres[i].time = times[i];
  }
  cout << "created normals and spheres" << endl;

  // next populate the grid with these sphere centres
  for (auto &sphere: spheres)
  {
    Vector3d rad(radius,radius,radius);
    Vector3d bMin = (sphere.pos - rad - boxMin)/voxelWidth;
    Vector3d bMax = (sphere.pos + rad - boxMin)/voxelWidth;
    for (int x = (int)bMin[0]; x<=(int)bMax[0]; x++)
      for (int y = (int)bMin[1]; y<=(int)bMax[1]; y++)
        for (int z = (int)bMin[2]; z<=(int)bMax[2]; z++)
          grid.cell(x,y,z).data.push_back(&sphere);
  }
  cout << "grid is populated" << endl;

  // now walk every ray through the grid and mark if transient
  for (int i = 0; i<(int)ends.size(); i++)
  {
    Vector3d dir = ends[i] - starts[i];
    Vector3d dirSign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (starts[i] - boxMin)/voxelWidth;
    Vector3d end = (ends[i] - boxMin)/voxelWidth;
    Vector3i index(start[0], start[1], start[2]);
    Vector3i endIndex(end[0], end[1], end[2]);
    cout << "ray " << i << " walk from " << index.transpose() << " to " << endIndex.transpose() << endl;
    for (;;)
    {
      auto &spheres = grid.cell(index[0], index[1], index[2]).data;
      for (int j = (int)spheres.size()-1; j>=0; j--)
      {
        Event *sphere = spheres[j];
        if (sphere->transient)
        {
          spheres[j] = spheres.back();
          spheres.pop_back();
          continue;          
        }
        if (abs(times[i] - sphere->time)<timeDelta)
          continue;
        // ray sphere intersection
        double d = (sphere->pos - starts[i]).dot(dir)/dir.squaredNorm();
        double dist2 = (sphere->pos - (starts[i] + dir*d)).squaredNorm();
        if (d>0.0 && d<1.0 && dist2<radius*radius) // I require the ray to get half way into the sphere to consider it transient
        {
          cout << "found transient with d: " << d << " distance: " << sqrt(dist2) << endl;
          sphere->transient = true;
          // remove it from the list, for speed.
          spheres[j] = spheres.back();
          spheres.pop_back();
        }
      }
      if (index == endIndex)
        break;
      Vector3d mid = boxMin + voxelWidth*Vector3d(index[0]+0.5, index[1]+0.5, index[2]+0.5);
      Vector3d nextBoundary = mid + 0.5*voxelWidth*dirSign;
      Vector3d delta = nextBoundary - starts[i];
      Vector3d d(delta[0]/dir[0], delta[1]/dir[1], delta[2]/dir[2]);
      if (d[0] < d[1] && d[0]<d[2])
        index[0] += dirSign[0];
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += dirSign[1];
      else
        index[2] += dirSign[2];
    }
  }

  cout << "generating two clouds from flags" << endl;
  // Lastly, generate the new ray clouds from this sphere information
  for (int i = 0; i<(int)spheres.size(); i++)
  {
    if (spheres[i].transient)
    {
      transient.starts.push_back(starts[i]);
      transient.ends.push_back(ends[i]);
      transient.times.push_back(times[i]);
    }
    else
    {
      fixed.starts.push_back(starts[i]);
      fixed.ends.push_back(ends[i]);
      fixed.times.push_back(times[i]);
    }
  }
}
