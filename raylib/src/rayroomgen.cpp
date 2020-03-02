#include "rayroomgen.h"
using namespace Ray;
using namespace std;
using namespace Eigen;

struct Cuboid
{
  Cuboid(const Vector3d &minB, const Vector3d &maxB) : minBound(minB), maxBound(maxB) {}
  Vector3d minBound, maxBound;
}

void rayIntersectBox(const Vector3d &start, const Vector3d &dir, const Cuboid &box, double &depth)
{
  double maxDepth = 0.0;
  Vector3d centre = (box.minBound + box.maxBound)/2.0;
  Vector3d extent = (box.maxBound - box.minBound)/2.0;
  Vector3d toCentre = centre - start;
  for (int ax = 0; ax<3; ax++)
  {
    double d = dir[a]>0.0 ? (toCentre[a] - extent[a])/dir[a] : (toCentre[a] + extent[a])/dir[a];
    maxDepth = max(maxDepth, d);
  }
  if (maxDepth < depth && maxDepth != 0.0)
    depth = maxDepth;
}

void rayIntersectNegativeBox(const Vector3d &start, const Vector3d &dir, const Cuboid &box, double &depth)
{
  // TODO: how on earth do I do this?
}

// A room with a door, window, table and cupboard
void RoomGen::generate()
{
  double pointDensity = 500.0; 
  double roomWidth = random(3.0, 6.0);
  double roomLength = random(3.0, 6.0);
  double roomHeight = random(2.75, 3.0);

  Vector3d floorCentre(random(-10.0, 10.0), random(-10.0, 10.0), random(-10.0, 10.0));
  double roomYaw = random(0.0, 2.0*pi);

  vector<Cuboid> negatives;
  Cuboid room(Vector3d(0,0,0), Vector3d(roomWidth, roomLength, roomHeight));
  negatives.push_back(room);

  double doorWidth = 0.7;
  double doorStart = random(0.0, roomWidth-doorWidth);
  double doorHeight = random(2.2, 2.6);
  Vector3d doorPos = Vector3d(doorStart, -0.2, 0.0);
  Cuboid door(doorPos, doorPos + Vector3d(doorWidth, 0.3, doorHeight));
  negatives.push_back(door);

  Cuboid outsideDoor(Vector3d(-20.0, -20.0, 0.0), Vector3d(20.0, -0.15, 20.0));
  negatives.push_back(outsideDoor);

  double windowWidth = random(0.7, 1.5);
  double windowHeight = random(0.9, 1.4);
  double windowStart = random(0.4, roomLength - windowWidth - 0.4);
  Vector3d windowPos = Vector3d(-0.2, windowStart, 1.0);
  Cuboid window(windowPos, windowPos + Vector3d(0.3, windowWidth, windowHeight));
  negatives.push_back(window);

  Cuboid outsideWindow(Vector3d(-20.0, -20.0, 0.0), Vector3d(-0.15, 20.0, 20.0));
  negatives.push_back(outsideWindow);


  double tableWidth = random(0.5, 1.5);
  double tableLength = random(0.5, 1.5);
  double tableHeight = random(0.5, 1.2);
  Vector3d tablePos(Vector3d(random(0.0, roomWidth - tableWidth - 1.0), Vector3d(0.0, roomLength - tableLength - 1.0), tableHeight);
  Cuboid tableTop(tablePos, tablePos + Vector3d(tableWidth, tableLength, 0.05));
  positives.push_back(tableTop);
  for (int x = 0; x<2; x++)
  {
    for (int y = 0; y<2; y++)
    {
      Vector3d pos(tablePos[0] + 0.05 + (tableWidth-0.15)*(double)x, tablePos[1] + 0.05 + (tableLength-0.15)*(double)y,0.0);
      Cuboid leg(pos, pos+Vector3d(0.05, 0.05, tableHeight));
      positives.push_back(leg);
    }
  }

  double cupboardWidth = random(1.0, 2.9);
  double cupboardStart = random(0.0, roomLength - cupboardWidth);
  double cupboardDepth = random(0.2, 0.8);
  double cupboardHeight = random(1.0, 2.5);
  Vector3d cupboardPos(roomWidth - cupboardDepth, cupboardStart, 0.0);
  Cuboic cupboard(cupboardPos, cupboardPos + Vector3d(cupboardDepth, cupboardWidth, cupboardHeight));
  positives.push_back(cupboard);


  // OK now we ray trace from some random location onto the set of cuboids...
  Vector3d start(random(1.4, roomWidth - 1.4), random(1.4, roomLength - 1.4), random(1.3, 2.0));
  Vector3d s = start - Vector3d(roomWidth, roomLength, 0.0)/2.0;
  Vector3d rayStart = Vector3d(sin(s[0]) + cos(s[1]), cos(s[0]) - sin(s[1]), 0.0) + floorCentre;
  int numRays = pointDensity * (roomWidth*roomLength)*2.0 + roomWidth*roomHeight*2.0 + roomLength*roomHeight*2.0;
  for (int i = 0; i<numRays; i++)
  {
    int face = rand()%3;
    Vector3d dir(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0));
    dir.normalize();
    const double maxRange = 20.0;
    double range = maxRange;
    for (auto &cuboid: positives)
      rayIntersectBox(cuboid, range);
    for (auto &cuboid: negatives)
      rayIntersectNegativeBox(cuboid, range);

    Vector3d end = start + range * dir;
    Vector3d s = end - Vector3d(roomWidth, roomLength, 0.0)/2.0;
    Vector3d rayEnd = Vector3d(sin(s[0]) + cos(s[1]), cos(s[0]) - sin(s[1]), 0.0) + floorCentre;
    rayStarts.push_back(rayStart);
    rayEnds.push_back(rayEnd);
  }
}

