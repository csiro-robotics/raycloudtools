#include "raycloud.h"
#include "raytreegen.h"
#include "rayforestgen.h"
#include "rayroomgen.h"
#include "rayterraingen.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Generates simple example ray clouds" << endl;
  cout << "usage:" << endl;
  cout << "raycreate room 3 - generates a room using the seed 3. Also:" << endl;
  cout << "          building" << endl;
  cout << "          tree" << endl;
  cout << "          forest" << endl;
  cout << "          terrain" << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();
  
  string type = argv[1];
  int seed = stoi(argv[2]);
  srand(seed);

  Cloud cloud;
  if (type == "room")
  {
    // create room
    RoomGen roomGen;
    roomGen.generate();
    cloud.starts = roomGen.rayStarts;
    cloud.ends = roomGen.rayEnds;
    double time = 0.0;
    double timeDelta = 0.01;
    for (int i = 0; i<(int)cloud.starts.size(); i++)
    {
      cloud.times.push_back(time);
      if (i == (int)cloud.starts.size()/2)
        time += 0.5;
      time += timeDelta;
    }
  }
  else if (type == "building")
  {
    // create building...
    cout << "Sorry, building generation not implemented yet" << endl;
  }
  else if (type == "tree" || type == "forest")
  {
    fillBranchAngleLookup();
    double density = 500.0;
    Vector3d boxMin(-2.0, -2.0, -0.025), boxMax(2.0,2.0,0.025);
    double time = 0.0;
    double timeDelta = 0.01;
    if (type == "tree")
    {
      TreeGen treeGen;
      treeGen.make(Vector3d(0,0,0), 0.1, 0.25);
      treeGen.generateRays(density);
      cloud.starts = treeGen.rayStarts;
      cloud.ends = treeGen.rayEnds;
      for (int i = 0; i<(int)cloud.starts.size(); i++)
      {
        cloud.times.push_back(time);
        time += timeDelta;
      }
    }
    else if (type == "forest")
    {
      ForestGen forestGen;
      forestGen.make(0.25);
      forestGen.generateRays(density);
      for (auto &tree: forestGen.trees)
      {
        cloud.starts.insert(cloud.starts.end(), tree.rayStarts.begin(), tree.rayStarts.end());
        cloud.ends.insert(cloud.ends.end(), tree.rayEnds.begin(), tree.rayEnds.end());
        for (int i = 0; i<(int)tree.rayEnds.size(); i++)
        {
          cloud.times.push_back(time);
          time += timeDelta;
        }
      }
      boxMin *= 2.5;
      boxMax *= 2.5;
    }
    int num = 0.25*density*(boxMax[0]-boxMin[0])*(boxMax[1]-boxMin[1]);
    for (int i = 0; i<num; i++)
    {
      Vector3d pos(random(boxMin[0], boxMax[0]), random(boxMin[1], boxMax[1]), random(boxMin[2], boxMax[2]));
      cloud.ends.push_back(pos);
      cloud.starts.push_back(pos + Vector3d(random(-0.1,0.1), random(-0.1,0.1), random(0.2,0.5)));
      cloud.times.push_back(time);
      time += timeDelta;
    }
  }
  else if (type == "terrain")
  {
    TerrainGen terrain;
    terrain.generate();
    cloud.starts = terrain.rayStarts;
    cloud.ends = terrain.rayEnds;
    double time = 0.0;
    double timeDelta = 0.01;
    for (int i = 0; i<(int)cloud.starts.size(); i++)
    {
      cloud.times.push_back(time);
      time += timeDelta;
    }
  }
  else
    usage();
  cloud.save(type + ".ply");

  return true;
}
