#include "raycloud.h"
#include "rayforestgen.h"
#include "rayroomgen.h"

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
  string seed = argv[2];
  srand(seed);

  Cloud cloud;
  if (type == "room")
  {
    // create room
    RoomGen roomGen;
    roomGen.generate();
    cloud.starts = roomGen.rayStarts;
    cloud.ends = roomGen.rayEnds;
  }
  else if (type == "building")
  {
    // create building...
    cout << "Sorry, building generation not implemented yet" << endl;
  }
  else if (type == "tree")
  {
    fillBranchAngleLookup();
    TreeGen treeGen;
    treeGen.make(Vector3d(0,0,0), 0.1, 0.25);
    treeSim.generateRays(500.0);
    cloud.starts = treeGen.rayStarts;
    cloud.ends = treeGen.rayEnds;
    double time = 0.0;
    double timeDelta = 0.01;
    for (int i = 0; i<(int)tree.rayEnds.size(); i++)
    {
      cloud.times.push_back(time);
      time += timeDelta;
    }
  }
  else if (type == "forest")
  {
    fillBranchAngleLookup();
    ForestGen forestGen;
    forestGen.make(0.25);
    forestGen.generateRays(500.0);
    
    double time = 0.0;
    double timeDelta = 0.01;
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
  }
  else
    usage();
  cloud.save(type + ".ply");

  return true;
}
