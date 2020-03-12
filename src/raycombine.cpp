#include "raycloud.h"
#include "raymesh.h"
#include "rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Combines multiple ray clouds as a union of the volumes. Outputs the combined cloud and the residual cloud of differences." << endl;
  cout << "usage:" << endl;
  cout << "raycombine raycloud1 raycloud2 ... raycloudN - combines into one cloud" << endl;
  cout << "                               --intersect - combines as an intersection rather than union." << endl;
  cout << "                               --concat    - combines as a simple concatenation, with all rays remaining." << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc < 2)
    usage();
  vector<string> files;
  int numFiles = argc-1;
  bool maximal = false;
  bool concatenate = false;
  if (string(argv[argc-1]) == "--intersect" || string(argv[argc-1]) == "-i")
  {
    numFiles--;
    maximal = true;
  }
  else if (string(argv[argc-1]) == "--concatenate" || string(argv[argc-1]) == "--concat" || string(argv[argc-1]) == "-c")
  {
    numFiles--;
    concatenate = true;
  }
  if (numFiles < 2)
    usage();
  for (int i = 0; i<numFiles; i++)
  {
    files.push_back(string(argv[i+1]));
    ifstream f(files.back().c_str());
    if (!f.good())
    {
      cout << "could not open file: " << files.back() << endl;
      usage();
    }
  }

  string fileStub = files[0];
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);

  vector<Cloud> clouds(files.size());
  for (int i = 0; i<(int)files.size(); i++)
    clouds[i].load(files[i]);

  Cloud combined;
  if (concatenate)
  {
    for (auto &cloud: clouds)
    {
      combined.starts.insert(combined.starts.end(), cloud.starts.begin(), cloud.starts.end());
      combined.ends.insert(combined.ends.end(), cloud.ends.begin(), cloud.ends.end());
      combined.times.insert(combined.times.end(), cloud.times.begin(), cloud.times.end());
    }
  }
  else
  {
    Cloud differences;
    combined.combine(clouds, differences, maximal);
    differences.save(fileStub + "_differences.ply");
  }
  combined.save(fileStub + "_combined.ply");
  return true;
}
