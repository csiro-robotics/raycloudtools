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
  if (string(argv[argc-1]) == "--max" || string(argv[argc-1]) == "-m")
  {
    numFiles--;
    maximal = true;
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

  vector<Cloud> clouds(files.size());
  for (int i = 0; i<(int)files.size(); i++)
    clouds[i].load(files[i]);

  Cloud combined;
  Cloud differences;
  combined.combine(clouds, differences, maximal);

  string fileStub = files[0];
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);

  combined.save(fileStub + "_combined.ply");
  differences.save(fileStub + "_differences.ply");
  return true;
}
