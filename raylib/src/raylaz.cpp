// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylaz.h"
#if defined(USE_LAS)
#include <liblas/reader.hpp>
#include <liblas/factory.hpp>
#include <liblas/point.hpp>
#endif
using namespace std;
using namespace Eigen;
using namespace RAY;

bool RAY::readLas(string fileName, vector<Vector3d> &positions, vector<double> &times, vector<RGBA> &colours, int decimation)
{
#if defined(USE_LAS)
  vector<double> intensities;
  cout << "readLas: filename: " << fileName << endl;
  
  std::ifstream ifs;
  ifs.open(fileName.c_str(), std::ios::in | std::ios::binary);

  if (!ifs.is_open()){
      cerr << "readLas: failed to open stream" << endl;
      return false;
  }
 
  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  liblas::Header header = reader.GetHeader();

  unsigned int size = header.GetPointRecordsCount();
 
  int numIntensities = 0;
  for(unsigned int i=0; i<size; i++)
  {
    reader.ReadNextPoint();
    liblas::Point point = reader.GetPoint();
    //we're downsampling to 1/decimation of the orginal file.
    if((i%decimation)==0)
    {
      Vector3d position;
      position[0] = point.GetX();
      position[1] = point.GetY();
      position[2] = point.GetZ();
      positions.push_back(position);
      times.push_back(point.GetTime());
      double intensity = point.GetIntensity();
      if (intensity > 0.0)
        numIntensities++;
      intensities.push_back(intensity);
    }
  }
  if (numIntensities == 0)
    for (auto &i: intensities)
      i = 1.0;
  redGreenBlueGradient(times, colours);
  for (int i = 0; i<(int)colours.size(); i++) // add intensity into alhpa channel
    colours[i].alpha = 255.0*clamped(intensities[i], 0.0, 1.0);
  
  cout << "loaded " << fileName << " with " << positions.size() << " points" << endl;
  return true;
#else
  cerr << "readLas: cannot read file as USE_LAS not enabled. Enable using: cmake .. -DUSE_LAS=true" << endl;
  return false;  
#endif
}
