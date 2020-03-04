#include "raylaz.h"

#include <liblas/reader.hpp>
#include <liblas/factory.hpp>
#include <liblas/point.hpp>
using namespace std;
using namespace Eigen;
using namespace RAY;

bool RAY::readLas(string fileName, vector<Vector3d> &positions, vector<double> &times, vector<double> &intensities, int decimation)
{
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
      intensities.push_back(point.GetIntensity());
    }
  }
  return true;
}

void RAY::writeLas(string fileName, const vector<Vector3d> &points, 
        const vector<double> &times, const vector<double> &intensities)
{
  cout << "saving LAZ file" << endl;
  
  liblas::Header header;
  header.SetDataFormatId(liblas::ePointFormat1); // Time only

  if (fileName.find(".las") == std::string::npos && fileName.find(".laz") == std::string::npos)
    fileName += ".las";    
  
  if (fileName.find(".laz") != std::string::npos)
    header.SetCompressed(true); 
  
  cout << "Saving points to " << fileName << endl;

  std::ofstream ofs;
  ofs.open(fileName.c_str(), std::ios::out | std::ios::binary);

  const double scale = 1e-4;
  header.SetScale(scale, scale, scale);

  liblas::Writer writer(ofs, header);

  liblas::Point point(&header); 
  point.SetHeader(&header);//TODO HACK Version 1.7.0 does not correctly resize the data. Commit 6e8657336ba445fcec3c9e70c2ebcd2e25af40b9 (1.8.0 3 July fixes it)
  for (unsigned int i = 0; i < points.size(); i++)
  {
    point.SetCoordinates(points[i][0], points[i][1], points[i][2]);
    point.SetIntensity(intensities[i]);
    if (!times.empty()) 
      point.SetTime(times[i]);
    writer.WritePoint(point);
  }
}

