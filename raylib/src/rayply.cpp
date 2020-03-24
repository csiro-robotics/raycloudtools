// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayply.h"
using namespace std;
using namespace Eigen;
using namespace RAY;

Vector3d redGreenBlue(double x)
{
  Vector3d res;
  res[0] = 1.0 - x;
  res[1] = 3.0*x*(1.0-x);
  res[2] = x;
  return res;
}

void redGreenBlueGradient(const vector<double> &values, vector<uint32_t> &gradient)
{
  gradient.resize(values.size());
  vector<Vector3d> rgb(values.size());
  for (int i= 0; i<(int)rgb.size(); i++)
    rgb[i] = redGreenBlue(fmod(values[i], 10.0)/10.0);
  for (unsigned int i = 0; i<rgb.size(); i++)
  {
    gradient[i] = 0;
    for (int j = 0; j<3; j++)
      gradient[i] += uint32_t(rgb[i][j] * 255.0) << (j*8);
    gradient[i] += 255<<24; // i.e. alpha is 255
  }
}

// Save the polygon file to disk
void RAY::writePly(const string &fileName, const vector<Vector3d> &starts, const vector<Vector3d> &ends, const vector<double> &times, const vector<double> &intensities, const vector<unsigned int> &colours)
{
  cout << "saving to " << fileName << ", " << ends.size() << " rays." << endl;
  
  vector<uint32_t> RGB(times.size());
  if (colours.size() > 0)
    RGB = colours;
  else
    redGreenBlueGradient(times, RGB);

  vector<Matrix<float, 10, 1> > vertices(ends.size()); // 4d to give space for colour
  for (unsigned int i = 0; i<ends.size(); i++)
  {
    Vector3d n = starts[i] - ends[i];
    union U // TODO: this is nasty, better to just make vertices an unsigned char vector
    {
      float f[2];
      double d;
    };
    U u;
    u.d = times[i];
    vertices[i] << (float)ends[i][0], (float)ends[i][1], (float)ends[i][2], (float)u.f[0], (float)u.f[1], (float)n[0], (float)n[1], (float)n[2], (float &)RGB[i], (float)intensities[i];
  }

  FILE *fid = fopen(fileName.c_str(), "w+");

  // do we put the starts in as normals?
  fprintf(fid, "ply\n"); 
  fprintf(fid, "format binary_little_endian 1.0\n"); 
  fprintf(fid, "comment SDK generated\n"); // TODO: add version here
  fprintf(fid, "element vertex %d\n", int(vertices.size())); 
  fprintf(fid, "property float x\n"); 
  fprintf(fid, "property float y\n"); 
  fprintf(fid, "property float z\n"); 
  fprintf(fid, "property double time\n"); 
  fprintf(fid, "property float nx\n"); 
  fprintf(fid, "property float ny\n"); 
  fprintf(fid, "property float nz\n"); 
  fprintf(fid, "property uchar red\n"); 
  fprintf(fid, "property uchar green\n"); 
  fprintf(fid, "property uchar blue\n"); 
  fprintf(fid, "property uchar alpha\n"); 
  fprintf(fid, "property float intensity\n"); 
  fprintf(fid, "element face 0\n"); 
  fprintf(fid, "property list uchar int vertex_indices\n"); 
  fprintf(fid, "end_header\n"); 

  fwrite(&vertices[0], sizeof(Matrix<float, 10, 1>), vertices.size(), fid);

  fclose(fid);  
}

bool RAY::readPly(const string &fileName, vector<Vector3d> &starts, vector<Vector3d> &ends, vector<double> &times, vector<double> &intensities, vector<uint32_t> &colours)
{
  ifstream input(fileName.c_str());
  if (!input.is_open())
  {
    cerr << "Couldn't open file: " << fileName << endl;
    return false;
  }
  string line;
  int rowSize = 0;
  int offset = -1, normalOffset = -1, intensityOffset = -1, timeOffset = -1, colourOffset = -1;
  bool timeIsFloat = false;
  while (line != "end_header\r" && line != "end_header")
  {
    getline(input, line);
    if (line.find("property float x") != string::npos)
      offset = rowSize;
    if (line.find("property float nx") != string::npos)
      normalOffset = rowSize;
    if (line.find("property float intensity") != string::npos)
      intensityOffset = rowSize;
    if (line.find("time") != string::npos)
    {
      timeOffset = rowSize;
      if (line.find("float") != string::npos)
        timeIsFloat = true;
    }
    if (line.find("property uchar red") != string::npos)
      colourOffset = rowSize;
    if (line.find("float") != string::npos)
      rowSize += sizeof(float);
    if (line.find("double") != string::npos)
      rowSize += sizeof(double);
    if (line.find("property uchar") != string::npos)
      rowSize += sizeof(unsigned char);
  }
  if (offset == -1)
  {
    cerr << "could not find position properties of file: " << fileName << endl;
    return false;
  }
  if (normalOffset == -1)
  {
    cerr << "could not find normal properties of file: " << fileName << endl;
    return false;
  }

  streampos start = input.tellg();
  input.seekg(0, input.end);
  int length = input.tellg() - start;
  input.seekg(start);
  int size = length / rowSize;
  vector<unsigned char> vertices(length);
  // read data as a block:
  input.read((char *)&vertices[0], length);
  bool warningIssued = false;
  int numBounded = 0;
  int numUnbounded = 0;
  for (int i = 0; i<size; i++)
  {
    Vector3f b = (Vector3f &)vertices[rowSize*i + offset];
    Vector3d end(b[0],b[1],b[2]);
    ends.push_back(Vector3d(end));
    if (timeOffset != -1)
    {
      double time;
      if (timeIsFloat)
        time = (double)((float &)vertices[rowSize*i + timeOffset]);
      else
        time = (double &)vertices[rowSize*i + timeOffset];
      if (times.size() > 0 && times.back() > time && !warningIssued)
      {
        cout << "warning, ray times are not in order. Many functions assume chronological ordering" << endl;
        cout << "i: " << i << ", offset: " << timeOffset << ", timeIsFloat" << timeIsFloat << ", last time: " << times.back() << ", this time: " << time << endl;
        warningIssued = true;
      }
      times.push_back(time);
    }
    Vector3f n = (Vector3f &)vertices[rowSize*i + normalOffset];
    starts.push_back(end + Vector3d(n[0], n[1], n[2]));
    if (intensityOffset != -1)
    {
      float intensity = (float &)vertices[rowSize*i + intensityOffset];
      intensities.push_back(intensity);
      if (intensity > 0)
        numBounded++;
      else
        numUnbounded++;
    }
    if (colourOffset != -1)
    {
      uint32_t colour = (uint32_t &)vertices[rowSize*i + colourOffset];
      colours.push_back(colour);
    }
  }
  if (intensities.size() == 0)
  {
    cout << "warning: no intensity information found in " << fileName << ", setting to 1, and assuming no out of range rays" << endl;
    intensities.resize(ends.size());
    for (auto &intensity: intensities)
      intensity = 1.0;
  }
  if (times.size() == 0)
  {
    cout << "warning: no time information found in " << fileName << ", setting times at 1 second intervals per ray" << endl;
    times.resize(ends.size());
    for (int i = 0; i<(int)times.size(); i++)
      times[i] = (double)i;
  }
  if (colours.size() == 0)
  {
    cout << "warning: no colour information found in " << fileName << ", setting colours red->green->blue based on time" << endl;
    redGreenBlueGradient(times, colours);
  }
  cout << "reading from " << fileName << ", " << size << " rays, of which " << numBounded << " bounded and " << numUnbounded << " unbounded" << endl;
  return true; 
}

void RAY::writePlyMesh(const string &fileName, const Mesh &mesh, bool flipNormals)
{
  cout << "saving to " << fileName << ", " << mesh.vertices.size() << " vertices." << endl;
  
  vector<Vector4f> vertices(mesh.vertices.size()); // 4d to give space for colour
  for (unsigned int i = 0; i<mesh.vertices.size(); i++)
    vertices[i] << (float)mesh.vertices[i][0], (float)mesh.vertices[i][1], (float)mesh.vertices[i][2], 1.0;


  FILE *fid = fopen(fileName.c_str(), "w+");
  fprintf(fid, "ply\n"); 
  fprintf(fid, "format binary_little_endian 1.0\n"); 
  fprintf(fid, "comment SDK generated\n"); // TODO: add version here
  fprintf(fid, "element vertex %d\n", int(vertices.size())); 
  fprintf(fid, "property float x\n"); 
  fprintf(fid, "property float y\n"); 
  fprintf(fid, "property float z\n"); 
  fprintf(fid, "property uchar red\n"); 
  fprintf(fid, "property uchar green\n"); 
  fprintf(fid, "property uchar blue\n"); 
  fprintf(fid, "property uchar alpha\n"); 
  fprintf(fid, "element face %d\n", (int)mesh.indexList.size()); 
  fprintf(fid, "property list int int vertex_indices\n"); 
  fprintf(fid, "end_header\n"); 

  fwrite(&vertices[0], sizeof(Vector4f), vertices.size(), fid);

  vector<Vector4i> triangles(mesh.indexList.size());
  if (flipNormals)
    for (int i = 0; i<(int)mesh.indexList.size(); i++)
      triangles[i] = Vector4i(3, mesh.indexList[i][2], mesh.indexList[i][1], mesh.indexList[i][0]);
  else
    for (int i = 0; i<(int)mesh.indexList.size(); i++)
      triangles[i] = Vector4i(3, mesh.indexList[i][0], mesh.indexList[i][1], mesh.indexList[i][2]);
  fwrite(&triangles[0], sizeof(Vector4i), triangles.size(), fid);
  fclose(fid);   
}


bool RAY::readPlyMesh(const string &file, Mesh &mesh)
{
  ifstream input(file.c_str());
  if (!input.is_open())
  {
    cerr << "Couldn't open file: " << file << endl;
    return false;
  }
  string line;
  int rowSize = 0;
  int numberOfFaces = 0;
  int numberOfVertices = 0;
  char char1[100], char2[100];
  while (line != "end_header\r" && line != "end_header")
  {
    getline(input, line);
    if (line.find("float") != string::npos)
    {
      rowSize += 4;
    }
    if (line.find("property uchar") != string::npos)
    {
      rowSize ++;
    }
    if (line.find("element vertex") != string::npos)
      sscanf(line.c_str(), "%s %s %d", char1, char2, &numberOfVertices);
    if (line.find("element face") != string::npos)
      sscanf(line.c_str(), "%s %s %d", char1, char2, &numberOfFaces);    
  }
  
  vector<Vector4f> vertices(numberOfVertices);
  // read data as a block:
  input.read((char *)&vertices[0], sizeof(Vector4f) * vertices.size());
  vector<Vector4i> triangles(numberOfFaces);
  input.read((char *)&triangles[0], sizeof(Vector4i) * triangles.size());
  
  mesh.vertices.resize(vertices.size());
  for (int i = 0; i<(int)vertices.size(); i++)
    mesh.vertices[i] = Vector3d(vertices[i][0], vertices[i][1], vertices[i][2]);
 
  mesh.indexList.resize(triangles.size());
  for (int i = 0; i<(int)triangles.size(); i++)
    mesh.indexList[i] = Vector3i(triangles[i][1], triangles[i][2], triangles[i][3]);
  cout << "reading from " << file << ", " << mesh.indexList.size() << " triangles." << endl;
  return true;
}


