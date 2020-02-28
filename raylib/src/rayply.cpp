#include "rayply.h"
using namespace std;
using namespace Eigen;
using namespace RAY;

// Save the polygon file to disk
void writePly(const string &fileName, const vector<Vector3d> &points, const vector<Vector3d> &rgb)
{
  cout << "saving to " << fileName << endl;
  
  vector<uint32_t> RGB(rgb.size());
  for (unsigned int i = 0; i<rgb.size(); i++)
  {
    RGB[i] = 0;
    for (int j = 0; j<3; j++)
      RGB[i] += uint32_t(rgb[i][j] * 255.0) << (j*8);
    RGB[i] += 255<<24; // i.e. alpha is 255
  }

  vector<Vector4f> vertices(points.size()); // 4d to give space for colour
  for (unsigned int i = 0; i<points.size(); i++)
  {
    vertices[i] << (float)points[i][0], (float)points[i][1], (float)points[i][2], (float &)RGB[i];
  }


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
  fprintf(fid, "element face 0\n"); 
  fprintf(fid, "property list uchar int vertex_indices\n"); 
  fprintf(fid, "end_header\n"); 

  fwrite(&vertices[0], sizeof(Vector4f), vertices.size(), fid);

  fclose(fid);  
}

void RAYwritePlyMesh(const string &fileName, const vector<Vector3d> &points, const vector<Vector3i> &indexList, bool flipNormals)
{
  cout << "saving to " << fileName << endl;
  
  vector<Vector4f> vertices(points.size()); // 4d to give space for colour
  for (unsigned int i = 0; i<points.size(); i++)
    vertices[i] << (float)points[i][0], (float)points[i][1], (float)points[i][2], 1.0;


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
  fprintf(fid, "element face %d\n", (int)indexList.size()); 
  fprintf(fid, "property list int int vertex_indices\n"); 
  fprintf(fid, "end_header\n"); 

  fwrite(&vertices[0], sizeof(Vector4f), vertices.size(), fid);

  vector<Vector4i> triangles(indexList.size());
  if (flipNormals)
    for (int i = 0; i<(int)indexList.size(); i++)
      triangles[i] = Vector4i(3, indexList[i][2], indexList[i][1], indexList[i][0]);
  else
    for (int i = 0; i<(int)indexList.size(); i++)
      triangles[i] = Vector4i(3, indexList[i][0], indexList[i][1], indexList[i][2]);
  fwrite(&triangles[0], sizeof(Vector4i), triangles.size(), fid);
  fclose(fid);   
}


bool readPlyMesh(const string &file, vector<Vector3d> &points, vector<Vector3i> &indexList)
{
  cout << "reading from " << file << endl;
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
 //   cout << line << endl;
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
  input.read((char *)&triangles[0], sizeof(Vector4i)*triangles.size());
  
  points.resize(vertices.size());
  for (int i = 0; i<(int)vertices.size(); i++)
    points[i] = Vector3d(vertices[i][0], vertices[i][1], vertices[i][2]);
 
  indexList.resize(triangles.size());
  for (int i = 0; i<(int)triangles.size(); i++)
    indexList[i] = Vector3i(triangles[i][1], triangles[i][2], triangles[i][3]);
  return true;
}

bool RAYreadPly(const string &file, vector<Vector3d> &points, vector<Vector3d> &normals, int &colourOffset)
{
  cout << "reading from " << file << endl;
  ifstream input(file.c_str());
  if (!input.is_open())
  {
    cerr << "Couldn't open file: " << file << endl;
    return false;
  }
  string line;
  int rowSize = 0;
  int normalOffset = 0;
  while (line != "end_header\r" && line != "end_header")
  {
    getline(input, line);
    if (line.find("property float nx") != string::npos)
    {
      normalOffset = rowSize;
    }
    if (line.find("float") != string::npos)
    {
      rowSize += 4;
    }
    if (line.find("property uchar red") != string::npos)
    {
      colourOffset = rowSize;
    }
    if (line.find("property uchar") != string::npos)
    {
      rowSize ++;
    }
 //   cout << line << endl;
  }

  streampos start = input.tellg();
  input.seekg(0, input.end);
  int length = input.tellg() - start;
  input.seekg(start);
  int size = length / rowSize;
  cout << "number of points: " << size << endl;
  vector<unsigned char> vertices(length);
  // read data as a block:
  input.read((char *)&vertices[0], length);
  for (int i = 0; i<size; i++)
  {
    Vector3f b = (Vector3f &)vertices[rowSize*i];
    points.push_back(Vector3d(b[0],b[1],b[2]));
    Vector3f n = (Vector3f &)vertices[rowSize*i + normalOffset];
    if (n[0] == n[0] && n[1] == n[1] && n[2] == n[2])
    {
      normals.push_back(Vector3d(n[0], n[1], n[2]));
    }
    else
    {
      normals.push_back(Vector3d(1, 0, 0));
    }
  }
  return true;
}



