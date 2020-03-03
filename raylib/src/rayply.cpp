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

// Save the polygon file to disk
void writePly(const std::string &fileName, const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times)
{
  cout << "saving to " << fileName << endl;
  
  vector<Vector3d> rgb(times.size());
  for (int i= 0; i<(int)rgb.size(); i++)
    rgb[i] = redGreenBlue(fmod(times[i], 10.0)/10.0);
  vector<uint32_t> RGB(rgb.size());
  for (unsigned int i = 0; i<rgb.size(); i++)
  {
    RGB[i] = 0;
    for (int j = 0; j<3; j++)
      RGB[i] += uint32_t(rgb[i][j] * 255.0) << (j*8);
    RGB[i] += 255<<24; // i.e. alpha is 255
  }

  vector<Matrix<float, 8, 1> > vertices(ends.size()); // 4d to give space for colour
  for (unsigned int i = 0; i<ends.size(); i++)
  {
    Vector3d n = starts[i] - ends[i];
    vertices[i] << (float)ends[i][0], (float)ends[i][1], (float)ends[i][2], (float)times[i], (float)n[0], (float)n[1], (float)n[2], (float &)RGB[i];
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
  fprintf(fid, "property float t\n"); 
  fprintf(fid, "property float nx\n"); 
  fprintf(fid, "property float ny\n"); 
  fprintf(fid, "property float nz\n"); 
  fprintf(fid, "property uchar red\n"); 
  fprintf(fid, "property uchar green\n"); 
  fprintf(fid, "property uchar blue\n"); 
  fprintf(fid, "property uchar alpha\n"); 
  fprintf(fid, "element face 0\n"); 
  fprintf(fid, "property list uchar int vertex_indices\n"); 
  fprintf(fid, "end_header\n"); 

  fwrite(&vertices[0], sizeof(Matrix<float, 8, 1>), vertices.size(), fid);

  fclose(fid);  
}

#if 0 // is this a better save PLY function?
void savePLY(const std::string &filename, const State &state, const std::vector<PoseEvent> &trajectory, const LaserData &laser, const PointAttributeOptions& options)
{
  struct Property
  {
    function<void(ofstream&)> writeHeader;
    function<void(ofstream&, const unsigned int&)> writeData;
   
    Property()=default;
    Property(function<void(ofstream&)> header, function<void(ofstream&, const unsigned int&)> data) : writeHeader(header), writeData(data) {};
  };
 
  vector<Property> propertyList;
 
  //Covariance
  if(options.saveCovariance && !laser.times.empty())
  {
    vector<double> covariances;
    for(auto&& covarianceEvent : state.covarianceTrajectory)
      covariances.push_back(covarianceEvent.covariance.norm());
    Lookup<double> covarianceLookup (components(state.covarianceTrajectory, time), covariances);
 
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property float covariance" << endl;
          },
          [&laser, covarianceLookup](ofstream& out, const unsigned int& i)
          {
            float value = covarianceLookup.linear(laser.times[i], false);
            out.write((const char*)&value, sizeof(float));
          }));
  }
 
  //Intensity
  if(options.saveIntensity && !laser.intensities.empty())
  {
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property ushort intensity" << endl;
          },
          [&laser](ofstream& out, const unsigned int& i)
          {
            out.write((const char*)&laser.intensities[i], sizeof(uint16_t));
          }));
  }
 
  //Normal
  if(options.saveNormals && !laser.times.empty() && !trajectory.empty())
  {
    vector<Vector3d> normals = generateNormals(laser.positions, components(Lookup<PoseEvent>(trajectory).linearInterpolate(laser.times), pose.position));
   
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property float nx" << endl;
            out << "property float ny" << endl;
            out << "property float nz" << endl;
          },
          [&laser, normals](ofstream& out, const unsigned int& i)
          {
            Vector3f normal = normals[i].cast<float>();
            out.write((const char*)normal.data(), sizeof(float)*3);
          }));
  }
 
  //Angle of incidence
  if(options.saveIncidenceAngle && !laser.times.empty() && !trajectory.empty())
  {
    vector<Vector3d> normals = generateNormals(laser.positions, components(Lookup<PoseEvent>(trajectory).linearInterpolate(laser.times), pose.position));
    Lookup<PoseEvent> trajLookup (trajectory);
   
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property float incidence_angle" << endl;
          },
          [&laser, normals, trajLookup](ofstream& out, const unsigned int& i)
          {
            float incidenceAngle = (180.0/M_PI) * acos(normals[i].dot((trajLookup.linearNormalized(laser.times[i]).pose.position - laser.positions[i]).normalized()));
            out.write((const char*)&incidenceAngle, sizeof(float));
          }));
  }
 
  //Position
  propertyList.push_back(Property(
        [](ofstream& out){
          out << "property double x" << endl;
          out << "property double y" << endl;
          out << "property double z" << endl;
        },
        [&laser](ofstream& out, const unsigned int& i)
        {
          out.write((const char*)laser.positions[i].data(), sizeof(double)*3);
        }));
 
  //Range
  if(options.saveRange && !laser.times.empty() && !trajectory.empty())
  {
    Lookup<PoseEvent> trajLookup (trajectory);
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property float range" << endl;
          },
          [&laser, trajLookup](ofstream& out, const unsigned int& i)
          {
            float range = (trajLookup.linearNormalized(laser.times[i]).pose.position - laser.positions[i]).norm();
            out.write((const char*)&range, sizeof(float));
          }));
  }
 
  //Return number
  if(options.saveReturnNumber && std::any_of(laser.flags.begin(), laser.flags.end(), [](const LaserData::Flags& flag){return !((bool)(flag & LaserData::Flag::SINGLE_RETURN));}))
  {
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property uchar normalized_return_number" << endl;
          },
          [&laser](ofstream& out, const unsigned int& i)
          {
            uint8_t value;
            if(laser.flags[i] & LaserData::Flag::SINGLE_RETURN)//Single return
              value = 0;
            else if(laser.flags[i] & LaserData::Flag::FIRST_RETURN)//First return
              value = 1;
            else if(laser.flags[i] & LaserData::Flag::LAST_RETURN)//Last return
              value = numeric_limits<uint8_t>::max();
            else//Mid return TODO Is the order required?
              value = numeric_limits<uint8_t>::max()/2;
           
            out.write((const char*)&value, sizeof(uint8_t));
          }));
  }
 
  //Ring number
  if(options.saveRingNumber && !laser.rings.empty())
  {
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property uchar ring_number" << endl;
          },
          [&laser](ofstream& out, const unsigned int& i)
          {
            out.write((const char*)&laser.rings[i], sizeof(uint8_t));
          }));
  }
 
  //Shape
  if(options.saveShapes && !laser.times.empty() && !trajectory.empty())
  {
    vector<Vector3d> shapes = generateShapes(laser.positions, components(Lookup<PoseEvent>(trajectory).linearInterpolate(laser.times), pose.position));
   
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property uchar planarity" << endl;
            out << "property uchar sphericity" << endl;
            out << "property uchar cylindricality" << endl;
          },
          [&laser, shapes](ofstream& out, const unsigned int& i)
          {
            for(unsigned int j = 0; j < 3; ++j)
            {
              uint8_t shape = clamped(shapes[i][j], 0.0, 1.0) * std::numeric_limits<uint8_t>::max();
              out.write((const char*)&shape, sizeof(uint8_t));
            }
          }));
  }
 
  //Strongest Return
  if(options.saveStrongestReturn && std::any_of(laser.flags.begin(), laser.flags.end(), [](const LaserData::Flags& flag){return (bool)(flag & LaserData::Flag::STRONGEST_RETURN);}))
  {
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property uchar strongest_return" << endl;
          },
          [&laser](ofstream& out, const unsigned int& i)
          {
            uint8_t value = (bool)(laser.flags[i] & LaserData::Flag::STRONGEST_RETURN);
            out.write((const char*)&value, sizeof(uint8_t));
          }));
  }
 
  //Time
  if(options.saveTime && !laser.times.empty())
  {
    propertyList.push_back(Property(
          [](ofstream& out){
            out << "property double time" << endl;
          },
          [&laser](ofstream& out, const unsigned int& i)
          {
            out.write((const char*)&laser.times[i], sizeof(double));
          }));
  }
 
  //PLY output
  {
    const vector<Property>& properties = propertyList;
   
    std::ofstream ply (filename, ios::binary | ios::out);
   
    ply << "ply" << endl;
    ply << "format binary_little_endian 1.0" << endl;
    ply << "comment nrrslam" << endl; //TODO version
    ply << "element vertex " << laser.size() << endl;
    for(auto&& property : properties)
      property.writeHeader(ply);
    ply << "element face 0" << endl;
    ply << "property list uchar int vertex_indices" << endl;
    ply << "end_header" << endl;
   
    //Data
    for(unsigned int i = 0; i < laser.size(); ++i)
    {
      for(auto&& property : properties)
        property.writeData(ply, i);
    }
  }
}
#endif

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


bool readPly(const std::string &fileName, std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times)
{
  cout << "reading from " << fileName << endl;
  ifstream input(fileName.c_str());
  if (!input.is_open())
  {
    cerr << "Couldn't open file: " << fileName << endl;
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
    if (line.find("property uchar") != string::npos)
    {
      rowSize++;
    }
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
    Vector4f b = (Vector4f &)vertices[rowSize*i];
    Vector3d end(b[0],b[1],b[2]);
    ends.push_back(Vector3d(end));
    times.push_back(b[3]);
    Vector3f n = (Vector3f &)vertices[rowSize*i + normalOffset];
    starts.push_back(end + Vector3d(n[0], n[1], n[2]));
  }
  return true; 
}



