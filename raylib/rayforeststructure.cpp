#include "rayforeststructure.h"
namespace ray
{
Eigen::Array<double, 6, 1> ForestStructure::getMoments() const
{
  Eigen::Array<double, 6, 1> moments;
  moments[0] = (double)trees.size();
  Eigen::Vector3d sum(0,0,0);
  Eigen::Vector3d sum_sqr(0,0,0);
  double rad = 0;
  double rad_sqr = 0.0;
  double volume = 0.0;
  for (auto &tree: trees)
  {
    Eigen::Vector3d p = tree.segments()[0].tip;
    sum += p;
    sum_sqr += Eigen::Vector3d(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
    double r = tree.segments()[0].radius;
    rad += r;
    rad_sqr += r*r;
    for (auto &segment: tree.segments())
    {
      if (segment.parent_id != -1)
      {
        double length = (segment.tip - tree.segments()[segment.parent_id].tip).norm();
        double r = segment.radius;
        volume += length * r * r;
      }
    }
  }
  moments[1] = sum.norm();
  moments[2] = sum_sqr.norm();
  moments[3] = rad;
  moments[4] = rad_sqr;
  moments[5] = volume * kPi;
  return moments;
}

// Parse the tree file into a vector of tree structures
bool ForestStructure::load(const std::string &filename)
{
  std::cout << "loading tree file: " << filename << std::endl;
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << std::endl;
    return false;
  }  
  std::string line;
  do // skip initial comment lines
  {
    std::getline(ifs, line);
  } while (line[0] == '#');

  const std::string mandatory_text = "x,y,z,radius";
  if (line.substr(0, mandatory_text.length()) != mandatory_text)
  {
    std::cerr << "Error: tree files must start with the mandatory format: " << mandatory_text << std::endl;
    return false;
  }

  // parse the attributes line. This line defines any additional attributes beyond mandatory_text
  int commas_per_segment = 1 + (int)std::count(line.begin(), line.end(), ',');
  line = line.substr(mandatory_text.length());
  if (line[0] == ',')
    line = line.substr(1);  
  std::istringstream ss(line);
  std::vector<std::string> attributes;
  bool has_parent_id = false;
  // parse each attribute
  for (int i = 4; i<commas_per_segment; i++)
  {
    std::string token;
    std::getline(ss, token, ',');
    if (token == "parent_id") // parent_id is a special case as it is an int rather than double
      has_parent_id = true;
    else
      attributes.push_back(token);
  }
  if (attributes.size() > 0)
  {
    // let the user know what attributes it found
    std::cout << "reading extra tree attributes: ";
    for (auto &at: attributes)
      std::cout << at << "; ";
    std::cout << std::endl;
  }

  int line_number = 0;
  // now parse the data, one line at a time
  while (!ifs.eof())
  {
    std::string line;
    std::getline(ifs, line);
    if (line.length() == 0 || line[0] == '#') // ignore empty and comment lines
      continue;
    // use commas to separate the line
    int num_commas = 1 + (int)std::count(line.begin(), line.end(), ',');
    TreeStructure tree;
    tree.attributes() = attributes;
    // make sure the number of commas is what we expect
    if (num_commas > commas_per_segment && !has_parent_id)
    {
      std::cerr << "Error: trees with multiple segments need parent_id field to connect them" << std::endl;
      return false;
    }
    if (num_commas%commas_per_segment != 0)
    {
      std::cerr << "Error: line " << line_number << " contains " << num_commas << " commas, which is not divisible by " << commas_per_segment << std::endl;
      return false;
    }

    line_number++;
    std::istringstream ss(line);
    // for each segment...
    for (int c = 0; c<num_commas; c += commas_per_segment)
    {
      TreeStructure::Segment segment; 
      // for each attribute of this segment...
      for (int i = 0; i<commas_per_segment; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i<3)
          segment.tip[i] = std::stod(token.c_str()); // special case for tip x,y,z
        else if (i==3)
          segment.radius = std::stod(token.c_str()); // special case for radius
        else if (i==4 && has_parent_id)
          segment.parent_id = std::stoi(token.c_str()); // special case for parent_id
        else
          segment.attributes.push_back(std::stod(token.c_str())); // user attributes
      }
      // generate the tree segment from the attributes
      tree.segments().push_back(segment);
    }
    // another error check
    if (has_parent_id && tree.segments().size() == 1)
    {
      std::cerr << "Error: format contains parent id, but trunk only (single segment) detected on line " << line_number << std::endl;
      return false;
    }
    trees.push_back(tree);
  }
  if (trees.empty())
    return false;
  if (trees[0].segments().empty())
    return false;

#if defined OUTPUT_MOMENTS // enabled to provide results for the unit tests
  Eigen::Array<double, 6, 1> mom = getMoments();
  std::cout << "stats: " << std::endl;
  for (int i = 0; i<mom.rows(); i++)
    std::cout << ", " << mom[i];
  std::cout << std::endl;
#endif // defined OUTPUT_MOMENTS    
  return true;
}

bool ForestStructure::save(const std::string &filename)
{
  if (trees.empty())
  {
    std::cerr << "No data to save to " << filename << std::endl;
    return false;
  }
  std::cout << "outputting tree file: " << filename << std::endl;
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree file:" << std::endl;
  ofs << "x,y,z,radius";
  if (trees[0].segments().size() > 1)
    ofs << ",parent_id";
  if (trees[0].attributes().size() > 0)
    std::cout << "Saving additional attributes: ";
  for (auto &att: trees[0].attributes())
  {
    ofs << "," << att;
    std::cout << att << ", ";
  }
  std::cout << std::endl;
  ofs << std::endl;
  // for each tree
  for (auto &tree: trees)
  {
    // for each segment in the tree
    for (size_t i = 0; i<tree.segments().size(); i++)
    {
      auto &segment = tree.segments()[i];
      // save the mandatory attributes
      ofs << segment.tip[0] << "," << segment.tip[1] << "," << segment.tip[2] << "," << segment.radius;
      if (tree.segments().size() > 1)
        ofs << "," << segment.parent_id; // save the special-case parent_id
      for (auto &att: segment.attributes) // save the user attributes
        ofs << "," << att;
      if (i != tree.segments().size()-1)
        ofs << ", ";
    }
    ofs << std::endl;
  }
  return true;      
}
} // namespace ray