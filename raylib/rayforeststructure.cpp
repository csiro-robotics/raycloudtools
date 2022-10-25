// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforeststructure.h"
// #define OUTPUT_MOMENTS  // used in unit tests
#include <unordered_map>

namespace ray
{
Eigen::Array<double, 9, 1> ForestStructure::getMoments() const
{
  Eigen::Array<double, 9, 1> moments;
  moments[0] = (double)trees.size();
  Eigen::Vector3d sum(0, 0, 0);
  Eigen::Vector3d sum_sqr(0, 0, 0);
  double rad = 0;
  double rad_sqr = 0.0;
  double volume = 0.0;
  for (auto &tree : trees)
  {
    Eigen::Vector3d p = tree.segments()[0].tip;
    sum += p;
    sum_sqr += Eigen::Vector3d(p[0] * p[0], p[1] * p[1], p[2] * p[2]);
    double r = tree.segments()[0].radius;
    rad += r;
    rad_sqr += r * r;
    for (auto &segment : tree.segments())
    {
      if (segment.parent_id != -1)
      {
        double length = (segment.tip - tree.segments()[segment.parent_id].tip).norm();
        double r = segment.radius;
        volume += length * r * r;
      }
    }
  }
  std::hash<std::string> hasher;

  size_t attribute_hash = 0;
  for (auto &attribute: trees[0].treeAttributeNames())
  {
    attribute_hash += hasher(attribute);
  }
  for (auto &attribute: trees[0].attributeNames())
  {
    attribute_hash += hasher(attribute);
  }
  double sum_attributes = 0.0; // a total over all the attributes
  for (auto &tree: trees)
  {
    for (auto &segment: tree.segments())
    {
      for (auto &att: segment.attributes)
      {
        sum_attributes += att;
      }
    }
  }
  sum_attributes /= static_cast<double>(trees.size());
  sum_attributes /= static_cast<double>(trees[0].segments().size());
  // keep the hash small enough to be represented uniquely by a double
  attribute_hash = attribute_hash % 100000; 
  moments[1] = sum.norm();
  moments[2] = sum_sqr.norm();
  moments[3] = rad;
  moments[4] = rad_sqr;
  moments[5] = volume * kPi;
  // we should care about attributes too:
  moments[6] = static_cast<double>(trees[0].treeAttributes().size());
  moments[7] = static_cast<double>(attribute_hash);
  moments[8] = sum_attributes;

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
  do  // skip initial comment lines
  {
    if (!std::getline(ifs, line))
    {
      break;
    }
  } while (line[0] == '#');
  if (!ifs) 
  {
    return false;
  }

  const std::string mandatory_text = "x,y,z,radius";
  size_t found = line.find(mandatory_text);
  if (found == std::string::npos)
  {
    std::cerr << "Error: tree files must declare the mandatory format: " << mandatory_text << std::endl;
    return false;
  }
  std::vector<std::string> attributes[2];
  std::string lines[2];
  if (found > 1)
  {
    lines[0] = line.substr(0, found-2);
    std::cout << "tree format: (" << lines[0] << ")" << std::endl; 
  }
  lines[1] = line.substr(found);
  std::cout << "branch format: (" << lines[1] << ")" << std::endl; 

  // check if there is a second x,y,z,radius on the same line, this will indicate a separate set of branch attributes
  int commas_per_segment = 0;
  int commas_per_tree = 0;

  if (!lines[0].empty())
  {
    // parse the attributes line. This line defines any additional attributes beyond mandatory_text
    commas_per_tree = 1 + (int)std::count(lines[0].begin(), lines[0].end(), ',');
//    std::cout << "commas per tree: " << commas_per_tree << std::endl;
    std::istringstream ss(lines[0]);
    // parse each attribute
    for (int i = 0; i < commas_per_tree; i++)
    {
      std::string token;
      std::getline(ss, token, ',');
      attributes[0].push_back(token);
    }
  }
  bool has_parent_id = false;
  {
    // parse the attributes line. This line defines any additional attributes beyond mandatory_text
    commas_per_segment = 1 + (int)std::count(lines[1].begin(), lines[1].end(), ',');
    lines[1] = lines[1].substr(mandatory_text.length());
    if (lines[1][0] == ',')
    {
      lines[1] = lines[1].substr(1);
    }
//    std::cout << "commas per segment: " << commas_per_segment << std::endl;
    std::istringstream ss(lines[1]);
    // parse each attribute
    for (int i = 4; i < commas_per_segment; i++)
    {
      std::string token;
      std::getline(ss, token, ',');
      if (token == "parent_id")  // parent_id is a special case as it is an int rather than double
      {
        has_parent_id = true;
      }
      else
      {
        attributes[1].push_back(token);
      }
    }
  }
  for (int att = 0; att<2; att++)
  {
    if (attributes[att].size() > 0)
    {
      // let the user know what attributes it found
      std::cout << "reading extra " << (att==0 ? "tree" : "branch") << " attributes: ";
      for (auto &at : attributes[att]) 
      {
        std::cout << at << "; ";
      }
      std::cout << std::endl;
    }
  }


  int line_number = 0;
  // now parse the data, one line at a time
  for (std::string line; std::getline(ifs, line); )
  {
    if (line.length() == 0 || line[0] == '#')  // ignore empty and comment lines
    {
      continue;
    }
    // use commas to separate the line
    int num_commas = 1 + (int)std::count(line.begin(), line.end(), ',');
    TreeStructure tree;
    tree.treeAttributeNames() = attributes[0];
    tree.attributeNames() = attributes[1];
    // make sure the number of commas is what we expect
    if ((num_commas-commas_per_tree) % commas_per_segment != 0)
    {
      std::cerr << "Error: line " << line_number << " contains " << num_commas << " commas, which is not " << commas_per_tree << " plus a multiple of " << commas_per_segment << std::endl;
      return false;
    }

    line_number++;
    std::istringstream ss(line);
    for (int i = 0; i < commas_per_tree; i++)
    {
      std::string token;
      if (!std::getline(ss, token, ','))
      {
        break;
      }
      tree.treeAttributes().push_back(std::stod(token.c_str()));  // whole tree attributes
    }    

    // for each segment...
    for (int c = commas_per_tree; c < num_commas; c += commas_per_segment)
    {
      TreeStructure::Segment segment;
      // for each attribute of this segment...
      for (int i = 0; i < commas_per_segment; i++)
      {
        std::string token;
        if (!std::getline(ss, token, ','))
        {
          break;
        }
        if (i < 3)
        {
          segment.tip[i] = std::stod(token.c_str());  // special case for tip x,y,z
        }
        else if (i == 3)
        {
          segment.radius = std::stod(token.c_str());  // special case for radius
        }
        else if (i == 4 && has_parent_id)
        {
          segment.parent_id = std::stoi(token.c_str());  // special case for parent_id
        }
        else
        {
          segment.attributes.push_back(std::stod(token.c_str()));  // user attributes
        }
      }
      // generate the tree segment from the attributes
      tree.segments().push_back(segment);
    }
    // another error check
    if (has_parent_id && tree.segments().size() == 1)
    {
      std::cerr << "Error: format contains parent id, but trunk only (single segment) detected on line " << line_number
                << std::endl;
      return false;
    }
    trees.push_back(tree);
  }
  if (trees.empty())
  {
    return false;
  }
  if (trees[0].segments().empty())
  {
    return false;
  }

#if defined OUTPUT_MOMENTS  // enabled to provide results for the unit tests
  Eigen::Array<double, 9, 1> mom = getMoments();
  std::cout << "load stats: " << std::endl;
  for (int i = 0; i < mom.rows(); i++) 
  {
    std::cout << ", " << mom[i];
  }
  std::cout << std::endl;
#endif  // defined OUTPUT_MOMENTS
  return true;
}

bool ForestStructure::save(const std::string &filename)
{
  if (trees.empty())
  {
    std::cerr << "No data to save to " << filename << std::endl;
    return false;
  }

#if defined OUTPUT_MOMENTS  // enabled to provide results for the unit tests
  Eigen::Array<double, 9, 1> mom = getMoments();
  std::cout << "save stats: " << std::endl;
  for (int i = 0; i < mom.rows(); i++) 
  {
    std::cout << ", " << mom[i];
  }
  std::cout << std::endl;
#endif  // defined OUTPUT_MOMENTS

  std::cout << "outputting tree file: " << filename << std::endl;
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << "# Tree file. Optional per-tree attributes (e.g. 'height,crown_radius, ') followed by 'x,y,z,radius' and any additional per-segment attributes:" << std::endl;
  for (auto &att : trees[0].treeAttributeNames())
  {
    ofs << att << ",";
  }
  if (!trees[0].treeAttributeNames().empty())
  {
    ofs << " ";
  }
  ofs << "x,y,z,radius";
  if (trees[0].segments().size() > 1)
  {
    ofs << ",parent_id";
  }
  if (trees[0].attributeNames().size() > 0)
  {
    std::cout << "Saving additional branch attributes: ";
  }
  for (auto &att : trees[0].attributeNames())
  {
    ofs << "," << att;
    std::cout << att << ", ";
  }
  std::cout << std::endl;
  ofs << std::endl;
  // for each tree
  for (auto &tree : trees)
  {
    // whole tree attributes:
    for (auto &att: tree.treeAttributes())
    {
      ofs << att << ",";
    }
    if (!tree.treeAttributes().empty())
    {
      ofs << " ";
    }
    // for each segment in the tree
    for (size_t i = 0; i < tree.segments().size(); i++)
    {
      auto &segment = tree.segments()[i];
      if (i > 0)
      {
        ofs << ", ";
      }
      // save the mandatory attributes
      ofs << segment.tip[0] << "," << segment.tip[1] << "," << segment.tip[2] << "," << segment.radius;
      if (tree.segments().size() > 1)
      {  
        ofs << "," << segment.parent_id;    // save the special-case parent_id
      }
      for (auto &att : segment.attributes)  // save the user attributes
      {  
        ofs << "," << att;
      }
    }
    ofs << std::endl;
  }
  return true;
}
}  // namespace ray