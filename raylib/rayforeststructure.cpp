// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforeststructure.h"
// #define OUTPUT_MOMENTS  // used in unit tests
#include <unordered_map>
#include <complex>

namespace ray
{
typedef std::complex<double> Comp;

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
  double sum_tree_attributes = 0.0; // a total over all the tree attributes
  for (auto &tree: trees)
  {
    for (auto &att: tree.treeAttributes())
    {
      sum_tree_attributes += att;
    }
  }
  sum_tree_attributes /= static_cast<double>(trees.size());
  double sum_attributes = 0.0; // a total over all the attributes
  int num_segments = 0;
  for (auto &tree: trees)
  {
    for (auto &segment: tree.segments())
    {
      num_segments++;
      for (auto &att: segment.attributes)
      {
        sum_attributes += att;
      }
    }
  }
  sum_attributes /= static_cast<double>(trees.size());
  sum_attributes /= static_cast<double>(num_segments);
  // keep the hash small enough to be represented uniquely by a double
  attribute_hash = attribute_hash % 100000; 
  moments[1] = sum.norm();
  moments[2] = sum_sqr.norm();
  moments[3] = rad;
  moments[4] = rad_sqr;
  moments[5] = volume * kPi;
  // we should care about attributes too:
  moments[6] = static_cast<double>(attribute_hash);
  moments[7] = sum_tree_attributes;
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
  comments.clear();
  do  // skip initial comment lines
  {
    if (!std::getline(ifs, line))
    {
      break;
    }
    if (line[0] == '#')
    {
      comments.push_back(line);
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
  }
  lines[1] = line.substr(found);

  // check if there is a second x,y,z,radius on the same line, this will indicate a separate set of branch attributes
  int commas_per_segment = 0;
  int commas_per_tree = 0;

  if (!lines[0].empty())
  {
    // parse the attributes line. This line defines any additional attributes beyond mandatory_text
    commas_per_tree = 1 + (int)std::count(lines[0].begin(), lines[0].end(), ',');
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
  ofs << std::setprecision(4) << std::fixed;
  if (comments.empty())
  {
    ofs << "# Tree file. Optional per-tree attributes (e.g. 'height,crown_radius, ') followed by 'x,y,z,radius' and any additional per-segment attributes:" << std::endl;
  }
  else
  {
    for (auto &line: comments)
    {  
      ofs << line << std::endl;
    }
  }
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

void ForestStructure::splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside)
{
  // first implementation is gonna be slow I guess... 
  // I could either grid up the cylinder indices, or I could grid up the point indices.... I wonder what is better...
  // points are easier in that they have no width... 
  double voxel_width = 0.2; // TODO: where to get grid cell size from
  Eigen::Vector3d minbound = cloud.calcMinBound();
  Grid<int> grid(minbound, cloud.calcMaxBound(), voxel_width); 
  // fill the acceleration structure
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    grid.insert(grid.index(cloud.ends[i]), (int)i);
  }

  // now find all ends that we can remove on a per-segment basis:
  std::vector<bool> remove(cloud.ends.size(), false);
  for (auto &tree: trees)
  {
    for (auto &segment: tree.segments())
    {
      if (segment.parent_id != -1)
      {
        Eigen::Vector3d pos1 = tree.segments()[segment.parent_id].tip;
        Eigen::Vector3d pos2 = segment.tip;
        double r = segment.radius;
        Eigen::Vector3d dir = (pos2 - pos1).normalized();
        pos1 -= dir*offset;
        pos2 += dir*offset;
        r += offset;
        double len = (pos2 - pos1).norm();
        Eigen::Vector3d one(1,1,1);
        Eigen::Vector3i minindex = grid.index(minVector(pos1, pos2) - r*one);
        Eigen::Vector3i maxindex = grid.index(maxVector(pos1, pos2) + r*one);
        for (int i = minindex[0]; i<=maxindex[0]; i++)
        {
          for (int j = minindex[1]; j<=maxindex[1]; j++)
          {
            for (int k = minindex[2]; k<=maxindex[2]; k++)
            {
              auto &cell = grid.cell(i,j,k);
              for (auto id: cell.data)
              {
                // intersect the point cloud.ends[id] with segment:
                const Eigen::Vector3d &p = cloud.ends[id];
                double d = (p - pos1).dot(dir);
                if (d < 0 || d > len)
                {
                  continue;
                }
                Eigen::Vector3d closest = pos1 + dir*d;
                double r2 = (p - closest).squaredNorm();
                if (r2 < r*r) // inside the cylinder
                {
                  remove[id] = true;
                }
              }
            }
          }
        }
      }
    }
  }

  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    Cloud &dest = remove[i] ? inside : outside;
    dest.addRay(cloud.starts[i], cloud.ends[i], cloud.times[i], cloud.colours[i]);
  }
}

// add a single section of a capsule. Each one is like a node in the polyline with a radius.
void addCapsulePiece(Mesh &mesh, int wind, const Eigen::Vector3d &pos, const Eigen::Vector3d &side1,
                     const Eigen::Vector3d &side2, double radius, const RGBA &rgba, bool cap_start, bool cap_end, bool add_uvs)
{
  // this is the square root of the volume of a cylinder divided by the volume of the 14 sided polyhderon
  // representing the cylinder (a kind of twisted hexagonal prism). This constant only works for this hexagonal polyhedron
  const double radius_scale = 1.07234; // to keep the volume about equal to that of the cylinder
  const int start_index = static_cast<int>(mesh.vertices().size());
  const Eigen::Vector3i start_indices(start_index, start_index, start_index);  // start indices
  std::vector<Eigen::Vector3i> &indices = mesh.indexList();
  std::vector<Eigen::Vector3cf> &uvs = mesh.uvList();
  std::vector<Eigen::Vector3d> vertices;
  vertices.reserve(7);
  Eigen::Vector3d dir = side2.cross(side1);
  if (cap_start)
    vertices.push_back(pos - radius_scale * radius * dir);

  // add the six vertices in the circumferential ring for this point along the branch
  Eigen::Vector3cf uv;
  for (int i = 0; i < 6; i++)
  {
    const double pi = 3.14156;
    double angle = (static_cast<double>(i) * 2.0 + static_cast<double>(wind)) * pi / 6.0;

    vertices.push_back(pos + radius_scale * radius * (side1 * std::sin(angle) + side2 * std::cos(angle)));
    // the indexing is a bit more complicated, to connect the vertices with triangles
    if (cap_start)
    {
      indices.push_back(start_indices + Eigen::Vector3i(0, 1 + i, 1 + ((i + 1) % 6)));
      if (add_uvs)
      {
        uv[0] = uv[1] = uv[2] = Comp(0.0f,0.0f);
        uvs.push_back(uv);
      }
    }
    else
    {
      if (add_uvs)
      {
        double I = ((double)i + 0.5*(double)wind) / 6.0;
        uv[2] = Comp(I + 0.0/6.0, 0.0);
        uv[1] = Comp(I + 0.5/6.0, 1.0);
        uv[0] = Comp(I + 1.0/6.0, 0.0);
        uvs.push_back(uv);
        uv[2] = Comp(I + 1.5/6.0, 1.0);
        uv[1] = Comp(I + 1.0/6.0, 0.0);
        uv[0] = Comp(I + 0.5/6.0, 1.0);
        uvs.push_back(uv);
      }
      indices.push_back(start_indices + Eigen::Vector3i(i - 6, i, ((i + 1) % 6) - 6));
      indices.push_back(start_indices + Eigen::Vector3i((i + 1) % 6, ((i + 1) % 6) - 6, i));
    }
  }
  if (cap_end)
  {
    vertices.push_back(pos + radius_scale * radius * dir);
    for (int i = 0; i < 6; i++)
    {
      if (add_uvs)
      {
        uv[0] = uv[1] = uv[2] = Comp(0.0f,1.0f);
        uvs.push_back(uv);
      }
      indices.push_back(start_indices + Eigen::Vector3i(6, (i + 1) % 6, i));
    }
  }
  // add these vertices into the mesh
  mesh.vertices().insert(mesh.vertices().end(), vertices.begin(), vertices.end());
  for (size_t i = 0; i < vertices.size(); i++)
  {
    mesh.colours().push_back(rgba);
  }
}

/// @brief This converts the piecewise cylindrical model into a smoother mesh than individual capsule meshes
///        Specifically, each branch (from its base up through the widest radius at each bifurcation) is a continuous
///        mesh with 6 vertices around its circumference. This is equivalent to the capsules being connected
///        wherever it is a continuation of the branch. The result is fewer triangles and a smoother result.
/// @param mesh the mesh object to generate into
/// @param red_id the first colour channel id, used to colour the trees
/// @param red_scale scale on the red colour component
/// @param green_scale scale on the green channel
/// @param blue_scale scale on the blue channel
void ForestStructure::generateSmoothMesh(Mesh &mesh, int red_id, double red_scale,
                        double green_scale, double blue_scale, bool add_uvs)
{
  for (const auto &tree : trees)
  {
    const auto &segments = tree.segments();
    // first generate the list of children for each segment
    std::vector<std::vector<int>> children(segments.size());
    for (size_t i = 0; i < segments.size(); i++)
    {
      const auto &segment = segments[i];
      int parent = segment.parent_id;
      if (parent != -1)
      {
        children[parent].push_back(static_cast<int>(i));
      }
    }
    // now generate the set of root segments
    std::vector<int> roots;
    for (int i = 1; i < static_cast<int>(segments.size()); i++)
    {
      if (segments[i].parent_id > 0)
      {
        break;
      }
      roots.push_back(i);
    }

    RGBA rgba;
    // for each root, we follow up through the largest child to make a contiguous branch
    for (size_t i = 0; i < roots.size(); i++)
    {
      int root_id = roots[i];
      Eigen::Vector3d normal(1, 2, 3);  // unspecial 'up' direction for placing vertices along the circumference

      // we iterate through this list and grow it at the same time
      std::vector<int> childlist = { root_id };
      int wind = 0;  // this is what rotates the vertices half a triangle width at each segment, to keep the triangles isoceles
      for (size_t j = 0; j < childlist.size(); j++)
      {
        int child_id = childlist[j];
        int par_id = segments[child_id].parent_id;
        // generate an orthogonal frame for each ring of vertices to sit on
        Eigen::Vector3d dir = (segments[child_id].tip - segments[par_id].tip).normalized();
        Eigen::Vector3d axis1 = normal.cross(dir).normalized();
        Eigen::Vector3d axis2 = axis1.cross(dir);
        rgba = RGBA::treetrunk();  // standardised colour in raycloudtools
        if (red_id != -1)               // use the per-segment colour if it exists (e.g. from treecolour)
        {
          rgba.red = uint8_t(std::min(red_scale * segments[child_id].attributes[red_id], 255.0));
          rgba.green = uint8_t(std::min(green_scale * segments[child_id].attributes[red_id + 1], 255.0));
          rgba.blue = uint8_t(std::min(blue_scale * segments[child_id].attributes[red_id + 2], 255.0));
        }

        if (child_id == root_id)  // add the base cap of the cylinder if we are at the root of the branch
        {
          addCapsulePiece(mesh, wind, segments[par_id].tip, axis1, axis2, segments[child_id].radius, rgba, true, false, add_uvs);
        }

        wind++;
        std::vector<int> kids = children[child_id];
        if (kids.empty())  // add the end cap of the cylinder if we are at the end of the whole branch
        {
          addCapsulePiece(mesh, wind, segments[child_id].tip, axis1, axis2, segments[child_id].radius, rgba, false, true, add_uvs);
          break;
        }
        // now find the maximum radius subbranch
        double max_rad = 0.0;
        int max_k = 0;
        for (int k = 0; k < static_cast<int>(kids.size()); k++)
        {
          double rad = segments[kids[k]].radius;
          if (rad > max_rad)
          {
            max_rad = rad;
            max_k = k;
          }
        }
        for (int k = 0; k < static_cast<int>(kids.size()); k++)
        {
          if (k != max_k)
          {
            roots.push_back(kids[k]);  // all other subbranches get added to the list, to be iterated over on their turn
          }
        }

        int next_id = kids[max_k];
        Eigen::Vector3d dir2 = (segments[next_id].tip - segments[child_id].tip).normalized();

        Eigen::Vector3d top_dir = (dir2 + dir).normalized();  // here we average the directions of the two segments
        // and generate an orthogonal basis for the ring of points on the branch
        Eigen::Vector3d mid_axis1 = normal.cross(top_dir).normalized();
        Eigen::Vector3d mid_axis2 = mid_axis1.cross(top_dir);
        normal = -mid_axis2;
        // add the ring of points
        addCapsulePiece(mesh, wind, segments[child_id].tip, mid_axis1, mid_axis2, segments[child_id].radius, rgba,
                        false, false, add_uvs);
        // add the biggest subbranch to the list, so we continue to build the branch
        childlist.push_back(kids[max_k]);
      }
    }
  }
}

void ForestStructure::reindex()
{
  for (auto &tree: trees)
  {
    tree.reindex();
  }
}


}  // namespace ray