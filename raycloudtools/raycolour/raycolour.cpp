// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/extraction/raysegment.h"
#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/rayparse.h"
#define STB_IMAGE_IMPLEMENTATION
#include "raylib/imageread.h"

#include <nabo/nabo.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Colour the ray cloud, and/or shade it" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycolour raycloud time          - colour by time (optional on all types)" << std::endl;
  std::cout << "                   height        - colour by height" << std::endl;
  std::cout << "                   shape         - colour by geometry shape (r,g,b: spherical, cylinderical, planar)" << std::endl;
  std::cout << "                   normal        - colour by normal" << std::endl;
  std::cout << "                   alpha         - colour by alpha channel (which typically represents intensity)" << std::endl;
  std::cout << "                   alpha 1       - set only alpha channel (zero represents unbounded rays)" << std::endl;
  std::cout << "                   1,1,1         - set r,g,b" << std::endl;
  std::cout << "                   branches      - red and green are lidar intensity and cylindricality respectively, greater for branches than for leaves" << std::endl;
  std::cout << "                   image planview.png - colour all points from image, stretched to fit the point bounds" << std::endl;
  std::cout << "                         --lit   - shaded (slow on large datasets)" << std::endl;
  // clang-format on
  exit(exit_code);
}

/// This is a 2D area measurement on a 3x3 matrix, it is the measure halfway between the trace (1D) and the determinant
/// (3D) of the matrix.
double areaMeasure(const Eigen::Matrix3d &mat)
{
  return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0) + mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0) +
         mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1);
}

// shortcut, to place the red green blue spectrum into the RGBA structure
void spectrumRGB(double value, ray::RGBA &colour)
{
  const Eigen::Vector3d col = ray::redGreenBlueSpectrum(value);
  colour.red = static_cast<uint8_t>(255.0 * col[0]);
  colour.green = static_cast<uint8_t>(255.0 * col[1]);
  colour.blue = static_cast<uint8_t>(255.0 * col[2]);
}

/// Function to colour the cloud from a horizontal projection of a supplied image, stretching to match the cloud bounds.
void colourFromImage(const std::string &cloud_file, const std::string &image_file, ray::CloudWriter &writer)
{
  ray::Cloud::Info info;
  if (!ray::Cloud::getInfo(cloud_file, info))
  {
    usage();
  }
  const ray::Cuboid bounds = info.ends_bound;
  stbi_set_flip_vertically_on_load(1);
  int width, height, num_channels;
  unsigned char *image_data = stbi_load(image_file.c_str(), &width, &height, &num_channels, 0);
  const double width_x = (bounds.max_bound_[0] - bounds.min_bound_[0]) / (double)width;
  const double width_y = (bounds.max_bound_[1] - bounds.min_bound_[1]) / (double)height;
  if (std::max(width_x, width_y) > 1.05 * std::min(width_x, width_y))
  {
    std::cout << "Warning: image aspect ratio does not match aspect match of point cloud ("
              << bounds.max_bound_[0] - bounds.min_bound_[0] << " x " << bounds.max_bound_[1] - bounds.min_bound_[1]
              << ", stretching to fit) " << std::endl;
  }

  auto colour_from_image = [&bounds, &writer, &image_data, width_x, width_y, width, height, num_channels](
                             std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                             std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      const int ind0 = static_cast<int>((ends[i][0] - bounds.min_bound_[0]) / width_x);
      const int ind1 = static_cast<int>((ends[i][1] - bounds.min_bound_[1]) / width_y);
      const int index = num_channels * (ind0 + width * ind1);
      colours[i].red = image_data[index];
      colours[i].green = image_data[index + 1];
      colours[i].blue = image_data[index + 2];
    }
    writer.writeChunk(starts, ends, times, colours);
  };

  if (!ray::Cloud::read(cloud_file, colour_from_image))
    usage();

  stbi_image_free(image_data);
}

// Colours the ray cloud based on the specified arguments
int rayColour(int argc, char *argv[])
{
  ray::FileArgument cloud_file, image_file;
  ray::KeyChoice colour_type({ "time", "height", "shape", "normal", "alpha", "branches" });
  ray::OptionalFlagArgument lit("lit", 'l');
  ray::Vector3dArgument col(0.0, 1.0);
  ray::DoubleArgument alpha(0.0, 1.0);
  ray::TextArgument alpha_text("alpha"), image_text("image");
  const bool standard_format = ray::parseCommandLine(argc, argv, { &cloud_file, &colour_type }, { &lit });
  const bool flat_colour = ray::parseCommandLine(argc, argv, { &cloud_file, &col }, { &lit });
  const bool flat_alpha = ray::parseCommandLine(argc, argv, { &cloud_file, &alpha_text, &alpha }, { &lit });
  const bool image_format = ray::parseCommandLine(argc, argv, { &cloud_file, &image_text, &image_file }, { &lit });
  if (!standard_format && !flat_colour && !flat_alpha && !image_format)
    usage();

  std::string in_file = cloud_file.name();
  const std::string out_file = cloud_file.nameStub() + "_coloured.ply";
  const std::string type = colour_type.selectedKey();
  uint8_t split_alpha = 100;

  if (type != "shape" && type != "normal" && type != "branches")  // chunk loading possible for simple cases
  {
    ray::CloudWriter writer;
    if (!writer.begin(out_file))
      usage();

    auto colour_rays = [flat_colour, flat_alpha, &type, &col, &alpha, &writer, &split_alpha](
                         std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<ray::RGBA> &colours) {
      if (flat_colour)
      {
        for (auto &colour : colours)
        {
          colour.red = (uint8_t)(255.0 * col.value()[0]);
          colour.green = (uint8_t)(255.0 * col.value()[1]);
          colour.blue = (uint8_t)(255.0 * col.value()[2]);
        }
      }
      else if (flat_alpha)
      {
        for (auto &colour : colours)
        {
          colour.alpha = (uint8_t)(255.0 * alpha.value());
        }
      }
      else  // standard_format
      {
        if (type == "time")
        {
          const double colour_repeat_period = 60.0;  // repeating per minute gives a quick way to assess the scan length
          for (size_t i = 0; i < ends.size(); i++)
          {
            spectrumRGB(times[i] / colour_repeat_period, colours[i]);
          }
        }
        else if (type == "height")
        {
          const double wavelength = 10.0;
          for (size_t i = 0; i < ends.size(); i++)
          {
            spectrumRGB(ends[i][2] / wavelength, colours[i]);
          }
        }
        else if (type == "alpha")
        {
          for (auto &colour : colours)
          {
            const Eigen::Vector3d col_vec = ray::redGreenBlueGradient(colour.alpha / 255.0);
            colour.red = uint8_t(255.0 * col_vec[0]);
            colour.green = uint8_t(255.0 * col_vec[1]);
            colour.blue = uint8_t(255.0 * col_vec[2]);
          }
        }
        else
          usage();
      }
      writer.writeChunk(starts, ends, times, colours);
    };

    if (image_format)
    {
      colourFromImage(cloud_file.name(), image_file.name(), writer);
    }
    else if (!ray::Cloud::read(cloud_file.name(), colour_rays))
    {
      usage();
    }
    writer.end();
    if (!lit.isSet())
      return 0;
    in_file = out_file;  // when lit we have to load again, from the saved output file
    std::cout << "reopening file for lighting..." << std::endl;
  }

  // The remainder cannot currently be done with chunk loading
  ray::Cloud cloud;
  if (!cloud.load(in_file))
    usage();

  // what I need is the normal, curvature, eigenvalues, per point.
  struct Data
  {
    Eigen::Vector3d normal;
    Eigen::Vector3d dimensions;
    double curvature;
  };
  std::vector<Data> data(cloud.ends.size());
  const int search_size = std::min(20, (int)cloud.ends.size() - 1);
  std::vector<Eigen::Vector3d> centroids;
  std::vector<Eigen::Vector3d> dimensions;
  std::vector<Eigen::Vector3d> normals;
  Eigen::MatrixXi indices;
  std::vector<Eigen::Vector3d> *cents = nullptr, *dims = nullptr, *norms = nullptr;
  std::vector<Eigen::Matrix3d> *mats = nullptr;
  Eigen::MatrixXi *inds = nullptr;
  double max_distance = 0.0;

  // what do we want to calculate...
  bool calc_surfels = true;
  if (type == "normal")
    norms = &normals;
  else if (type == "shape")
    dims = &dimensions;
  else if (type == "branches")
  {
    inds = &indices;
    dims = &dimensions;
  }
  else
    calc_surfels = lit.isSet();
  if (lit.isSet())
  {
    norms = &normals;
    inds = &indices;
    cents = &centroids;
  }

  if (calc_surfels)
    cloud.getSurfels(search_size, cents, norms, dims, mats, inds, max_distance, false);
  if (type == "shape")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      const double sphericity = dimensions[i][0] / dimensions[i][2];
      const double cylindricality = 1.0 - dimensions[i][1] / dimensions[i][2];
      const double planarity = 1.0 - dimensions[i][0] / dimensions[i][1];
      cloud.colours[i].red = (uint8_t)(255.0 * sphericity);
      cloud.colours[i].green = (uint8_t)(255.0 * cylindricality);
      cloud.colours[i].blue = (uint8_t)(255.0 * planarity);
    }
  }
  else if (type == "normal")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      cloud.colours[i].red = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][0]));
      cloud.colours[i].green = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][1]));
      cloud.colours[i].blue = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][2]));
    }
  }
  // colour in order to distinguish branches.
  // The red channel is a function of the lidar return intensiity, which is typically higher on
  // branches than on leaves.
  // The green channel is a measure of the cylindricality of the neighbourhood of points
  // The resulting colour can be used to segment out branches by thresholding using
  //  raysplit cloud.ply colour x,y,0 for a choice of x, y
  else if (type == "branches")
  {
    std::vector<uint8_t> cols;
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      // 1. red is median alpha value, rescaled
      // we use the median of the neighbour points to be robust to noise
      cols.clear();
      cols.push_back(cloud.colours[i].alpha);
      for (int j = 0; j < 4 && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
      {
        cols.push_back(cloud.colours[indices(j, i)].alpha);
      }
      if (cols.size() == 1)
        cloud.colours[i].red = cloud.colours[i].alpha;
      else
        cloud.colours[i].red = (uint8_t)ray::median(cols);

      // we also compensate for a change in intensity with scale
      double range = (cloud.ends[i] - cloud.starts[i]).norm();
      double half_range = 100.0;
      double red = (double)cloud.colours[i].red / (1.0 + range / half_range);
      double scale = 2.0;
      cloud.colours[i].red =
        (uint8_t)std::max(0, std::min(127 + ((int)(0.5 + red * scale) - (int)(split_alpha * scale)), 255));

      // 2. green is cylindricality
      Eigen::Vector3d mean = cloud.ends[i];  // centroid
      int num = 1;
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
      {
        mean += cloud.ends[indices(j, i)];
        num++;
      }
      mean /= (double)num;
      // get teh scatter matrix of the neighbourhood of points
      Eigen::Matrix3d scatter = (cloud.ends[i] - mean) * (cloud.ends[i] - mean).transpose();
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
      {
        Eigen::Vector3d v = cloud.ends[indices(j, i)] - mean;
        scatter += v * v.transpose();
      }
      scatter /= (double)num;
      // if you divide sqrt(area) by the trace, you get a dimensionless value that
      // is large for disks and spheres, but small for lines/cylinders. So 1 minus this
      // value gives a measure of cylindricality
      double cylind = 1.0 - 3.0 * std::sqrt(areaMeasure(scatter) / 3.0) / scatter.trace();
      // the above measure is smoother in parameter space than previous methods,
      // and avoids the need to do an Eigendecomposition.
      cloud.colours[i].green = (uint8_t)std::max(0.0, 255.0 * std::min(cylind, 1.0));

      // 3. blue is nothing
      cloud.colours[i].blue = 0;
    }
  }
  if (lit.isSet())
  {
    std::vector<double> curvatures(cloud.ends.size());
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0, sum_yy = 0, n = 0;
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
      {
        int id = indices(j, i);
        Eigen::Vector3d flat = cloud.ends[id] - centroids[i];
        double y = flat.dot(normals[i]);
        flat -= y * normals[i];
        double x = flat.squaredNorm();
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        n++;
      }
      const double den = n * sum_xx - sum_x * sum_x;
      if (abs(den) < 1e-8)
        curvatures[i] = 0.0;
      else
        curvatures[i] = (n * sum_xy - sum_x * sum_y) / den;
    }
    Eigen::Vector3d light_dir = Eigen::Vector3d(0.2, 0.4, 1.0).normalized();
    double curve_scale = 4.0;
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      const double scale1 = 0.5 + 0.5 * normals[i].dot(light_dir);
      const double scale2 = 0.5 - 0.5 * curvatures[i] / curve_scale;
      const double s = 0.25 + 0.75 * ray::clamped((scale1 + scale2) / 2.0, 0.0, 1.0);
      cloud.colours[i].red = (uint8_t)((double)cloud.colours[i].red * s);
      cloud.colours[i].green = (uint8_t)((double)cloud.colours[i].green * s);
      cloud.colours[i].blue = (uint8_t)((double)cloud.colours[i].blue * s);
    }
  }
  cloud.save(out_file);

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayColour, argc, argv);
}