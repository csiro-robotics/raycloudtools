// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout << "Colour the ray cloud, and/or shade it" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycolour raycloud time          - colour by time (optional on all types)" << std::endl;
  std::cout << "                   height        - colour by height" << std::endl;
  std::cout << "                   shape         - colour by geometry shape (r,g,b: spherical, cylinderical, planar)" << std::endl;
  std::cout << "                   normal        - colour by normal" << std::endl;
  std::cout << "                   alpha         - colour by alpha channel (which typically represents intensity)" << std::endl;
  std::cout << "                   alpha 1       - set only alpha channel (zero represents unbounded rays)" << std::endl;
  std::cout << "                   1,1,1         - set r,g,b" << std::endl;
  std::cout << "                         --lit   - shaded (slow on large datasets)" << std::endl;
  std::cout << "                   foliage       - uses lidar intensity to split around best split point. r,g,b are different blur levels" << std::endl;
  exit(exit_code);
}

// shortcut, to place the red green blue spectrum into the RGBA structure
void spectrumRGB(double value, ray::RGBA &colour)
{
  const Eigen::Vector3d col = ray::redGreenBlueSpectrum(value);
  colour.red = static_cast<uint8_t>(255.0 * col[0]);
  colour.green = static_cast<uint8_t>(255.0 * col[1]);
  colour.blue = static_cast<uint8_t>(255.0 * col[2]);  
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::KeyChoice colour_type({"time", "height", "shape", "normal", "alpha", "foliage"});
  ray::OptionalFlagArgument lit("lit", 'l');
  ray::Vector3dArgument col(0.0, 1.0);
  ray::DoubleArgument alpha(0.0, 1.0);
  ray::TextArgument alpha_text("alpha");
  const bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &colour_type}, {&lit});
  const bool flat_colour = ray::parseCommandLine(argc, argv, {&cloud_file, &col}, {&lit});
  const bool flat_alpha = ray::parseCommandLine(argc, argv, {&cloud_file, &alpha_text, &alpha}, {&lit});
  if (!standard_format && !flat_colour && !flat_alpha)
    usage();
  
  const std::string type = colour_type.selectedKey();
  std::string in_file = cloud_file.name();
  const std::string out_file = cloud_file.nameStub() + "_coloured.ply";
  uint8_t split_alpha = 90; // was 185

  if (type != "shape" && type != "normal" && type != "foliage") // chunk loading possible for simple cases
  {
    ray::CloudWriter writer;
    if (!writer.begin(out_file))
      usage();

    if (type == "foliage")
    {
      int ticks_per_bucket = 10;
      int num_buckets = 1 + 255 / ticks_per_bucket;
      std::vector<int> histogram(num_buckets, 0);

      auto fill_histogram = [&histogram, &ticks_per_bucket](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &, std::vector<double> &, std::vector<ray::RGBA> &colours)
      {
        for (auto &colour: colours)
        {
          int bucket = (int)colour.alpha / ticks_per_bucket;
          histogram[bucket]++;
        }
      };
      if (!ray::Cloud::read(cloud_file.name(), fill_histogram))
        usage();    

      for (auto &hist: histogram)
        std::cout << hist << ", ";
      std::cout << std::endl;

      int starti = split_alpha / ticks_per_bucket;
      int i;
      for (i = starti; i < num_buckets-1 && histogram[i+1]<histogram[i]; i++);
      int j;
      for (j = starti; j > 0 && histogram[j-1]<histogram[j]; j--);
      if (histogram[j] < histogram[i])
        i = j;
      
   //   split_alpha = (uint8_t)(ticks_per_bucket*i + ticks_per_bucket/2);
      std::cout << "foliage / branch alpha split value: " << (int)split_alpha << std::endl;
    }

    auto colour_rays = [flat_colour, flat_alpha, &type, &col, &alpha, &writer, &split_alpha]
      (std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, 
       std::vector<double> &times, std::vector<ray::RGBA> &colours)
    {
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
      else // standard_format
      {
        if (type == "time")
        {
          const double colour_repeat_period = 60.0; // repeating per minute gives a quick way to assess the scan length
          for (size_t i = 0; i<ends.size(); i++)
          {
            spectrumRGB(times[i]/colour_repeat_period, colours[i]);
          }
        }
        else if (type == "height")
        {
          const double wavelength = 10.0;
          for (size_t i = 0; i<ends.size(); i++)
          {
            spectrumRGB(ends[i][2]/wavelength, colours[i]);
          }
        }
        else if (type == "alpha")
        {
          for (auto &colour: colours)
          {
            const Eigen::Vector3d col_vec = ray::redGreenBlueGradient(colour.alpha/255.0);
            colour.red = uint8_t(255.0 * col_vec[0]);
            colour.green = uint8_t(255.0 * col_vec[1]);
            colour.blue = uint8_t(255.0 * col_vec[2]);
          }
        }
        else if (type == "foliage")
        {
          for (auto &colour: colours)
          {
            int scale = 4;
            int shade = 127 + ((int)colour.alpha - (int)split_alpha) * scale;
            colour.red = colour.green = colour.blue = (uint8_t)std::max(0, std::min(shade, 255));
 //           colour.alpha = 255;
          }
        }
        else
          usage();
      }
      writer.writeChunk(starts, ends, times, colours);
    };

    if (!ray::Cloud::read(cloud_file.name(), colour_rays))
      usage();
    writer.end();
    if (!lit.isSet())
      return 0;
    in_file = out_file; // when lit we have to load again, from the saved output file
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
  int search_size = 20;
  if (type == "foliage")
    search_size = 5;
  std::vector<Eigen::Vector3d> centroids;
  std::vector<Eigen::Vector3d> dimensions;
  std::vector<Eigen::Vector3d> normals;
  Eigen::MatrixXi indices;
  std::vector<Eigen::Vector3d> *cents = NULL, *dims = NULL, *norms = NULL;
  std::vector<Eigen::Matrix3d> *mats = NULL;
  Eigen::MatrixXi *inds = NULL;
  double max_distance = 0.0;

  // what do we want to calculate...
  bool calc_surfels = true;
  if (type == "normal")
    norms = &normals;
  else if (type == "shape")
    dims = &dimensions;
  else if (type == "foliage")
  {
    inds = &indices;
    max_distance = 0.2;
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
    cloud.getSurfels(search_size, cents, norms, dims, mats, inds, max_distance);
  if (type == "shape")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      const double sphericity = dimensions[i][0] / dimensions[i][2];
      const double cylindricality = 1.0 - dimensions[i][1] / dimensions[i][2];
      const double planarity = 1.0 - dimensions[i][0] / dimensions[i][1];
      cloud.colours[i].red = (uint8_t)(255.0 * (0.3 + 0.7 * sphericity));
      cloud.colours[i].green = (uint8_t)(255.0 * (0.3 + 0.7 * cylindricality));
      cloud.colours[i].blue = (uint8_t)(255.0 * (0.3 + 0.7 * planarity));
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
  else if (type == "foliage")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      double alpha = cloud.colours[i].alpha;
      double num = 1;
      for (int j = 0; j<search_size && indices(j, i) > -1; j++)
      {
        alpha += (double)cloud.colours[indices(j, i)].alpha;
        num++;
      }
      cloud.colours[i].red = (uint8_t)(alpha / num);
    }    
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      double alpha = cloud.colours[i].red;
      double num = 1;
      for (int j = 0; j<search_size && indices(j, i) > -1; j++)
      {
        alpha += (double)cloud.colours[indices(j, i)].red;
        num++;
      }
      cloud.colours[i].green = (uint8_t)(alpha / num);
    }   
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      double alpha = cloud.colours[i].green;
      double num = 1;
      for (int j = 0; j<search_size && indices(j, i) > -1; j++)
      {
        alpha += (double)cloud.colours[indices(j, i)].green;
        num++;
      }
      cloud.colours[i].blue = (uint8_t)(alpha / num);
    }  
    for (auto &colour: cloud.colours)
    {
      int scale = 3;
      colour.red = (uint8_t)std::max(0, std::min(127 + ((int)colour.red - split_alpha)*scale, 255));
      colour.green = (uint8_t)std::max(0, std::min(127 + ((int)colour.green - split_alpha)*scale, 255));
      colour.blue = (uint8_t)std::max(0, std::min(127 + ((int)colour.blue - split_alpha)*scale, 255));
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
      for (int j = 0; j < search_size && indices(j, i) > -1; j++)
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
