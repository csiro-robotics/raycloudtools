// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/rayply.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imagewrite.h"

#define DENSITY_MIN_RAYS 10 // larger is more accurate but more blurred. 0 for no adaptive blending

void usage(int exit_code = 1)
{
  std::cout << "Render a ray cloud as an image, from a specified viewpoint" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrender raycloudfile.ply top ends        - render from the top (plan view) the end points" << std::endl;
  std::cout << "                           left            - facing negative x axis" << std::endl;
  std::cout << "                           right           - facing positive x axis" << std::endl;
  std::cout << "                           front           - facing negative y axis" << std::endl;
  std::cout << "                           back            - facing positive y axis" << std::endl;
  std::cout << "                               mean        - mean colour on axis" << std::endl;
  std::cout << "                               sum         - sum colours (globally scaled to colour range)"
    << std::endl;
  std::cout << "                               starts      - render the ray start points" << std::endl;
  std::cout << "                               rays        - render the full set of rays" << std::endl;
  std::cout << "                               density     - shade according to estimated density within pixel"
    << std::endl;
  std::cout << "                               density_rgb - r->g->b colour by estimated density"
    << std::endl;
  std::cout << "                     --pixel_width 0.1     - optional pixel width in m" << std::endl;
  std::cout << "                     --output name.png     - optional output file name. " << std::endl;
  std::cout << "                                             Supports .png, .tga, .hdr, .jpg, .bmp" << std::endl;
  std::cout << "Default output is raycloudfile.png" << std::endl;
  exit(exit_code);
}

struct VoxelGrid
{
  static const int minVoxelHits = 2;
  static constexpr double sphericalDistributionScale = 2.0; // average area scale due to a spherical uniform distribution of leave angles relative to the rays
  
  void calculateDensities(const std::string &file_name);
  Eigen::Vector3d getCentre(const Eigen::Vector3d &pos);
  Eigen::Vector3d getCentre(const Eigen::Vector3i &inds);
  int getIndex(const Eigen::Vector3d &pos);
  int getIndex(const Eigen::Vector3i &inds);

  ray::Cuboid bounds;
  
  struct Voxel
  {
    Voxel(){ numHits = numRays = pathLength = 0.0; }
    float numHits;
    float numRays;
    float pathLength;
    inline double density() const 
    { 
      if (numRays <= minVoxelHits)
        return 0.0;
      return sphericalDistributionScale * (numRays-1.0) * numHits / (1e-10 + numRays*pathLength); 
    } 
    void operator +=(const Voxel &other)
    {
      numHits += other.numHits;
      numRays += other.numRays; 
      pathLength += other.pathLength; 
    }
    Voxel operator *(float scale)
    {
      Voxel voxel;
      voxel.numHits = numHits * scale;
      voxel.numRays = numRays * scale;
      voxel.pathLength = pathLength * scale;
      return voxel;
    }
  };

  std::vector<Voxel> voxels;
  double voxelWidth;
  Eigen::Vector3i voxelDims;
};

// density is the probability of hitting something per metre depth
void VoxelGrid::calculateDensities(const std::string &file_name)
{
  auto calculate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &, std::vector<ray::RGBA> &colours)
  {
    for (size_t i = 0; i<ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end   = ends[i];
      bounds.clipRay(start, end);

      // now walk the voxels
      Eigen::Vector3d dir = end - start;
      Eigen::Vector3d source = (start - bounds.min_bound_)/voxelWidth;
      Eigen::Vector3d p = source;
      Eigen::Vector3d target = (end - bounds.min_bound_)/voxelWidth;
      double length = dir.norm();
      double maxDist = (target - source).norm();
      Eigen::Vector3i inds = p.cast<int>();
      double depth = 0;
      do
      {
        int axis = 0;
        double minL = 1e10;
        for (int k = 0; k<3; k++)
        {
          double l = (dir[k] > 0 ? ceil(p[k]) - p[k] : p[k] - floor(p[k])) * length / abs(dir[k]);
          if (l < minL)
          {
            minL = l;
            axis = k;
          }
        }
        depth += minL + 1e-9;
        inds[axis] += dir[axis] > 0 ? 1 : -1;
        if (inds[axis] < 0 || inds[axis] >= voxelDims[axis])
          break;
        p = source + depth * dir / length;
        int j = getIndex(inds);
        if (colours[i].alpha > 0 && depth > maxDist)
        {
          double d = minL + maxDist - depth;
          voxels[j].pathLength += static_cast<float>(d*voxelWidth);
          voxels[j].numHits++;
          voxels[j].numRays++;
        }
        else
        {
          voxels[j].pathLength += static_cast<float>(minL*voxelWidth); 
          voxels[j].numRays++;
        }
      } while (depth <= maxDist);
    }
  };
  ray::readPly(file_name, true, calculate, 0);
}

Eigen::Vector3d VoxelGrid::getCentre(const Eigen::Vector3d &pos)
{
  Eigen::Vector3d source = (pos - bounds.min_bound_)/voxelWidth;
  Eigen::Vector3i inds = source.cast<int>();
  return getCentre(inds);
}

Eigen::Vector3d VoxelGrid::getCentre(const Eigen::Vector3i &inds)
{
  Eigen::Vector3d res(inds[0], inds[1], inds[2]);
  res += Eigen::Vector3d(0.5,0.5,0.5);
  res = res*voxelWidth + bounds.min_bound_;
  return res;  
}

int VoxelGrid::getIndex(const Eigen::Vector3d &pos)
{
  Eigen::Vector3d source = (pos - bounds.min_bound_)/voxelWidth;
  Eigen::Vector3i inds = source.cast<int>();
  return getIndex(inds);
}

int VoxelGrid::getIndex(const Eigen::Vector3i &inds)
{
  if (inds[0]>=0 && inds[0]<voxelDims[0]
   && inds[1]>=0 && inds[1]<voxelDims[1]
   && inds[2]>=0 && inds[2]<voxelDims[2])
    return inds[0] + inds[1]*voxelDims[0] + inds[2] * voxelDims[0]*voxelDims[1];
  return 0; // error value
}

int main(int argc, char *argv[])
{
  ray::KeyChoice viewpoint({"top", "left", "right", "front", "back"});
  ray::KeyChoice style({"ends", "mean", "sum", "starts", "rays", "density", "density_rgb"});
  ray::DoubleArgument pixel_width(0.0001, 1000.0);
  ray::FileArgument cloud_file, image_file;
  ray::OptionalKeyValueArgument pixel_width_option("pixel_width", 'p', &pixel_width);
  ray::OptionalKeyValueArgument output_file_option("output", 'o', &image_file);
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &viewpoint, &style}, {&pixel_width_option, &output_file_option}))
    usage();
  if (!output_file_option.isSet())
    image_file.name() = cloud_file.nameStub() + ".png";

  ray::Cloud::CloudInfo info;
  if (!ray::Cloud::getInfo(cloud_file.name(), info))
    usage();
  ray::Cuboid bounds = info.ends_bound; // exclude the unbounded ray lengths (e.g. up into the sky)
  double pix_width = pixel_width.value();
  if (!pixel_width_option.isSet())
    pix_width = ray::Cloud::estimatePointSpacing(cloud_file.name(), bounds, info.num_bounded);
  if (pix_width <= 0.0)
    usage();

  Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  int axis = 0;
  if (viewpoint.selectedKey() == "top")
    axis = 2;
  else if (viewpoint.selectedKey() == "front" || viewpoint.selectedKey() == "back")
    axis = 1;
  double dir = 1;
  if (viewpoint.selectedKey() == "left" || viewpoint.selectedKey() == "front")
    dir = -1;
  bool flip_x = viewpoint.selectedKey() == "left" || viewpoint.selectedKey() == "back";
  
  int x_axes[] = {1, 0, 0};
  int y_axes[] = {2, 2, 1};
  int ax1 = x_axes[axis];
  int ax2 = y_axes[axis];
  int width  = 1 + static_cast<int>(extent[ax1] / pix_width);
  int height = 1 + static_cast<int>(extent[ax2] / pix_width);
  int depth  = 1 + static_cast<int>(extent[axis] / pix_width);
  std::cout << "outputting " << width << "x" << height << " image" << std::endl;

  std::vector<Eigen::Vector4d> pixels(width * height); 
  memset(&pixels[0], 0, sizeof(Eigen::Vector4d) * width*height);
  if (style.selectedKey() == "density" || style.selectedKey() == "density_rgb") // special (and complicatd) case
  {
    Eigen::Vector3i dims = (extent/pix_width).cast<int>() + Eigen::Vector3i(1,1,1);
    #if DENSITY_MIN_RAYS > 0
    dims += Eigen::Vector3i(1,1,1); // so that we have extra space to convolve
    #endif
    VoxelGrid grid;
    grid.bounds = bounds;
    grid.bounds.min_bound_ -= Eigen::Vector3d(pix_width, pix_width, pix_width);
    grid.voxelWidth = pix_width;
    grid.voxelDims = dims;
    grid.voxels.resize(grid.voxelDims[0]*grid.voxelDims[1]*grid.voxelDims[2]);
    grid.calculateDensities(cloud_file.name());

    #if DENSITY_MIN_RAYS > 0
    int X = 1;
    int Y = dims[0];
    int Z = dims[0]*dims[1];
    auto &voxels = grid.voxels;
    VoxelGrid::Voxel neighbours;
    double num_hit_points = 0.0;
    double num_hit_points_unsatisfied = 0.0;
    // This simple 3x3x3 convolution needs to be a bit sneaky to avoid having to double the memory cost.
    // well, not that sneaky, just output into relative cell -1,-1,-1
    
    for (int x = 1; x<grid.voxelDims[0]-1; x++)
    {
      for (int y = 1; y<grid.voxelDims[1]-1; y++)
      {
        for (int z = 1; z<grid.voxelDims[2]-1; z++)
        {
          int ind = grid.getIndex(Eigen::Vector3i(x,y,z));
          if (voxels[ind].numHits > 0)
            num_hit_points++;
          float needed = DENSITY_MIN_RAYS - voxels[ind].numRays;
          VoxelGrid::Voxel corner_vox = voxels[ind - X - Y - Z];
          voxels[ind - X - Y - Z] = voxels[ind]; // move centre up to corner 
          VoxelGrid::Voxel &voxel = voxels[ind - X - Y - Z]; 
          if (needed < 0.0)
            continue;
          neighbours  = voxels[ind-X];
          neighbours += voxels[ind+X];
          neighbours += voxels[ind-Y];
          neighbours += voxels[ind+Y];
          neighbours += voxels[ind-Z];
          neighbours += voxels[ind+Z];
          if (neighbours.numRays >= needed)
          {
            voxel += neighbours * (needed/neighbours.numRays); // add minimal amount to reach DENSITY_MIN_RAYS
            continue;
          }
          voxel += neighbours;
          needed -= neighbours.numRays;

          neighbours  = voxels[ind-X-Y];
          neighbours += voxels[ind-X+Y];
          neighbours += voxels[ind+X-Y];
          neighbours += voxels[ind+X+Y];

          neighbours += voxels[ind-X-Z];
          neighbours += voxels[ind-X+Z];
          neighbours += voxels[ind+X-Z];
          neighbours += voxels[ind+X+Z];

          neighbours += voxels[ind-Y-Z];
          neighbours += voxels[ind-Y+Z];
          neighbours += voxels[ind+Y-Z];
          neighbours += voxels[ind+Y+Z];
          if (neighbours.numRays >= needed)
          {
            voxel += neighbours * (needed/neighbours.numRays); // add minimal amount to reach DENSITY_MIN_RAYS
            continue;
          }
          voxel += neighbours;
          needed -= neighbours.numRays;

          neighbours  = corner_vox;          
          neighbours += voxels[ind-X-Y+Z];          
          neighbours += voxels[ind-X+Y-Z];          
          neighbours += voxels[ind+X-Y-Z];          
          neighbours += voxels[ind-X+Y+Z];          
          neighbours += voxels[ind+X-Y+Z];          
          neighbours += voxels[ind+X+Y-Z];          
          neighbours += voxels[ind+X+Y+Z];     
          if (neighbours.numRays >= needed)
          {
            voxel += neighbours * (needed/neighbours.numRays); // add minimal amount to reach DENSITY_MIN_RAYS
            continue;
          }
          voxel += neighbours;    
          if (voxels[ind].numHits > 0)
            num_hit_points_unsatisfied++;
        }
      }
    }
    double percentage = 100.0*num_hit_points_unsatisfied/num_hit_points;
    std::cout << "Density calculation: " << percentage << "% of voxels had insufficient (<" 
      << DENSITY_MIN_RAYS << ") rays within them" << std::endl;
    if (percentage > 50.0)
    {
      std::cout << "This is high. Consider using a larger pixel size, or a denser cloud, or reducing DENSITY_MIN_RAYS, for consistent results"
        << std::endl;
    }
    else if (percentage < 1.0)
    {
      std::cout << "This is low enough that you could get more fidelity from using a smaller pixel size" << std::endl;
      std::cout << "or more accuracy by increasing DENSITY_MIN_RAYS" << std::endl;
    }
    #endif

    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        double total_density = 0.0;
        for (int z = 0; z< depth; z++)
        {
          Eigen::Vector3i ind;
          ind[axis] = z;
          ind[ax1] = x;
          ind[ax2] = y;
          total_density += grid.voxels[grid.getIndex(ind)].density();
        }
        pixels[x + width * y] = Eigen::Vector4d(total_density, total_density, total_density, total_density);
      }
    }
  }
  else
  {
    auto render = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &, std::vector<ray::RGBA> &colours)
    {
      for (size_t i = 0; i<ends.size(); i++)
      {
        ray::RGBA &colour = colours[i];
        if (colour.alpha == 0)
          continue;
        Eigen::Vector3d col = Eigen::Vector3d(colour.red, colour.green, colour.blue)/255.0;
        Eigen::Vector3d point = style.selectedID() == 3 ? starts[i] : ends[i];
        Eigen::Vector3d pos = (point - bounds.min_bound_) / pix_width;
        Eigen::Vector3i p = (pos).cast<int>();
        int x = p[ax1], y = p[ax2];
        Eigen::Vector4d &pix = pixels[x + width*y];
        switch (style.selectedID())
        {
          case 0: // ends
          case 3: // starts
            if (pos[axis]*dir > pix[3]*dir || pix[3] == 0.0)
              pix = Eigen::Vector4d(col[0], col[1], col[2], pos[axis]);
            break;
          case 1: // mean
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case 2: // sum
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case 4: // rays
          {
            // factor this out (same as in rayalignment.cpp)
            Eigen::Vector3d cloud_start = starts[i];
            Eigen::Vector3d cloud_end = ends[i];
            bounds.clipRay(cloud_start, cloud_end);

            Eigen::Vector3d start = (cloud_start - bounds.min_bound_) / pix_width;
            Eigen::Vector3d end = (cloud_end - bounds.min_bound_) / pix_width;
            Eigen::Vector3d dir = cloud_end - cloud_start;
            Eigen::Vector3d dir_sign(ray::sgn(dir[0]), ray::sgn(dir[1]), ray::sgn(dir[2]));

            Eigen::Vector3i start_index(start.cast<int>());
            Eigen::Vector3i end_index(end.cast<int>());
            double length_sqr = (end_index - start_index).squaredNorm();
            Eigen::Vector3i index = start_index;
            while ((index - start_index).squaredNorm() <= length_sqr + 1e-10)
            {
              if (index[ax1] >= 0 && index[ax1] < width && index[ax2] >= 0 && index[ax2] < height)
                pixels[index[ax1] + width*index[ax2]] += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
              Eigen::Vector3d mid = bounds.min_bound_ + pix_width * Eigen::Vector3d(index[0] + 0.5, index[1] + 0.5, index[2] + 0.5);
              Eigen::Vector3d next_boundary = mid + 0.5 * pix_width * dir_sign;
              Eigen::Vector3d delta = next_boundary - cloud_start;
              Eigen::Vector3d d(delta[0] / dir[0], delta[1] / dir[1], delta[2] / dir[2]);
              if (d[ax1] < d[ax2])
                index[ax1] += dir_sign.cast<int>()[ax1];
              else
                index[ax2] += dir_sign.cast<int>()[ax2];
            }
            break;
          }
          default:
            break;
        }
      }
    };
    ray::readPly(cloud_file.name(), true, render, 0);
  }

  double max_val = 1.0;
  bool is_hdr = image_file.nameExt() == "hdr";
  if (!is_hdr) // limited range, so work out a decent maximum value, I'm using mean + two standard deviations:
  {
    double sum = 0.0;
    double num = 0.0;
    for (auto &pixel: pixels)
    {
      sum += pixel[3];
      if (pixel[3] > 0.0)
        num++;
    }
    double mean = sum / num;
    double sum_sqr = 0.0;
    for (auto &pixel: pixels)
    {
      if (pixel[3] > 0.0)
        sum_sqr += ray::sqr(pixel[3] - mean);
    }
    double standard_deviation = std::sqrt(sum_sqr / num);
    std::cout << "mean: " << mean << ", sd: " << standard_deviation << std::endl;
    max_val = mean + 2.0*standard_deviation;
  }

  std::vector<ray::RGBA> pixel_colours;
  std::vector<float> float_pixel_colours;
  if (is_hdr)
    float_pixel_colours.resize(3 * width * height);
  else
    pixel_colours.resize(width*height);

  for (int x = 0; x < width; x++)
  {
    int indx = flip_x ? width - 1 - x : x;
    for (int y = 0; y < height; y++)
    {
      Eigen::Vector4d colour = pixels[x + width*y];
      Eigen::Vector3d col3d(colour[0], colour[1], colour[2]);
      uint8_t alpha = colour[3] == 0.0 ? 0 : 255;
      switch (style.selectedID())
      {
        case 1: // mean
        case 4: // rays
          col3d /= colour[3];
          break;
        case 2: // sum
        case 5: // density
          col3d /= max_val;
          break;
        case 6: // density_rgb
        {
          if (is_hdr)
            col3d = colour[0] * ray::redGreenBlueSpectrum(std::log10(std::max(1e-6, colour[0])));
          else 
          {
            double shade = colour[0] / max_val;
            col3d = ray::redGreenBlueGradient(shade);
            if (shade < 0.05)
              col3d *= 20.0*shade;
          }
          break;
        }
        default:
          break;
      }
      int ind = indx + width *y;
      if (is_hdr)
      {
        float_pixel_colours[3*ind + 0] = (float)col3d[0];
        float_pixel_colours[3*ind + 1] = (float)col3d[1];
        float_pixel_colours[3*ind + 2] = (float)col3d[2];
      }
      else 
      {
        ray::RGBA col;
        col.red   = uint8_t(std::min(255.0*col3d[0], 255.0));
        col.green = uint8_t(std::min(255.0*col3d[1], 255.0));
        col.blue  = uint8_t(std::min(255.0*col3d[2], 255.0));
        col.alpha = alpha;
        pixel_colours[ind] = col;
      }
    }
  }
  std::cout << "outputting image: " << image_file.name() << std::endl;
  const char *image_name = image_file.name().c_str();
  stbi_flip_vertically_on_write(1);
  if (image_file.nameExt() == "png")
    stbi_write_png(image_name, width, height, 4, (void *)&pixel_colours[0], 4 * width);
  else if (image_file.nameExt() == "bmp")
    stbi_write_bmp(image_name, width, height, 4, (void *)&pixel_colours[0]);
  else if (image_file.nameExt() == "tga")
    stbi_write_tga(image_name, width, height, 4, (void *)&pixel_colours[0]);
  else if (image_file.nameExt() == "png")
    stbi_write_jpg(image_name, width, height, 4, (void *)&pixel_colours[0], 100); // maximal quality
  else if (image_file.nameExt() == "hdr")
    stbi_write_hdr(image_name, width, height, 3, &float_pixel_colours[0]);
  else
  {
    std::cerr << "Error: image format " << image_file.nameExt() << " not known" << std::endl;
    usage();
  }

  return 0;
}
