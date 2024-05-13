// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayrenderer.h"
#include "imagewrite.h"
#include "raycloud.h"
#include "raylib/raylibconfig.h"
#include "rayparse.h"
#if RAYLIB_WITH_TIFF   // build option to support outputting to geotif (.tif) format
#include "geotiffio.h" /* for GeoTIFF */
#include "xtiffio.h"   /* for TIFF */
#endif
#include <fstream>
#include "rayunused.h"

#define DENSITY_MIN_RAYS 10  // larger is more accurate but more blurred. 0 for no adaptive blending

namespace ray
{
#if RAYLIB_WITH_TIFF
// save to geotif format using floating-point per-channel colour data. This function passes a projection file in order
// to geolocate the image
bool writeGeoTiffFloat(const std::string &filename, int x, int y, const float *data, double pixel_width, bool scalar,
                       const std::string &projection_file, double origin_x, double origin_y)
{
  /* Open TIFF descriptor to write GeoTIFF tags */
  TIFF *tif = XTIFFOpen(filename.c_str(), "w");
  if (!tif)
  {
    return false;
  }

  /* Open GTIF Key parser */
  GTIF *gtif = GTIFNew(tif);
  if (!gtif)
  {
    return false;
  }

  if (x < 0 || y < 0)
  {
    std::cerr << "Bad image size: " << x << ", " << y << std::endl;
    return false;
  }
  const uint32_t w = (uint32_t)x;
  const uint32_t h = (uint32_t)y;
  const int channels = scalar ? 2 : 4;

  /* Set up standard TIFF file */
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
  /* set other TIFF tags and write out image ... */
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, scalar ? PHOTOMETRIC_MINISBLACK : PHOTOMETRIC_RGB);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, channels);
  const uint16_t ex_samp[] = { EXTRASAMPLE_ASSOCALPHA };
  TIFFSetField(tif, TIFFTAG_EXTRASAMPLES, 1, &ex_samp);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  // now go line by line to write out the image data
  for (uint32_t row = 0; row < h; row++)
  {
    std::vector<float> pdst(channels * w);

    // moving the data from the dib to a row structure that
    // can be used by the tiff library
    for (uint32_t col = 0; col < w; col++)
    {
      const uint32_t index = 3 * ((h - 1 - row) * w + col);
      const float shade = (data[index + 0] + data[index + 1] + data[index + 2]) / 3.0f;
      if (scalar)
      {
        pdst[2 * col] = shade;
        pdst[2 * col + 1] = shade == 0.0 ? 0.0 : 255.0;
      }
      else
      {
        pdst[4 * col + 0] = data[index + 0];
        pdst[4 * col + 1] = data[index + 1];
        pdst[4 * col + 2] = data[index + 2];
        pdst[4 * col + 3] = shade == 0.0 ? 0.0 : 255.0;
      }
    }

    // now actually write the row data
    TIFFWriteScanline(tif, &pdst[0], row, 0);
  }

  // read in the projection parameters
  if (!projection_file.empty())
  {
    std::ifstream ifs(projection_file.c_str(), std::ios::in);
    if (ifs.fail())
    {
      std::cerr << "cannot open file " << projection_file << std::endl;
      return false;
    }
    std::string line;
    if (!getline(ifs, line))
    {
      return false;
    }
    // the set of keys in the key-value pairs that we are parsing
    const std::vector<std::string> keys = { "+proj", "+ellps", "+datum", "+units", "+lat_0", "+lon_0", "+x_0", "+y_0" };
    std::vector<std::string> values;
    for (const auto &key : keys)
    {
      std::string::size_type found = line.find(key);
      if (found == std::string::npos)  // error checking
      {
        if (key == "+ellps")
        {
          std::cout << "No ellps field found in proj file, setting it equal to the datum." << std::endl;
          values.push_back("");
          continue;
        }
        std::cerr << "Error: cannot find key: " << key << " in the projection file: " << projection_file << std::endl;
        return false;
      }
      // generate the list of values that correspond to the list of keys
      found += key.length() + 1;
      std::string::size_type space = line.find(" ", found);
      if (space == std::string::npos)
        space = line.length() - 1;
      values.push_back(line.substr(found, space - found));
    }
    if (values[1].empty())  // if ellipsoid type not specified, we take it to be the same as the datum
    {
      values[1] = values[2];
    }
    double coord_lat = 0.0;
    if (!values[4].empty())
    {
      coord_lat = std::stod(values[4]);  // latitude
    }
    double coord_long = 0.0;
    if (!values[5].empty())
    {
      coord_long = std::stod(values[5]);  // longitude
    }
    Eigen::Vector2d geo_offset(0, 0);
    if (!values[6].empty())
    {
      geo_offset[0] = std::stod(values[6]);  // offset in m
    }
    if (!values[7].empty())
    {
      geo_offset[1] = std::stod(values[7]);  // offset in m
    }
    std::cout << "geooffset: " << geo_offset << ", geokey: " << values[1] << ", datum: " << values[2]
              << ", coord_long: " << coord_long << std::endl;

    const double scales[3] = { pixel_width, pixel_width, pixel_width };
    TIFFSetField(tif, TIFFTAG_GEOPIXELSCALE, 3, scales);  // set the width of a pixel

    // Set GeoTIFF information
    // We are only supporting a limited set of projection types, so here we assume standard settings
    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelTypeProjected);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);

    GTIFKeySet(gtif, ProjLinearUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);
    GTIFKeySet(gtif, VerticalUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);

    GTIFKeySet(gtif, ProjectionGeoKey, TYPE_SHORT, 1, KvUserDefined);
    GTIFKeySet(gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1, CT_Orthographic);

    // describe the coordinates of the image corners
    const double tiepoints[6] = { 0, 0, 0, origin_x + geo_offset[0], origin_y + geo_offset[1], 0 };
    TIFFSetField(tif, TIFFTAG_GEOTIEPOINTS, 6, tiepoints);

    if (values[1] == "WGS84")  // we support WGS84 by name
    {
      GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);
    }
    else  // all other Geographic type geokeys we parse directly by their number
    {
      std::stringstream ss(values[1]);
      int geokey = 0;
      ss >> geokey;
      if (!ss.fail())  // we are using a direct number here, so
      {
        GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, geokey);
      }
      else
      {
        std::cout << "unknown geographic projection type: " << values[1] << std::endl;
        return false;
      }
    }
    if (values[2] == "WGS84")  // we support the datum type by name
    {
      GTIFKeySet(gtif, GeogGeodeticDatumGeoKey, TYPE_SHORT, 1, Datum_WGS84);
    }
    else if (!values[2].empty())
    {
      std::cout << "unknown geodetic datum: " << values[2] << std::endl;
      return false;
    }

    GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, KvUserDefined);
    GTIFKeySet(gtif, ProjectionGeoKey, TYPE_SHORT, 1, KvUserDefined);
    if (values[0] != "ortho")  // we only support ortho projection
    {
      std::cout << "unknown projection type: " << values[0] << std::endl;
      return false;
    }
    if (values[3] != "m")  // we only support metres as the units
    {
      std::cout << "unknown unit type: " << values[3] << std::endl;
      return false;
    }
    GTIFKeySet(gtif, ProjCenterLongGeoKey, TYPE_DOUBLE, 1, coord_long);
    GTIFKeySet(gtif, ProjCenterLatGeoKey, TYPE_DOUBLE, 1, coord_lat);

    // Store the keys into the TIFF Tags
    GTIFWriteKeys(gtif);
  }

  // get rid of the key parser
  GTIFFree(gtif);

  // save and close the TIFF file descriptor
  XTIFFClose(tif);

  return true;
}
#endif

/// Calculate the surface area per cubic metre within each voxel of the grid. Assuming an unbiased distribution
/// of surface angles.
void DensityGrid::calculateDensities(const std::string &file_name)
{
  auto calculate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                       std::vector<RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end))
      {
        continue; // ray is outside of bounds
      }

      // now walk the voxels
      Eigen::Vector3d dir = end - start;
      const Eigen::Vector3d source = (start - bounds_.min_bound_) / voxel_width_;
      const Eigen::Vector3d target = (end - bounds_.min_bound_) / voxel_width_;
      const double length = dir.norm();
      for (int a = 0; a<3; a++)
      {
        if (dir[a] == 0.0)
        {
          dir[a] = 1e-10; // prevent division by 0
        }
      }
      const double eps = 1e-9;  // to stay away from edge cases
      const double maxDist = (target - source).norm();

      // cached values to speed up the loop below
      Eigen::Vector3i adds;
      Eigen::Vector3d offsets;
      for (int k = 0; k < 3; ++k)
      {
        if (dir[k] > 0.0)
        {
          adds[k] = 1;
          offsets[k] = 0.5;
        }
        else
        {
          adds[k] = -1;
          offsets[k] = -0.5;
        }
      }

      Eigen::Vector3d p = source;  // our moving variable as we walk over the grid
      Eigen::Vector3i inds = p.cast<int>();

      double depth = 0;
      // walk over the grid, one voxel at a time.
      do
      {
        double ls[3] = { (round(p[0] + offsets[0]) - p[0]) / dir[0], (round(p[1] + offsets[1]) - p[1]) / dir[1],
                         (round(p[2] + offsets[2]) - p[2]) / dir[2] };
        int axis = (ls[0] < ls[1] && ls[0] < ls[2]) ? 0 : (ls[1] < ls[2] ? 1 : 2);
        inds[axis] += adds[axis];
        if (inds[axis] < 0 || inds[axis] >= voxel_dims_[axis])
        {
          break;
        }
        double minL = ls[axis] * length;
        depth += minL + eps;
        p = source + dir * (depth / length);
        int index = getIndex(inds);
        if (depth > maxDist)
        {
          double length_in_voxel = minL + maxDist - depth;
          if (colours[i].alpha > 0)
          {
            voxels_[index].addHitRay(static_cast<float>(length_in_voxel * voxel_width_));
          }
          else
          {
            voxels_[index].addMissRay(static_cast<float>(length_in_voxel * voxel_width_));
          }
        }
        else
        {
          voxels_[index].addMissRay(static_cast<float>(minL * voxel_width_));
        }
      } while (depth <= maxDist);
    }
  };
  Cloud::read(file_name, calculate);
}

// This is a form of windowed average over the Moore neighbourhood (3x3x3) window.
void DensityGrid::addNeighbourPriors()
{
#if DENSITY_MIN_RAYS > 0
  const int X = 1;
  const int Y = voxel_dims_[0];
  const int Z = voxel_dims_[0] * voxel_dims_[1];
  DensityGrid::Voxel neighbours;
  double num_hit_points = 0.0;
  double num_hit_points_unsatisfied = 0.0;

  // This simple 3x3x3 convolution needs to be a bit sneaky to avoid having to double the memory cost.
  // well, not that sneaky, we just shift the output -1,-1,-1 for each cell
  for (int x = 1; x < voxel_dims_[0] - 1; x++)
  {
    for (int y = 1; y < voxel_dims_[1] - 1; y++)
    {
      for (int z = 1; z < voxel_dims_[2] - 1; z++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        if (voxels_[ind].numHits() > 0)
          num_hit_points++;
        float needed = DENSITY_MIN_RAYS - voxels_[ind].numRays();
        const DensityGrid::Voxel corner_vox = voxels_[ind - X - Y - Z];
        voxels_[ind - X - Y - Z] = voxels_[ind];  // move centre up to corner
        DensityGrid::Voxel &voxel = voxels_[ind - X - Y - Z];
        if (needed < 0.0)
          continue;
        neighbours = voxels_[ind - X];
        neighbours += voxels_[ind + X];
        neighbours += voxels_[ind - Y];
        neighbours += voxels_[ind + Y];
        neighbours += voxels_[ind - Z];
        neighbours += voxels_[ind + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = voxels_[ind - X - Y];
        neighbours += voxels_[ind - X + Y];
        neighbours += voxels_[ind + X - Y];
        neighbours += voxels_[ind + X + Y];

        neighbours += voxels_[ind - X - Z];
        neighbours += voxels_[ind - X + Z];
        neighbours += voxels_[ind + X - Z];
        neighbours += voxels_[ind + X + Z];

        neighbours += voxels_[ind - Y - Z];
        neighbours += voxels_[ind - Y + Z];
        neighbours += voxels_[ind + Y - Z];
        neighbours += voxels_[ind + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = corner_vox;
        neighbours += voxels_[ind - X - Y + Z];
        neighbours += voxels_[ind - X + Y - Z];
        neighbours += voxels_[ind + X - Y - Z];
        neighbours += voxels_[ind - X + Y + Z];
        neighbours += voxels_[ind + X - Y + Z];
        neighbours += voxels_[ind + X + Y - Z];
        neighbours += voxels_[ind + X + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        if (voxels_[ind].numHits() > 0)
          num_hit_points_unsatisfied++;
      }
    }
  }
  const double percentage = 100.0 * num_hit_points_unsatisfied / num_hit_points;
  std::cout << "Density calculation: " << percentage << "% of voxels had insufficient (<" << DENSITY_MIN_RAYS
            << ") rays within them" << std::endl;
  if (percentage > 50.0)
  {
    std::cout << "This is high. Consider using a larger pixel size, or a denser cloud, or reducing DENSITY_MIN_RAYS, "
                 "for consistent results"
              << std::endl;
  }
  else if (percentage < 1.0)
  {
    std::cout << "This is low enough that you could get more fidelity from using a smaller pixel size" << std::endl;
    std::cout << "or more accuracy by increasing DENSITY_MIN_RAYS" << std::endl;
  }
#endif
}

bool renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction, RenderStyle style,
                 double pix_width, const std::string &image_file, const std::string &projection_file, bool mark_origin,
                 const std::string *const transform_file)
{
  // convert the view direction into useable parameters
  int axis = 0;
  if (view_direction == ViewDirection::Top)
    axis = 2;
  else if (view_direction == ViewDirection::Front || view_direction == ViewDirection::Back)
    axis = 1;
  double dir = 1;
  if (view_direction == ViewDirection::Left || view_direction == ViewDirection::Front)
    dir = -1;
  const bool flip_x = view_direction == ViewDirection::Left || view_direction == ViewDirection::Back;

  // pull out the main image axes (ax1,ax2 are the horiz,vertical axes)
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  // for each view axis (side,top,front = 0,1,2) we need to have an image x axis, and y axis.
  // e.g. x_axes[axis] is the 3D axis to use (x,y,z = 0,1,2) for the image horizontal direction
  const std::array<int, 3> x_axes = { 1, 0, 0 };
  const std::array<int, 3> y_axes = { 2, 2, 1 };
  const int ax1 = x_axes[axis];
  const int ax2 = y_axes[axis];
  const int width = 1 + static_cast<int>(extent[ax1] / pix_width);
  const int height = 1 + static_cast<int>(extent[ax2] / pix_width);
  const int depth = 1 + static_cast<int>(extent[axis] / pix_width);
  std::cout << "outputting " << width << "x" << height << " image" << std::endl;

  try  // there is a possibility of running out of memory here. So provide a helpful message rather than just asserting
  {
    // accumulated colour buffer
    std::vector<Eigen::Vector4d> pixels(width * height, Eigen::Vector4d(0, 0, 0, 0));
    // density calculation is a special case
    if (style == RenderStyle::Density || style == RenderStyle::Density_rgb)
    {
      Eigen::Vector3i dims = (extent / pix_width).cast<int>() + Eigen::Vector3i(1, 1, 1);
#if DENSITY_MIN_RAYS > 0
      dims += Eigen::Vector3i(1, 1, 1);  // so that we have extra space to convolve
#endif
      Cuboid grid_bounds = bounds;
      grid_bounds.min_bound_ -= Eigen::Vector3d(pix_width, pix_width, pix_width);
      DensityGrid grid(grid_bounds, pix_width, dims);

      grid.calculateDensities(cloud_file);

      grid.addNeighbourPriors();

      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          double total_density = 0.0;
          for (int z = 0; z < depth; z++)
          {
            Eigen::Vector3i ind;
            ind[axis] = z;
            ind[ax1] = x;
            ind[ax2] = y;
            total_density += grid.voxels()[grid.getIndex(ind)].density();
          }
          pixels[x + width * y] = Eigen::Vector4d(total_density, total_density, total_density, total_density);
        }
      }
    }
    else  // otherwise we use a common algorithm, specialising on render style only per-ray
    {
      // this lambda expression lets us chunk load the ray cloud file, so we don't run out of RAM
      auto render = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<RGBA> &colours) {
        for (size_t i = 0; i < ends.size(); i++)
        {
          const RGBA &colour = colours[i];
          if (colour.alpha == 0)
            continue;
          const Eigen::Vector3d col = Eigen::Vector3d(colour.red, colour.green, colour.blue) / 255.0;
          const Eigen::Vector3d point = style == RenderStyle::Starts ? starts[i] : ends[i];
          const Eigen::Vector3d pos = (point - bounds.min_bound_) / pix_width;
          const Eigen::Vector3i p = (pos).cast<int>();
          const int x = p[ax1], y = p[ax2];
          // using 4 dimensions helps us to accumulate colours in a greater variety of ways
          Eigen::Vector4d &pix = pixels[x + width * y];
          switch (style)  // render the image according to the chosen style
          {
          case RenderStyle::Ends:
          case RenderStyle::Starts:
          case RenderStyle::Height:
            if (pos[axis] * dir > pix[3] * dir || pix[3] == 0.0)  // using 0.0 precisely as a flag here
            {
              pix = Eigen::Vector4d(col[0], col[1], col[2], pos[axis]);
            }
            break;
          case RenderStyle::Mean:
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case RenderStyle::Sum:
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case RenderStyle::Rays:
          {
            Eigen::Vector3d cloud_start = starts[i];
            Eigen::Vector3d cloud_end = ends[i];
            // clip to within the image (since we exclude unbounded rays from the image bounds)
            if (!bounds.clipRay(cloud_start, cloud_end))
            {
              continue;
            }
            Eigen::Vector3d start = (cloud_start - bounds.min_bound_) / pix_width;
            Eigen::Vector3d end = (cloud_end - bounds.min_bound_) / pix_width;
            const Eigen::Vector3d dir = cloud_end - cloud_start;

            // fast approximate 2D line rendering requires picking the long axis to iterate along
            const bool x_long = std::abs(dir[ax1]) > std::abs(dir[ax2]);
            const int axis_long = x_long ? ax1 : ax2;
            const int axis_short = x_long ? ax2 : ax1;
            const int width_long = x_long ? 1 : width;
            const int width_short = x_long ? width : 1;

            const double gradient = dir[axis_short] / dir[axis_long];
            if (dir[axis_long] < 0.0)
              std::swap(start, end);  // this lets us iterate from low up to high values
            const int start_long = static_cast<int>(start[axis_long]);
            const int end_long = static_cast<int>(end[axis_long]);
            // place a pixel at the height of each midpoint (of the pixel) in the long axis
            const double start_mid_point = 0.5 + static_cast<double>(start_long);
            double height = start[axis_short] + (start_mid_point - start[axis_long]) * gradient;
            for (int l = start_long; l <= end_long; l++, height += gradient)
            {
              const int s = static_cast<int>(height);
              pixels[width_long * l + width_short * s] += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            }
            break;
          }
          default:
            break;
          }
        }
      };
      if (!Cloud::read(cloud_file, render))
        return false;
    }

    double max_val = 1.0;
    double min_val = 0.0;
    const std::string image_ext = getFileNameExtension(image_file);
    const bool is_hdr = image_ext == "hdr" || image_ext == "tif";
    if (!is_hdr)  // limited range, so work out a sensible maximum value, I'm using mean + two standard deviations:
    {
      double sum = 0.0;
      double num = 0.0;
      for (auto &pixel : pixels)
      {
        sum += pixel[3];
        if (pixel[3] > 0.0)
          num++;
      }
      double mean = sum / num;
      double sum_sqr = 0.0;
      for (auto &pixel : pixels)
      {
        if (pixel[3] > 0.0)
          sum_sqr += sqr(pixel[3] - mean);
      }
      const double standard_deviation = std::sqrt(sum_sqr / num);
      max_val = mean + 2.0 * standard_deviation;
      min_val = mean - 2.0 * standard_deviation;
    }

    // The final pixel buffer
    std::vector<RGBA> pixel_colours;
    std::vector<float> float_pixel_colours;
    if (is_hdr)
      float_pixel_colours.resize(3 * width * height);
    else
      pixel_colours.resize(width * height);

    for (int x = 0; x < width; x++)
    {
      const int indx = flip_x ? width - 1 - x : x;  // possible horizontal flip, depending on view direction
      for (int y = 0; y < height; y++)
      {
        const Eigen::Vector4d colour = pixels[x + width * y];
        Eigen::Vector3d col3d(colour[0], colour[1], colour[2]);
        const uint8_t alpha = colour[3] == 0.0 ? 0 : 255;  // 'punch-through' alpha
        switch (style)  // convert to the colour data structure based on the chosen style
        {
        case RenderStyle::Mean:
        case RenderStyle::Rays:
          col3d /= colour[3];  // simple mean
          break;
        case RenderStyle::Height:
        {
          double shade =
            dir == 1.0 ? (colour[3] - min_val) / (max_val - min_val) : (colour[3] - max_val) / (min_val - max_val);
          col3d = Eigen::Vector3d(shade, shade, shade);
          break;
        }
        case RenderStyle::Sum:
        case RenderStyle::Density:
          col3d /= max_val;  // rescale to within limited colour range
          break;
        case RenderStyle::Density_rgb:
        {
          if (is_hdr)
            col3d = colour[0] * redGreenBlueSpectrum(std::log10(std::max(1e-6, colour[0])));
          else
          {
            double shade = colour[0] / max_val;
            col3d = redGreenBlueGradient(shade);
            if (shade < 0.05)
              col3d *= 20.0 * shade;  // this blends the lowest densities down to black
          }
          break;
        }
        default:
          break;
        }
        const int ind = indx + width * y;
        if (is_hdr)
        {
          float_pixel_colours[3 * ind + 0] = (float)col3d[0];
          float_pixel_colours[3 * ind + 1] = (float)col3d[1];
          float_pixel_colours[3 * ind + 2] = (float)col3d[2];
        }
        else
        {
          RGBA col;
          col.red = uint8_t(std::max(0.0, std::min(255.0 * col3d[0], 255.0)));
          col.green = uint8_t(std::max(0.0, std::min(255.0 * col3d[1], 255.0)));
          col.blue = uint8_t(std::max(0.0, std::min(255.0 * col3d[2], 255.0)));
          col.alpha = alpha;
          pixel_colours[ind] = col;
        }
      }
    }
    if (mark_origin)  // an option to mark the lidar origin in the image
    {
      if (pixel_colours.empty())
      {
        std::cout << "warning: mark origin not implemented for hdr images" << std::endl;
      }
      else
      {
        const Eigen::Vector3d pos = -bounds.min_bound_ / pix_width;
        std::cout << "origin pixel location: " << pos[0] << ", " << pos[1] << std::endl;
        const Eigen::Vector3i p = pos.cast<int>();
        const int x = p[ax1], y = p[ax2];
#define DRAW_ARROW  // render the origin as an arrow, which therefore defines the x direction in the lidar frame
#if defined DRAW_ARROW
        for (int xx = x - 2; xx <= x + 10; xx++)
        {
          std::cout << "xx: " << xx << std::endl;
          for (int yy = y - 2; yy <= y + 2; yy++)
          {
            if (xx >= 6 && std::abs(yy - y) > ((x + 10) - xx) / 2)
              continue;
            if (xx >= 0 && xx < width && yy >= 0 && yy < height)
            {
              const int indx = flip_x ? width - 1 - xx : xx;
              RGBA &col = pixel_colours[indx + width * yy];
              col.red = col.green = 0;
              col.blue = col.alpha = 255;
            }
          }
        }
        std::cout << "done" << std::endl;
#endif
        if (x >= 0 && x < width && y >= 0 && y < height)
        {
          const int indx = flip_x ? width - 1 - x : x;  // possible horizontal flip, depending on view direction
          // using 4 dimensions helps us to accumulate colours in a greater variety of ways
          RGBA &col = pixel_colours[indx + width * y];
          col.red = col.blue = col.alpha = 255;
          col.green = 0;
          // we leave alpha alone as it might be needed to indicate the presence of points
        }
        else
        {
          std::cerr << "error: the origin cannot be marked on this image as it is not within the image bounds"
                    << std::endl;
        }
      }
    }
    // option to output the transformation of the image
    if (transform_file != nullptr)
    {
      // Compute transform
      const double scale = pix_width;
      const double translate_x = bounds.min_bound_.x();
      const double translate_y = bounds.max_bound_.y();
      const Eigen::Matrix3d transform =
        (Eigen::Translation2d(translate_x, translate_y) * Eigen::Scaling(scale, -scale)).matrix();

      // Write transform
      std::cout << "outputting transform: " << *transform_file << std::endl;
      std::ofstream ofs;
      ofs.open(*transform_file, std::ios::out);
      if (ofs.fail())
      {
        std::cerr << "Error: cannot open " << *transform_file << " for writing." << std::endl;
        return false;
      }
      ofs << "# Generated by rayrender." << std::endl;
      ofs << "# For a given pixel:" << std::endl;
      ofs << "#   P_pixel = [x_pixel, y_pixel, 1]" << std::endl;
      ofs << "# The position of the centre of the pixel in the ray cloud's frame can be" << std::endl;
      ofs << "# computed as:" << std::endl;
      ofs << "#   P_raycloud = [x_raycloud, y_raycloud, _]" << std::endl;
      ofs << "#   P_raycloud = T * P_pixel" << std::endl;
      ofs << "# Where T is the 3*3 transformation matrix defined in this file:" << std::endl;
      ofs << "#   T = [" << std::endl;
      ofs << "#     transform[0], transform[1], transform[2];" << std::endl;
      ofs << "#     transform[3], transform[4], transform[5];" << std::endl;
      ofs << "#     transform[6], transform[7], transform[8];" << std::endl;
      ofs << "#   ]" << std::endl;
      ofs << "# All z information is lost in rayrender." << std::endl;
      ofs << "transform: [" << std::endl;
      ofs << "  " << transform(0, 0) << ", " << transform(0, 1) << ", " << transform(0, 2) << "," << std::endl;
      ofs << "  " << transform(1, 0) << ", " << transform(1, 1) << ", " << transform(1, 2) << "," << std::endl;
      ofs << "  " << transform(2, 0) << ", " << transform(2, 1) << ", " << transform(2, 2) << "," << std::endl;
      ofs << "]" << std::endl;
      ofs.close();
    }
    std::cout << "outputting image: " << image_file << std::endl;

    // write the image depending on the file format
    const char *image_name = image_file.c_str();
    stbi_flip_vertically_on_write(1);
    if (image_ext == "png")
      stbi_write_png(image_name, width, height, 4, (void *)&pixel_colours[0], 4 * width);
    else if (image_ext == "bmp")
      stbi_write_bmp(image_name, width, height, 4, (void *)&pixel_colours[0]);
    else if (image_ext == "tga")
      stbi_write_tga(image_name, width, height, 4, (void *)&pixel_colours[0]);
    else if (image_ext == "jpg")
      stbi_write_jpg(image_name, width, height, 4, (void *)&pixel_colours[0], 100);  // 100 is maximal quality
    else if (image_ext == "hdr")
      stbi_write_hdr(image_name, width, height, 3, &float_pixel_colours[0]);
#if RAYLIB_WITH_TIFF
    else if (image_ext == "tif")
    {
      // obtain the origin offsets
      const Eigen::Vector3d origin(0, 0, 0);
      const Eigen::Vector3d pos = -(origin - bounds.min_bound_);
      const double x = pos[ax1], y = pos[ax2] + static_cast<double>(height) * pix_width;
      // generate the geotiff file
      writeGeoTiffFloat(image_file, width, height, &float_pixel_colours[0], pix_width, false, projection_file, x, y);
    }
#endif
    else
    {
      std::cerr << "Error: image format " << image_ext << " not supported" << std::endl;
      return false;
    }
  }
  catch (std::bad_alloc const &)  // catch any memory allocation problems in generating large images
  {
    std::cout << "Not enough memory to process the " << width << "x" << height << " image." << std::endl;
    std::cout << "The --pixel_width option can be used to reduce the resolution." << std::endl;
  }
#if !RAYLIB_WITH_TIFF
  RAYLIB_UNUSED(projection_file);
#endif
  return true;
}

}  // namespace ray
