// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymerger.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/rayprogressthread.h"
#include "raylib/raythreads.h"
#include "raylib/raycloudwriter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Combines multiple ray clouds. Clouds are not moved but rays are omitted in the combined cloud according to the merge type specified." << std::endl;
  std::cout << "Outputs the combined cloud and the residual cloud of differences." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycombine all raycloud1 raycloud2 ... raycloudN   - concatenate all the rays in the _combined.ply cloud ('all' is optional)" << std::endl;
  std::cout << "           min raycloud1 ... raycloudN 20 rays - combines into one cloud with minimal objects at differences" << std::endl;
  std::cout << "                                                 20 is the number of pass through rays to define " << std::endl;
  std::cout << "           max    - maximal objects included. This is a form of volume intersection (rather than min: union)." << std::endl;
  std::cout << "           oldest - keeps the oldest geometry when there is a difference in later ray clouds." << std::endl;
  std::cout << "           newest - uses the newest geometry when there is a difference in newer ray clouds." << std::endl;
  std::cout << "           order  - conflicts are resolved in argument order, with the first taking priority." << std::endl;
  std::cout << "raycombine basecloud min raycloud1 raycloud2 20 rays - 3-way merge, choses the changed geometry (from basecloud) at any differences. " << std::endl;
  std::cout << "                                                       For merge conflicts it uses the specified merge type." << std::endl;
  std::cout << "        --output raycloud_combined.ply               - optionally specify the output file name." << std::endl;
  // clang-format on
  exit(exit_code);
}

// Combines multiple clouds together
int rayCombine(int argc, char *argv[])
{
  ray::KeyChoice merge_type({ "min", "max", "oldest", "newest", "order" });
  ray::FileArgumentList cloud_files(2);
  ray::DoubleArgument num_rays(0.0, 100.0);
  ray::TextArgument rays_text("rays"), all_text("all");

  // Below: false = allow unusual file extensions, for auto-merging, which occurs on non-standard temporary file names
  ray::FileArgument base_cloud(false), cloud_1(false), cloud_2(false), output_file(false);
  ray::OptionalKeyValueArgument output("output", 'o', &output_file);

  // three-way merge option
  bool standard_format = ray::parseCommandLine(argc, argv, { &merge_type, &cloud_files, &num_rays, &rays_text }, { &output });
  bool concatenate_all = ray::parseCommandLine(argc, argv, { &all_text, &cloud_files }, { &output });
  bool threeway = ray::parseCommandLine(
    argc, argv, { &base_cloud, &merge_type, &cloud_1, &cloud_2, &num_rays, &rays_text }, { &output });
  bool threeway_concatenate =
    ray::parseCommandLine(argc, argv, { &base_cloud, &all_text, &cloud_1, &cloud_2 }, { &output });
  if (!standard_format && !concatenate_all && !threeway && !threeway_concatenate)
  {
    concatenate_all = ray::parseCommandLine(argc, argv, { &cloud_files }, { &output }); // a bit more ambiguous, so only try if the other formats failed
    if (!concatenate_all)
    {
      usage();
    }
  }

  // we know there is at least one file, as we specified a minimum number in FileArgumentList
  std::string file_stub =
    (threeway || threeway_concatenate) ? base_cloud.nameStub() : cloud_files.files()[0].nameStub();

  std::vector<ray::Cloud> clouds;
  if (threeway || threeway_concatenate)
  {
    clouds.resize(2);
    if (!clouds[0].load(cloud_1.name(), false))
      usage();
    if (!clouds[1].load(cloud_2.name(), false))
      usage();
  }
  else if (!concatenate_all)
  {
    clouds.resize(cloud_files.files().size());
    for (int i = 0; i < (int)cloud_files.files().size(); i++)
      if (!clouds[i].load(cloud_files.files()[i].name()))
        usage();
  }

  ray::Threads::init();
  ray::MergerConfig config;
  config.voxel_size = 0.0;  // Infer voxel size
  config.num_rays_filter_threshold = num_rays.value();
  config.merge_type = ray::MergeType::Mininum;

  if (merge_type.selectedKey() == "order")
  {
    config.merge_type = ray::MergeType::Order;
  }
  if (merge_type.selectedKey() == "oldest")
  {
    config.merge_type = ray::MergeType::Oldest;
  }
  if (merge_type.selectedKey() == "newest")
  {
    config.merge_type = ray::MergeType::Newest;
  }
  if (merge_type.selectedKey() == "min")
  {
    config.merge_type = ray::MergeType::Mininum;
  }
  if (merge_type.selectedKey() == "max")
  {
    config.merge_type = ray::MergeType::Maximum;
  }
  if (threeway_concatenate || concatenate_all)
  {
    config.merge_type = ray::MergeType::All;
  }
  std::string combined_file = output.isSet() ? output_file.name() : file_stub + "_combined.ply";
  if (concatenate_all)
  {
    ray::CloudWriter writer;
    if (!writer.begin(combined_file))
      usage();

    // By maintaining these buffers below, we avoid almost all memory fragmentation
    auto concatenate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                        std::vector<double> &times, std::vector<ray::RGBA> &colours) 
    {
      ray::Cloud chunk;
      chunk.starts = starts;
      chunk.ends = ends;
      chunk.colours = colours;
      chunk.times = times;
      writer.writeChunk(chunk);
    };
    for (int i = 0; i < (int)cloud_files.files().size(); i++)
    {
      if (!ray::Cloud::read(cloud_files.files()[i].name(), concatenate))
        usage();
    }
    writer.end();    
    return 0;
  }

  ray::Merger merger(config);
  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);
  ray::Cloud concatenated_cloud;
  const ray::Cloud *fixed_cloud = &merger.fixedCloud();

  if (threeway || threeway_concatenate)
  {
    ray::Cloud base_cloud;
    if (!base_cloud.load(argv[1], false))
      usage();
    merger.mergeThreeWay(base_cloud, clouds[0], clouds[1], &progress);
  }
  else
  {
    merger.mergeMultiple(clouds, &progress);
    std::cout << merger.differenceCloud().rayCount() << " transients, " << merger.fixedCloud().rayCount()
              << " fixed rays." << std::endl;
    merger.differenceCloud().save(file_stub + "_differences.ply");
  }

  progress_thread.join();
  fixed_cloud->save(combined_file);
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayCombine, argc, argv);
}
