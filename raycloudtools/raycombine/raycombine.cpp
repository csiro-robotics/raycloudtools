// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raydebugdraw.h"
#include "raylib/raymerger.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayprogressthread.h"
#include "raylib/raythreads.h"
#include "raylib/rayparse.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout
    << "Combines multiple ray clouds. Clouds are not moved but rays are omitted in the combined cloud according to "
       "the merge type specified."
    << std::endl;
  std::cout << "Outputs the combined cloud and the residual cloud of differences." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout
    << "raycombine min raycloud1 raycloud2 ... raycloudN 20 rays - combines into one cloud with minimal objects at "
       "differences"
    << std::endl;
  std::cout
    << "                                                           20 is the number of pass through rays to define "
       "a difference"
    << std::endl;
  std::cout
    << "           max    - maximal objects included. This is a form of volume intersection (rather than min: union)."
    << std::endl;
  std::cout << "           oldest - keeps the oldest geometry when there is a difference in later ray clouds."
            << std::endl;
  std::cout << "           newest - uses the newest geometry when there is a difference in newer ray clouds."
            << std::endl;
  std::cout
    << "           all    - combines as a simple concatenation, with all rays remaining (don't include 'xx rays')."
    << std::endl;
  std::cout << "raycombine basecloud min raycloud1 raycloud2 20 rays - 3-way merge, choses the changed geometry (from "
               "basecloud) at any differences. "
            << std::endl;
  std::cout << "For merge conflicts it uses the specified merge type." << std::endl;
  std::cout << "        --output raycloud_combined.ply               - optionally specify the output file name." << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in timevpn-new.csiro.au
int main(int argc, char *argv[])
{
  ray::KeyChoice merge_type({"min", "max", "oldest", "newest"});
  ray::FileArgumentList cloud_files(2);
  ray::DoubleArgument num_rays(0.0, 100.0);
  ray::TextArgument rays_text("rays"), all_text("all");

  // Below: false = allow unusual file extensions, for auto-merging, which occurs on non-standard temporary file names
  ray::FileArgument base_cloud(false), cloud_1(false), cloud_2(false), output_file(false); 
  ray::OptionalKeyValueArgument output("output", 'o', &output_file);

  // three-way merge option
  bool standard_format = ray::parseCommandLine(argc, argv, {&merge_type, &cloud_files, &num_rays, &rays_text}, {&output});
  bool concatenate = ray::parseCommandLine(argc, argv, {&all_text, &cloud_files}, {&output}); 
  bool threeway = ray::parseCommandLine(argc, argv, {&base_cloud, &merge_type, &cloud_1, &cloud_2, &num_rays, &rays_text},
                                                    {&output});
  bool threeway_concatenate = ray::parseCommandLine(argc, argv, {&base_cloud, &all_text, &cloud_1, &cloud_2}, {&output});
  if (!standard_format && !concatenate && !threeway && !threeway_concatenate)
    usage();

  // we know there is at least one file, as we specified a minimum number in FileArgumentList
  std::string file_stub = (threeway || threeway_concatenate) ? base_cloud.nameStub() : cloud_files.files()[0].nameStub(); 

  std::vector<ray::Cloud> clouds;
  if (threeway || threeway_concatenate)
  {
    clouds.resize(2);
    if (!clouds[0].load(cloud_1.name(), false))
      usage();
    if (!clouds[1].load(cloud_2.name(), false))
      usage();
  }
  else
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
  if (concatenate || threeway_concatenate)
  {
    config.merge_type = ray::MergeType::All;
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
  else if (concatenate)
  {
    fixed_cloud = &concatenated_cloud;
    for (auto &cloud : clouds)
    {
      concatenated_cloud.starts.insert(concatenated_cloud.starts.end(), cloud.starts.begin(), cloud.starts.end());
      concatenated_cloud.ends.insert(concatenated_cloud.ends.end(), cloud.ends.begin(), cloud.ends.end());
      concatenated_cloud.times.insert(concatenated_cloud.times.end(), cloud.times.begin(), cloud.times.end());
      concatenated_cloud.colours.insert(concatenated_cloud.colours.end(), cloud.colours.begin(), cloud.colours.end());
    }
  }
  else
  {
    merger.mergeMultiple(clouds, &progress);
    std::cout << merger.differenceCloud().rayCount() << " transients, " << merger.fixedCloud().rayCount()
              << " fixed rays." << std::endl;
    merger.differenceCloud().save(file_stub + "_differences.ply");
  }

  progress_thread.join();

  if (output.isSet())
    fixed_cloud->save(output_file.name());
  else
    fixed_cloud->save(file_stub + "_combined.ply");
  return 0;
}
