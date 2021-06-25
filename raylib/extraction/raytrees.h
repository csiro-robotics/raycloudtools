// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYEXTRACT_TREES_H
#define RAYLIB_RAYEXTRACT_TREES_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"


namespace ray
{

struct Trees
{
  Trees(const Cloud &cloud, bool verbose);
  bool save(const std::string &filename);
};

} // namespace ray
#endif // RAYLIB_RAYEXTRACT_TREES_H