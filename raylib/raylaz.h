// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLAZ_H
#define RAYLIB_RAYLAZ_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

namespace RAY
{
bool RAYLIB_EXPORT readLas(std::string fileName, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times, std::vector<RGBA> &colours, int decimation);
}

#endif // RAYLIB_RAYLAZ_H
