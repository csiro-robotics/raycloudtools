// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"

namespace RAY
{
bool readLas(std::string fileName, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times, std::vector<double> &intensities, int decimation);
}
