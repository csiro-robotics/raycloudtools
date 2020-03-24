// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raymesh.h"

namespace RAY
{
bool readPly(const std::string &fileName, std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<double> &intensities, std::vector<uint32_t> &colours);
bool readPlyMesh(const std::string &file, Mesh &mesh);

void writePly(const std::string &fileName, const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times, const std::vector<double> &intensities, const std::vector<uint32_t> &colours);
void writePlyMesh(const std::string &fileNameRaw, const Mesh &mesh, bool flipNormals = false);
}