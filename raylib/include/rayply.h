#pragma once
#include "rayutils.h"

namespace RAY
{
bool readPly(const std::string &file, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3d> &normals, int &colourOffset);
bool readPlyMesh(const std::string &file, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3i> &indexList);

// Save the polygon file to disk
void writePly(const std::string &fileNameRaw, const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3d> &rgb);
void writePlyMesh(const std::string &fileNameRaw, const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> &indexList, bool flipNormals = false);
}