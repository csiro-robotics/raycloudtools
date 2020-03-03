#pragma once
#include "rayutils.h"

namespace RAY
{
bool readPly(const std::string &file, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3d> &normals, int &colourOffset);
bool readPly(const std::string &fileName, std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times);
bool readPlyMesh(const std::string &file, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3i> &indexList);

void writePly(const std::string &fileName, const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times);
void writePlyMesh(const std::string &fileNameRaw, const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> &indexList, bool flipNormals = false);
}