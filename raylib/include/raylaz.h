#pragma once
#include "rayutils.h"

bool RAYreadLas(std::string fileName, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times, std::vector<double> &intensities, int decimation);
void RAYwriteLas(std::string fileName, const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times, const std::vector<double> &intensities);

