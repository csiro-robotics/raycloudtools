#include "raydebugdraw.h"

#include "rayunused.h"

using namespace ray;

namespace ray
{
struct DebugDrawDetail
{
  std::string fixed_frame_id;
};
}  // namespace ray

std::unique_ptr<DebugDraw> DebugDraw::s_instance;

DebugDraw::DebugDraw(const std::string &fixed_frame_id)
  : imp_(new DebugDrawDetail)
{
  imp_->fixed_frame_id = fixed_frame_id;
}

DebugDraw::~DebugDraw() = default;

DebugDraw *DebugDraw::init(int argc, char *argv[], const char *context, bool ros_init)
{
  RAYLIB_UNUSED(argc);
  RAYLIB_UNUSED(argv);
  RAYLIB_UNUSED(context);
  RAYLIB_UNUSED(ros_init);
  if (!s_instance)
  {
    s_instance = std::make_unique<DebugDraw>();
  }

  return instance();
}

DebugDraw *DebugDraw::instance()
{
  return s_instance.get();
}

void DebugDraw::drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shade, int id)
{
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(point_shade);
  RAYLIB_UNUSED(id);
}

void DebugDraw::drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends)
{
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
}

void DebugDraw::drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                              const std::vector<double> &radii, int id)
{
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(id);
}

void DebugDraw::drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses,
                               const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id)
{
  RAYLIB_UNUSED(centres);
  RAYLIB_UNUSED(poses);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(colour);
  RAYLIB_UNUSED(id);
}
