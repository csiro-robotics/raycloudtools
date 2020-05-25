#include "raydebugdraw.h"

#include "rayunused.h"

namespace ray
{
struct DebugDrawDetail
{
  std::string fixed_frame_id;
};
}  // namespace ray

std::unique_ptr<ray::DebugDraw> ray::DebugDraw::s_instance;

ray::DebugDraw::DebugDraw(const std::string &fixed_frame_id)
  : imp_(new DebugDrawDetail)
{
  imp_->fixed_frame_id = fixed_frame_id;
}

ray::DebugDraw::~DebugDraw() = default;

ray::DebugDraw *ray::DebugDraw::init(int argc, char *argv[], const char *context, bool ros_init)
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

ray::DebugDraw *ray::DebugDraw::instance()
{
  return s_instance.get();
}

void ray::DebugDraw::drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shade, int id)
{
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(point_shade);
  RAYLIB_UNUSED(id);
}

void ray::DebugDraw::drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends)
{
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
}

void ray::DebugDraw::drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                              const std::vector<double> &radii, int id)
{
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(id);
}

void ray::DebugDraw::drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses,
                               const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id)
{
  RAYLIB_UNUSED(centres);
  RAYLIB_UNUSED(poses);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(colour);
  RAYLIB_UNUSED(id);
}
