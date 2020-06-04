// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "raydebugdraw.h"

#if RAYLIB_WITH_3ES

#include "rayunused.h"

#include <3esconnectionmonitor.h>
#include <3esmessages.h>
#include <3esserver.h>

#include <shapes/3esshapes.h>

using namespace ray;

namespace ray
{
struct DebugDrawDetail
{
  tes::Server *server = nullptr;
  std::string fixed_frame_id;
};
}  // namespace ray

namespace
{
void updateTes(tes::Server &server)
{
  server.updateTransfers(0);
  server.updateFrame(0.0f);
  tes::ConnectionMonitor *connections = server.connectionMonitor();
  if (connections->mode() == tes::ConnectionMonitor::Synchronous)
  {
    connections->monitorConnections();
  }
  connections->commitConnections();
}
}  // namespace

std::unique_ptr<DebugDraw> DebugDraw::s_instance;

DebugDraw::DebugDraw(const std::string &fixed_frame_id)
  : imp_(new DebugDrawDetail)
{
  tes::ServerInfoMessage info;
  // TODO: (KS) handle difference reference frames.
  tes::initDefaultServerInfo(&info);
  tes::ServerSettings settings(tes::SF_Default | tes::SF_Compress);
  imp_->server = tes::Server::create(settings, &info);
  imp_->fixed_frame_id = fixed_frame_id;

  imp_->server->connectionMonitor()->openFileStream("raylib.3es");

  tes::ConnectionMonitor *connections = imp_->server->connectionMonitor();
  connections->waitForConnection(200);
  if (connections->mode() == tes::ConnectionMonitor::Synchronous)
  {
    connections->monitorConnections();
  }
  connections->commitConnections();
}

DebugDraw::~DebugDraw()
{
  if (imp_->server)
  {
    imp_->server->dispose();
    imp_->server = nullptr;
  }
}

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

void DebugDraw::drawray::Cloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shade, int id)
{
  if (points.empty())
  {
    return;
  }

  // TODO: (KS) use point_shade to apply colour.
  RAYLIB_UNUSED(point_shade);

  std::vector<float> points_single(points.size() * 3);
  tes::Vector3d reference_vertex(points[0].x(), points[0].y(), points[0].z());

  // Convert all vertices to single precision relative to the first vertex.
  size_t next_points_write_index = 0;
  for (size_t i = 0; i < points.size(); ++i)
  {
    points_single[next_points_write_index++] = float(points[i].x() - reference_vertex.x);
    points_single[next_points_write_index++] = float(points[i].y() - reference_vertex.y);
    points_single[next_points_write_index++] = float(points[i].z() - reference_vertex.z);
  }

  tes::MeshShape points_shape(tes::DtPoints, points_single.data(), points.size(), sizeof(*points_single.data()) * 3, id,
                              reference_vertex);
  // Replace any existing shape with the same ID.
  points_shape.setFlags(points_shape.flags() | tes::OFReplace);
  // Send the shape.
  imp_->server->create(points_shape);

  // Send update.
  updateTes(*imp_->server);
}

void DebugDraw::drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends)
{
  if (starts.empty())
  {
    return;
  }

  // Convert to sinle precision.
  tes::Vector3d reference_vertex(starts[0].x(), starts[0].y(), starts[0].z());
  std::vector<float> vertices(starts.size() * 2 * 3);

  // Convert all vertices to single precision relative to the first vertex.
  size_t next_points_write_index = 0;
  for (size_t i = 0; i < starts.size(); ++i)
  {
    vertices[next_points_write_index++] = float(starts[i].x() - reference_vertex.x);
    vertices[next_points_write_index++] = float(starts[i].y() - reference_vertex.y);
    vertices[next_points_write_index++] = float(starts[i].z() - reference_vertex.z);
    vertices[next_points_write_index++] = float(ends[i].x() - reference_vertex.z);
    vertices[next_points_write_index++] = float(ends[i].y() - reference_vertex.z);
    vertices[next_points_write_index++] = float(ends[i].z() - reference_vertex.z);
  }

  // TODO: (KS) there's no as there is for points, so this will create a transient shape (single update).
  tes::MeshShape lines_shape(tes::DtLines, vertices.data(), vertices.size() / 3, sizeof(*vertices.data()) * 3, 0u,
                             reference_vertex);
  lines_shape.setColour(tes::Colour(178, 128, 76));
  // Replace any existing shape with the same ID.
  lines_shape.setFlags(lines_shape.flags() | tes::OFReplace);
  // Send the shape.
  imp_->server->create(lines_shape);

  // Send update.
  updateTes(*imp_->server);
}

void DebugDraw::drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                              const std::vector<double> &radii, int id)
{
  if (starts.empty())
  {
    return;
  }

  // TODO: (KS) handle requests for drawing more than tes::MultiShape::ShapeCountLimit items.
  tes::Vector3d reference_pos(starts[0].x(), starts[0].y(), starts[0].z());
  std::vector<tes::Cylinder> cylinders(starts.size());
  std::vector<tes::Shape *> shape_ptrs(starts.size());

  const tes::Colour colour = (id == 0) ? tes::Colour(204, 178, 102) : tes::Colour(128, 76, 102);
  for (size_t i = 0; i < starts.size(); ++i)
  {
    const tes::Vector3d start(starts[i].x(), starts[i].y(), starts[i].z());
    const tes::Vector3d end(ends[i].x(), ends[i].y(), ends[i].z());
    tes::Vector3d cylinder_axis = end - start;
    const tes::Vector3d cylinder_centre = 0.5 * (end + start) - reference_pos;

    const double length = cylinder_axis.normalise();
    cylinders[i] = tes::Cylinder(id, cylinder_centre, cylinder_axis, float(radii[i]), float(length));
    cylinders[i].setColour(colour);

    shape_ptrs[i] = &cylinders[i];
  }


  tes::MultiShape cylinders_multi_shape(shape_ptrs.data(), shape_ptrs.size(), reference_pos);
  // Replace any existing shape with the same ID.
  cylinders_multi_shape.setFlags(cylinders_multi_shape.flags() | tes::OFReplace);
  // Send the shape.
  imp_->server->create(cylinders_multi_shape);

  // Send update.
  updateTes(*imp_->server);
}

void DebugDraw::drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses,
                               const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id)
{
  if (centres.empty())
  {
    return;
  }

  // TODO: (KS) handle requests for drawing more than tes::MultiShape::ShapeCountLimit items.
  std::vector<tes::Sphere> elliptoids(centres.size());
  std::vector<tes::Shape *> shape_ptrs(centres.size());

  const tes::Colour tes_colour(float(colour.x()), float(colour.y()), float(colour.z()));
  for (size_t i = 0; i < centres.size(); ++i)
  {
    const tes::Vector3d centre(centres[i].x(), centres[i].y(), centres[i].z());
    const tes::Vector3d scale(radii[i].x(), radii[i].y(), radii[i].z());
    const tes::Matrix3d pose(poses[i](0, 0), poses[i](0, 1), poses[i](0, 2), poses[i](1, 0), poses[i](1, 1),
                             poses[i](1, 2), poses[i](2, 0), poses[i](2, 1), poses[i](2, 2));


    elliptoids[i] = tes::Sphere(id, centre);
    elliptoids[i].setRotation(rotationToQuaternion(pose));
    elliptoids[i].setScale(scale);
    elliptoids[i].setColour(tes_colour);

    shape_ptrs[i] = &elliptoids[i];
  }


  tes::MultiShape elliptoids_multi_shape(shape_ptrs.data(), shape_ptrs.size());
  // Replace any existing shape with the same ID.
  elliptoids_multi_shape.setFlags(elliptoids_multi_shape.flags() | tes::OFReplace);
  // Send the shape.
  imp_->server->create(elliptoids_multi_shape);

  // Send update.
  updateTes(*imp_->server);
}

#endif  // RAYLIB_WITH_3ES
