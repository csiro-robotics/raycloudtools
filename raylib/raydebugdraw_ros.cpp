#include "raydebugdraw.h"

#if RAYLIB_WITH_ROS

#include <eigen3/Eigen/Geometry>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/MarkerArray.h>

using namespace ray;

namespace ray
{
struct DebugDrawDetail
{
  ros::NodeHandle n;
  ros::Publisher cloud_publisher[2];
  ros::Publisher line_publisher;
  ros::Publisher cylinder_publisher[2];
  ros::Publisher ellipsoid_publisher[6];
  ros::Publisher cylinderPublisher;
  ros::Publisher ringPublisher;
  std::string fixed_frame_id;
};
}  // namespace ray

std::unique_ptr<DebugDraw> DebugDraw::s_instance;

DebugDraw::DebugDraw(const std::string &fixed_frame_id)
  : imp_(new DebugDrawDetail)
{
  imp_->cloud_publisher[0] = imp_->n.advertise<sensor_msgs::PointCloud2>("point_cloud1", 3, true);
  imp_->cloud_publisher[1] = imp_->n.advertise<sensor_msgs::PointCloud2>("point_cloud2", 3, true);
  imp_->line_publisher = imp_->n.advertise<visualization_msgs::Marker>("lines", 3, true);
  imp_->cylinder_publisher[0] = imp_->n.advertise<visualization_msgs::MarkerArray>("cylinders1", 3, true);
  imp_->cylinder_publisher[1] = imp_->n.advertise<visualization_msgs::MarkerArray>("cylinders2", 3, true);
  imp_->ellipsoid_publisher[0] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids", 3, true);
  imp_->ellipsoid_publisher[1] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids2", 3, true);
  imp_->ellipsoid_publisher[2] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids3", 3, true);
  imp_->ellipsoid_publisher[3] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids4", 3, true);
  imp_->ellipsoid_publisher[4] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids5", 3, true);
  imp_->ellipsoid_publisher[5] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids6", 3, true);
  imp_->cylinderPublisher = imp_->n.advertise<visualization_msgs::Marker>("cylinders", 3, true);
  imp_->ringPublisher = imp_->n.advertise<visualization_msgs::Marker>("rings", 3, true);
  imp_->fixed_frame_id = fixed_frame_id;
}

DebugDraw::~DebugDraw() = default;

DebugDraw *DebugDraw::init(int argc, char *argv[], const char *context, bool ros_init)
{
  if (!s_instance)
  {
    if (ros_init)
    {
      ros::init(argc, argv, context);
    }
    s_instance = std::make_unique<DebugDraw>();
  }

  return instance();
}

DebugDraw *DebugDraw::instance()
{
  return s_instance.get();
}

void setField2(sensor_msgs::PointField &field, const std::string &name, int offset, uint8_t type, int count)
{
  field.name = name;
  field.offset = offset;
  field.datatype = type;
  field.count = count;
}

void DebugDraw::drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shade, int id)
{
  sensor_msgs::PointCloud2 point_cloud;
  point_cloud.header.frame_id = imp_->fixed_frame_id;
  point_cloud.header.stamp = ros::Time();
  unsigned int point_step = 0;

  sensor_msgs::PointField x;
  setField2(x, "x", point_step, sensor_msgs::PointField::FLOAT32, 1);
  point_step += unsigned(sizeof(float));
  point_cloud.fields.push_back(x);

  sensor_msgs::PointField y;
  setField2(y, "y", point_step, sensor_msgs::PointField::FLOAT32, 1);
  point_step += unsigned(sizeof(float));
  point_cloud.fields.push_back(y);

  sensor_msgs::PointField z;
  setField2(z, "z", point_step, sensor_msgs::PointField::FLOAT32, 1);
  point_step += unsigned(sizeof(float));
  point_cloud.fields.push_back(z);

  sensor_msgs::PointField time;
  bool draw_time = true;
  if (draw_time)
  {
    setField2(time, "time", point_step, sensor_msgs::PointField::FLOAT64, 1);
    point_step += unsigned(sizeof(double));
    point_cloud.fields.push_back(time);
  }

  point_cloud.is_bigendian = false;
  point_cloud.is_dense = false;
  point_cloud.point_step = point_step;
  point_cloud.height = 1;
  point_cloud.width = unsigned(points.size());
  if (point_cloud.width <= 0)
    return;
  point_cloud.row_step = point_cloud.point_step * point_cloud.width;

  point_cloud.data.resize(point_cloud.row_step);
  for (unsigned int i = 0; i < point_cloud.width; ++i)
  {
    unsigned int point_index = i;
    unsigned int data_index = i * point_cloud.point_step;

    *((float *)&point_cloud.data[data_index + x.offset]) = (float)points[point_index][0];
    *((float *)&point_cloud.data[data_index + y.offset]) = (float)points[point_index][1];
    *((float *)&point_cloud.data[data_index + z.offset]) = (float)points[point_index][2];

    if (draw_time)
      *((double *)&point_cloud.data[data_index + time.offset]) = (float)point_shade[point_index];
  }

  if (point_cloud.width > 0)
    imp_->cloud_publisher[id].publish(point_cloud);
}

void DebugDraw::drawTrunks(const std::vector<Trunk> &trunks)
{
  visualization_msgs::Marker marker;

  marker.header.frame_id = imp_->fixed_frame_id;
  marker.header.stamp = ros::Time::now();
  marker.ns = "cylinder_marker";
  marker.type = visualization_msgs::Marker::LINE_LIST;
  marker.pose.orientation.x = 0.0;
  marker.pose.orientation.y = 0.0;
  marker.pose.orientation.z = 0.0;
  marker.pose.orientation.w = 1.0; 
  marker.color.r = 1.0;
  marker.color.g = 0.5;
  marker.color.b = 0.0;
  marker.color.a = 1.0;
  marker.scale.x = 0.02;
  marker.scale.y = 0.02;
  marker.scale.z = 0.02;
  marker.id = 0;
//  marker.lifetime = ros::Duration();
  marker.action = visualization_msgs::Marker::ADD;

  for (int i = 0; i<(int)trunks.size(); i++)
  {
    for (double angle = 0; angle < 2.0*kPi; angle += kPi/6.0)
    {
      for (int height = 0; height<2; height++)
      {
        double h = (double)height - 0.5;
        for (double dir = -1; dir < 1.1; dir += 2.0)
        {
          for (int j = 0; j<2; j++)
          {
            double ang = angle + (double)j * kPi/6.0;
            geometry_msgs::Point p;
            p.x = trunks[i].centre[0] + std::sin(ang) * (trunks[i].radius + dir*trunks[i].thickness) + h*trunks[i].lean[0]*trunks[i].length;
            p.y = trunks[i].centre[1] + std::cos(ang) * (trunks[i].radius + dir*trunks[i].thickness) + h*trunks[i].lean[1]*trunks[i].length;
            p.z = trunks[i].centre[2] + h*trunks[i].length;
            marker.points.push_back(p);
          }
        }
      }
    }
    for (double angle = 0; angle < 2.0*kPi; angle += kPi/3.0)
    {
      for (int height = 0; height<2; height++)
      {
        double h = (double)height - 0.5;
        double ang = angle;
        geometry_msgs::Point p;
        p.x = trunks[i].centre[0] + std::sin(ang) * (trunks[i].radius + trunks[i].thickness) + h*trunks[i].lean[0]*trunks[i].length;
        p.y = trunks[i].centre[1] + std::cos(ang) * (trunks[i].radius + trunks[i].thickness) + h*trunks[i].lean[1]*trunks[i].length;
        p.z = trunks[i].centre[2] + h*trunks[i].length;
        marker.points.push_back(p);
      }
    }
  }  
  imp_->cylinderPublisher.publish(marker);
}

void DebugDraw::drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends)
{
  visualization_msgs::Marker points;
  points.header.frame_id = imp_->fixed_frame_id;
  points.header.stamp = ros::Time::now();
  points.ns = "lines";
  points.action = visualization_msgs::Marker::ADD;
  points.pose.orientation.x = 0.0;
  points.pose.orientation.y = 0.0;
  points.pose.orientation.z = 0.0;
  points.pose.orientation.w = 1.0;

  points.id = 0;

  points.type = visualization_msgs::Marker::LINE_LIST;

  points.scale.x = 0.01;
  points.scale.y = 0.01;

  // points green
  points.color.r = 0.7f;
  points.color.g = 0.5f;
  points.color.b = 0.3f;
  points.color.a = 1.0f;

  for (unsigned int i = 0; i < starts.size(); i++)
  {
    geometry_msgs::Point p;
    p.x = starts[i][0];
    p.y = starts[i][1];
    p.z = starts[i][2];
    points.points.push_back(p);

    geometry_msgs::Point n;
    n.x = ends[i][0];
    n.y = ends[i][1];
    n.z = ends[i][2];
    points.points.push_back(n);
  }

  // Publish the marker
  imp_->line_publisher.publish(points);
}

void DebugDraw::drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                              const std::vector<double> &radii, int id)
{
  visualization_msgs::MarkerArray marker_array;
  for (int i = 0; i < (int)starts.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = imp_->fixed_frame_id;
    marker.id = i;
    marker.type = marker.CYLINDER;
    marker.action = marker.ADD;
    marker.scale.x = marker.scale.y = 2.0 * radii[i];
    marker.scale.z = (starts[i] - ends[i]).norm();
    if (id == 0)
    {
      marker.color.a = 1.0f;
      marker.color.r = 0.8f;
      marker.color.g = 0.7f;
      marker.color.b = 0.4f;
    }
    else
    {
      marker.color.a = 0.5f;
      marker.color.r = 0.5f;
      marker.color.g = 0.3f;
      marker.color.b = 0.4f;
    }

    Eigen::Vector3d dir = (starts[i] - ends[i]).normalized();
    Eigen::Vector3d ax = dir.cross(Eigen::Vector3d(0, 0, 1));
    double angle = atan2(ax.norm(), dir[2]);
    Eigen::Vector3d rot_vector = ax.normalized() * -angle;
    Eigen::AngleAxisd aa(rot_vector.norm(), rot_vector.normalized());
    Eigen::Quaterniond q(aa);
    marker.pose.orientation.w = q.w();
    marker.pose.orientation.x = q.x();
    marker.pose.orientation.y = q.y();
    marker.pose.orientation.z = q.z();
    Eigen::Vector3d mid = (starts[i] + ends[i]) / 2.0;
    marker.pose.position.x = mid[0];
    marker.pose.position.y = mid[1];
    marker.pose.position.z = mid[2];

    marker_array.markers.push_back(marker);
  }
  imp_->cylinder_publisher[id].publish(marker_array);
}

void DebugDraw::drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses,
                               const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id)
{
  visualization_msgs::MarkerArray marker_array;
  for (int i = 0; i < (int)centres.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = imp_->fixed_frame_id;
    marker.id = i;
    marker.type = marker.SPHERE;
    marker.action = marker.ADD;

    marker.scale.x = 2*radii[i][0];
    marker.scale.y = 2*radii[i][1];
    marker.scale.z = 2*radii[i][2];

    marker.color.a = 1.0;
    marker.color.r = float(colour[0]);
    marker.color.g = float(colour[1]);
    marker.color.b = float(colour[2]);

    Eigen::Quaterniond q(poses[i]);

    q.normalize();
    marker.pose.orientation.w = q.w();
    marker.pose.orientation.x = q.x();
    marker.pose.orientation.y = q.y();
    marker.pose.orientation.z = q.z();
    Eigen::Vector3d mid = centres[i];
    marker.pose.position.x = mid[0];
    marker.pose.position.y = mid[1];
    marker.pose.position.z = mid[2];

    marker_array.markers.push_back(marker);
  }

  imp_->ellipsoid_publisher[id].publish(marker_array);
}

#endif  // RAYLIB_WITH_ROS
