#include "raydebugdraw.h"

#include "rayunused.h"

#if RAYLIB_WITH_ROS
#include <eigen3/Eigen/Geometry>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/MarkerArray.h>
#endif

using namespace RAY;
using namespace std;
using namespace Eigen;

namespace RAY
{
struct DebugDrawDetail
{
#if RAYLIB_WITH_ROS
  ros::NodeHandle n;
  ros::Publisher cloudPublisher[2];
  ros::Publisher linePublisher;
  ros::Publisher cylinderPublisher[2];
  ros::Publisher ellipsoidPublisher[6];
  std::string fixedFrameId;
#endif // RAYLIB_WITH_ROS
};
}  // namespace RAY

std::unique_ptr<DebugDraw> DebugDraw::instance_;

DebugDraw::DebugDraw(const string& fixedFrameId)
: imp_(new DebugDrawDetail)
{
  #if !RAYLIB_WITH_ROS
  RAYLIB_UNUSED(fixedFrameId);
  #else
  imp_->cloudPublisher[0] = imp_->n.advertise<sensor_msgs::PointCloud2>("point_cloud1", 3, true);
  imp_->cloudPublisher[1] = imp_->n.advertise<sensor_msgs::PointCloud2>("point_cloud2", 3, true);
  imp_->linePublisher = imp_->n.advertise<visualization_msgs::Marker>("lines", 3, true);
  imp_->cylinderPublisher[0] = imp_->n.advertise<visualization_msgs::MarkerArray>("cylinders1", 3, true);
  imp_->cylinderPublisher[1] = imp_->n.advertise<visualization_msgs::MarkerArray>("cylinders2", 3, true);
  imp_->ellipsoidPublisher[0] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids", 3, true);
  imp_->ellipsoidPublisher[1] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids2", 3, true);
  imp_->ellipsoidPublisher[2] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids3", 3, true);
  imp_->ellipsoidPublisher[3] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids4", 3, true);
  imp_->ellipsoidPublisher[4] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids5", 3, true);
  imp_->ellipsoidPublisher[5] = imp_->n.advertise<visualization_msgs::MarkerArray>("ellipsoids6", 3, true);
  imp_->fixedFrameId = fixedFrameId;
  #endif
}

DebugDraw::~DebugDraw() = default;

DebugDraw *DebugDraw::init(int argc, char *argv[], const char *context, bool rosInit)
{
  if (!instance_)
  {
#if !RAYLIB_WITH_ROS
    RAYLIB_UNUSED(context);
    RAYLIB_UNUSED(argc);
    RAYLIB_UNUSED(argv);
    RAYLIB_UNUSED(rosInit);
#else  // RAYLIB_WITH_ROS
    if (rosInit)
    {
      ros::init(argc, argv, context);
    }
#endif  // RAYLIB_WITH_ROS
    instance_ = std::make_unique<DebugDraw>();
  }

  return instance();
}

DebugDraw *DebugDraw::instance()
{
  return instance_.get();
}

#if RAYLIB_WITH_ROS
void setField2(sensor_msgs::PointField &field, const string &name, int offset, uint8_t type, int count)
{
  field.name = name;
  field.offset = offset;
  field.datatype = type;
  field.count = count;
}
#endif

void DebugDraw::drawCloud(const vector<Vector3d> &points, const vector<double> &pointShade, int id)
{
  #if !RAYLIB_WITH_ROS
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(pointShade);
  RAYLIB_UNUSED(id);
  #else
  sensor_msgs::PointCloud2 pointCloud;
  pointCloud.header.frame_id = 3;
  pointCloud.header.stamp = ros::Time();
  unsigned int pointStep = 0;
  
  sensor_msgs::PointField x;
  setField2(x, "x", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += unsigned(sizeof(float));
  pointCloud.fields.push_back(x);
  
  sensor_msgs::PointField y;
  setField2(y, "y", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += unsigned(sizeof(float));
  pointCloud.fields.push_back(y);
  
  sensor_msgs::PointField z;
  setField2(z, "z", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += unsigned(sizeof(float));
  pointCloud.fields.push_back(z);
  
  sensor_msgs::PointField time;
  bool drawTime = true;
  if(drawTime)
  {
    setField2(time, "time", pointStep, sensor_msgs::PointField::FLOAT64, 1);
    pointStep += unsigned(sizeof(double));
    pointCloud.fields.push_back(time);
  }
  
  pointCloud.is_bigendian = false;
  pointCloud.is_dense = false;
  pointCloud.point_step = pointStep;
  pointCloud.height = 1; 
  pointCloud.width = unsigned(points.size());
  if (pointCloud.width <= 0)
    return;
  pointCloud.row_step = pointCloud.point_step * pointCloud.width; 

  pointCloud.data.resize(pointCloud.row_step);
  for (unsigned int i = 0; i < pointCloud.width; ++i)
  {
    unsigned int pointIndex = i;
    unsigned int dataIndex = i*pointCloud.point_step;
    
    *((float *)&pointCloud.data[dataIndex+x.offset]) = (float)points[pointIndex][0];
    *((float *)&pointCloud.data[dataIndex+y.offset]) = (float)points[pointIndex][1];
    *((float *)&pointCloud.data[dataIndex+z.offset]) = (float)points[pointIndex][2];
    
    if(drawTime)
      *((double *)&pointCloud.data[dataIndex+time.offset]) = (float)pointShade[pointIndex];
  }
  
  if (pointCloud.width > 0)
    imp_->cloudPublisher[id].publish(pointCloud);
  #endif
}

void DebugDraw::drawLines(const vector<Vector3d> &starts, const vector<Vector3d> &ends)
{
  #if !RAYLIB_WITH_ROS
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  #else
  visualization_msgs::Marker points;
  points.header.frame_id = imp_->fixedFrameId;
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
  
  for (unsigned int i = 0; i<starts.size(); i++)
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
  imp_->linePublisher.publish(points);
  #endif
}

void DebugDraw::drawCylinders(const vector<Vector3d> &starts, const vector<Vector3d> &ends, const vector<double> &radii, int id)
{
  #if !RAYLIB_WITH_ROS
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(id);
  #else
  visualization_msgs::MarkerArray markerArray;
  for (int i = 0; i<(int)starts.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = imp_->fixedFrameId;
    marker.id = i;
    marker.type = marker.CYLINDER;
    marker.action = marker.ADD;
    marker.scale.x = marker.scale.y = 2.0*radii[i]; 
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

    Vector3d dir = (starts[i] - ends[i]).normalized();
    Vector3d ax = dir.cross(Vector3d(0,0,1));
    double angle = atan2(ax.norm(), dir[2]);
    Vector3d rotVector = ax.normalized() * -angle;
    AngleAxisd aa(rotVector.norm(), rotVector.normalized());
    Quaterniond q(aa);
    marker.pose.orientation.w = q.w();
    marker.pose.orientation.x = q.x();
    marker.pose.orientation.y = q.y();
    marker.pose.orientation.z = q.z();
    Vector3d mid = (starts[i] + ends[i])/2.0;
    marker.pose.position.x = mid[0];
    marker.pose.position.y = mid[1];
    marker.pose.position.z = mid[2];

    markerArray.markers.push_back(marker);
  }
  imp_->cylinderPublisher[id].publish(markerArray);
  #endif
}

void DebugDraw::drawEllipsoids(const vector<Vector3d> &centres, const vector<Matrix3d> &poses, const vector<Vector3d> &radii, const Vector3d &colour, int id)
{
  #if !RAYLIB_WITH_ROS
  RAYLIB_UNUSED(centres);
  RAYLIB_UNUSED(poses);
  RAYLIB_UNUSED(radii);
  RAYLIB_UNUSED(colour);
  RAYLIB_UNUSED(id);
  #else
  visualization_msgs::MarkerArray markerArray;
  for (int i = 0; i<(int)centres.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = imp_->fixedFrameId;
    marker.id = i;
    marker.type = marker.SPHERE;
    marker.action = marker.ADD;

    double q1 = radii[i][0]/radii[i][1];
    double q2 = radii[i][1]/radii[i][2];

    int ind = q1<q2 ? 0 : 2;
    marker.scale.z = radii[i][ind]; 
    marker.scale.x = radii[i][(ind + 1)%3]; 
    marker.scale.y = radii[i][(ind+2)%3]; 
    marker.color.a = 1.0;
    marker.color.r = float(colour[0]);
    marker.color.g = float(colour[1]);
    marker.color.b = float(colour[2]);

    Vector3d len = poses[i].col(ind);
    Vector3d ax = len.cross(Vector3d(0,0,1));
    double angle = atan2(ax.norm(), len[2]);
    Vector3d rotVector = ax.normalized() * -angle;
    Quaterniond q(AngleAxisd(rotVector.norm(), rotVector.normalized()));  
    
    q.normalize();
    if (!(abs(q.squaredNorm()-1.0) < 0.001))
    {
      cout << "quat " << i << " unnormalized: " << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << endl;
      cout << "poses: " << poses[i] << endl;
    }
    marker.pose.orientation.w = q.w();
    marker.pose.orientation.x = q.x();
    marker.pose.orientation.y = q.y();
    marker.pose.orientation.z = q.z();
    Vector3d mid = centres[i];
    if (!(mid[0] == mid[0]))
    {
      cout << "bad centre: " << mid.transpose() << endl;
    }
    marker.pose.position.x = mid[0];
    marker.pose.position.y = mid[1];
    marker.pose.position.z = mid[2];

    markerArray.markers.push_back(marker);
  }

  imp_->ellipsoidPublisher[id].publish(markerArray);
  #endif
}
