#include "raydraw.h"
#include <eigen3/Eigen/Geometry>
#include <visualization_msgs/MarkerArray.h>
#include <sensor_msgs/PointCloud2.h>
using namespace RAY;
using namespace std;
using namespace Eigen;

DebugDraw::DebugDraw(const string& fixedFrameId)
{
  fixedFrameId_ = fixedFrameId;
  cloudPublisher[0] = n.advertise<sensor_msgs::PointCloud2>("point_cloud1", 3, true);
  cloudPublisher[1] = n.advertise<sensor_msgs::PointCloud2>("point_cloud2", 3, true);
  linePublisher = n.advertise<visualization_msgs::Marker>("lines", 3, true);
  cylinderPublisher[0] = n.advertise<visualization_msgs::MarkerArray>("cylinders1", 3, true);
  cylinderPublisher[1] = n.advertise<visualization_msgs::MarkerArray>("cylinders2", 3, true);
  ellipsoidPublisher[0] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids", 3, true);
  ellipsoidPublisher[1] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids2", 3, true);
  ellipsoidPublisher[2] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids3", 3, true);
  ellipsoidPublisher[3] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids4", 3, true);
  ellipsoidPublisher[4] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids5", 3, true);
  ellipsoidPublisher[5] = n.advertise<visualization_msgs::MarkerArray>("ellipsoids6", 3, true);
}

void setField2(sensor_msgs::PointField &field, const string &name, int offset, int type, int count)
{
  field.name = name;
  field.offset = offset;
  field.datatype = type;
  field.count = count;
}

void DebugDraw::drawCloud(const vector<Vector3d> &points, const vector<double> &pointShade, int id)
{
  sensor_msgs::PointCloud2 pointCloud;
  pointCloud.header.frame_id = fixedFrameId_;
  pointCloud.header.stamp = ros::Time();
  unsigned int pointStep = 0;
  
  sensor_msgs::PointField x;
  setField2(x, "x", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += sizeof(float);
  pointCloud.fields.push_back(x);
  
  sensor_msgs::PointField y;
  setField2(y, "y", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += sizeof(float);
  pointCloud.fields.push_back(y);
  
  sensor_msgs::PointField z;
  setField2(z, "z", pointStep, sensor_msgs::PointField::FLOAT32, 1);
  pointStep += sizeof(float);
  pointCloud.fields.push_back(z);
  
  sensor_msgs::PointField time;
  bool drawTime = true;
  if(drawTime)
  {
    setField2(time, "time", pointStep, sensor_msgs::PointField::FLOAT64, 1);
    pointStep += sizeof(double);
    pointCloud.fields.push_back(time);
  }
  
  pointCloud.is_bigendian = false;
  pointCloud.is_dense = false;
  pointCloud.point_step = pointStep;
  pointCloud.height = 1; 
  pointCloud.width = points.size();
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
    cloudPublisher[id].publish(pointCloud);
}

void DebugDraw::drawLines(const vector<Vector3d> &starts, const vector<Vector3d> &ends)
{
  visualization_msgs::Marker points;
  points.header.frame_id = fixedFrameId_;
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
  points.color.r = 0.7;
  points.color.g = 0.5;
  points.color.b = 0.3;
  points.color.a = 1.0;
  
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
  linePublisher.publish(points);
}

void DebugDraw::drawCylinders(const vector<Vector3d> &starts, const vector<Vector3d> &ends, const vector<double> &radii, int id)
{
  visualization_msgs::MarkerArray markerArray;
  for (int i = 0; i<(int)starts.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = fixedFrameId_;
    marker.id = i;
    marker.type = marker.CYLINDER;
    marker.action = marker.ADD;
    marker.scale.x = marker.scale.y = 2.0*radii[i]; 
    marker.scale.z = (starts[i] - ends[i]).norm();
    if (id == 0)
    {
      marker.color.a = 1.0;
      marker.color.r = 0.8;
      marker.color.g = 0.7;
      marker.color.b = 0.4;
    }
    else
    {
      marker.color.a = 0.5;
      marker.color.r = 0.5;
      marker.color.g = 0.3;
      marker.color.b = 0.4;
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
  cylinderPublisher[id].publish(markerArray);
}

void DebugDraw::drawEllipsoids(const vector<Vector3d> &centres, const vector<Matrix3d> &poses, const vector<Vector3d> &radii, const Vector3d &colour, int id)
{
  visualization_msgs::MarkerArray markerArray;
  for (int i = 0; i<(int)centres.size(); i++)
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = fixedFrameId_;
    marker.id = i;
    marker.type = marker.SPHERE;
    marker.action = marker.ADD;
    marker.scale.z = 2.0*radii[i][0]; 
    double wid = (radii[i][2] + radii[i][1])/2.0;
    marker.scale.x = 2.0*wid; 
    marker.scale.y = 2.0*wid; 
    marker.color.a = 1.0;
    marker.color.r = colour[0];
    marker.color.g = colour[1];
    marker.color.b = colour[2];

    Vector3d len = poses[i].col(0);

    Vector3d ax = len.cross(Vector3d(0,0,1));
    double angle = atan2(ax.norm(), len[2]);
    Vector3d rotVector = ax.normalized() * -angle;
    Quaterniond q(AngleAxisd(rotVector.norm(), rotVector.normalized()));  
    
       // Quat q(poses[i]);
  //  Quaterniond q(poses[i]);
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

  ellipsoidPublisher[id].publish(markerArray);
}
