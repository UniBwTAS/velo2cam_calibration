/*
  velo2cam_calibration - Automatic calibration algorithm for extrinsic parameters of a stereo camera and a velodyne
  Copyright (C) 2017-2018 Jorge Beltran, Carlos Guindel

  This file is part of velo2cam_calibration.

  velo2cam_calibration is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  velo2cam_calibration is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with velo2cam_calibration.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
  laser_pattern: Find the circle centers in the laser cloud
*/

#define PCL_NO_PRECOMPILE
#define DEBUG 1

#include <dynamic_reconfigure/server.h>
#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <message_filters/synchronizer.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/common/eigen.h>
#include <pcl/common/geometry.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/filter.h>
#include <pcl/filters/impl/passthrough.hpp>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_msgs/ModelCoefficients.h>
#include <pcl_msgs/PointIndices.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include "ros/package.h"

#include <velo2cam_calibration/ClusterCentroids.h>
#include <velo2cam_calibration/LaserConfig.h>
#include "velo2cam_utils.h"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;
using namespace sensor_msgs;

double circles_horizontal_distance_, circles_vertical_distance_, circles_diagonal_, circles_perimeter_, board_width_, board_height_;

class Square
{
  private:
    // SquareVertex _vertices [4];
    pcl::PointXYZ _center;
    std::vector<pcl::PointXYZ> _candidates;
    pcl::PointXYZ _top_left, _top_right, _bot_left, _bot_right;

  public:
    Square(std::vector<pcl::PointXYZ> candidates)
    {
        _candidates = candidates;

        // Compute candidates centroid
        for (int i = 0; i < candidates.size(); ++i)
        {
            _center.x += candidates[i].x;
            _center.y += candidates[i].y;
            _center.z += candidates[i].z;
        }

        _center.x /= candidates.size();
        _center.y /= candidates.size();
        _center.z /= candidates.size();

        std::sort(candidates.begin(), candidates.end(), [](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
            return a.x > b.x;
        });

        _top_left = candidates[0];
        _top_right = candidates[1];
        _bot_left = candidates[2];
        _bot_right = candidates[3];

        if (_top_left.y < _top_right.y)
        {
            std::swap(_top_left, _top_right);
        }

        if (_bot_left.y < _bot_right.y)
        {
            std::swap(_bot_left, _bot_right);
        }
    }

    float distance(pcl::PointXYZ pt1, pcl::PointXYZ pt2)
    {
        return sqrt(pow(pt1.x - pt2.x, 2) + pow(pt1.y - pt2.y, 2) + pow(pt1.z - pt2.z, 2));
    }

    float perimeter()
    {
        float perimeter = 0;
        for (int i = 0; i < 4; ++i)
        {
            perimeter += distance(_candidates[i], _candidates[(i + 1) % 4]);
        }
        return perimeter;
    }

    pcl::PointXYZ at(int i)
    {
        assert(0 <= i && i < 4);
        return _candidates[i];
    }

    bool is_valid()
    {
        // Check if candidates are at +-5cm of target's diagonal/2 to their centroid
        for (int i = 0; i < _candidates.size(); ++i)
        {
            float d = distance(_center, _candidates[i]);
            if (fabs(d - circles_diagonal_ / 2.) > 0.05)
            {
                return false;
            }
        }

        if (DEBUG)
            ROS_INFO_STREAM("[Laser] top_left: " << _top_left << ", top_right: " << _top_right << ", bot_left: " << _bot_left
                                         << ", bot_right: " << _bot_right);

        // Check perimeter?
        if (fabs(perimeter() - circles_perimeter_) > 0.5)
        {
            if (DEBUG)
                ROS_INFO_STREAM("[Laser] Perimeter does not fit: " << perimeter());
            return false;
        }

        // Check width + height?
        if (fabs(distance(_top_left, _top_right) - circles_horizontal_distance_) > 0.1 ||
            fabs(distance(_bot_left, _bot_right) - circles_horizontal_distance_) > 0.1)
        {
            if (DEBUG)
                ROS_INFO_STREAM("[Laser] Width does not fit");
            return false;
        }
        if (fabs(distance(_top_left, _bot_left) - circles_vertical_distance_) > 0.1 ||
            fabs(distance(_top_right, _bot_right) - circles_vertical_distance_) > 0.1)
        {
            if (DEBUG)
                ROS_INFO_STREAM("[Laser] Height does not fit");
            return false;
        }

        return true;
    }
};

// class SquareVertex{
// private:
//   SquareVertex* _parent_neighbour;
//   SquareVertex* _child_neighbour;
//   pcl::PointXYZ _point;

// public:
//   SquareVertex(pcl::PointXYZ pt){
//     _point = pt;
//     _parent_neighbour = NULL;
//     _child_neighbour = NULL;
//   }

//   void set_parent(SquareVertex* vertex){
//     _parent_neighbour = vertex;
//   }

//   void set_child(SquareVertex* vertex){
//     _child_neighbour = vertex;
//   }
// };

void comb(int N, int K, std::vector<std::vector<int>>& groups)
{

    int upper_factorial = 1;
    int lower_factorial = 1;

    for (int i = 0; i < K; i++)
    {
        upper_factorial *= (N - i);
        lower_factorial *= (K - i);
    }
    int n_permutations = upper_factorial / lower_factorial;

    if (DEBUG)
        ROS_INFO("[Laser] %d centers found. Iterating over %d possible set of candidates", N, n_permutations);

    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0);      // N-K trailing 0's

    // print integers and permute bitmask
    do
    {
        std::vector<int> group;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
            {
                group.push_back(i);
            }
        }
        groups.push_back(group);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    assert(groups.size() == n_permutations);
}

ros::Publisher plane_cloud_pub, cumulative_pub, center_pc_pub, centers_pub, pattern_pub, coeff_pub, aux_pub,
    auxpoint_pub, debug_pub;
int nFrames; // Used for resetting center computation
pcl::PointCloud<pcl::PointXYZ>::Ptr cumulative_cloud;

// Dynamic parameters
double edge_threshold_, max_plane_distance_;
double circle_radius_;
int rings_count_;
Eigen::Vector3f plane_normal_axis_;
double plane_max_angle_to_upright_, plane_angle_threshold_to_normal_, plane_max_width_deviation_, plane_max_height_deviation_, plane_distance_threshold_;
int min_num_points_per_ring_, min_num_rings_;
double cluster_size_;
int clouds_proc_ = 0, clouds_used_ = 0;
int min_centers_found_;
bool accumulate_;

void callback(const PointCloud2::ConstPtr& laser_cloud)
{
    if (DEBUG)
        ROS_INFO("[Laser] Processing cloud...");

    pcl::PointCloud<Velodyne::Point>::Ptr velocloud(new pcl::PointCloud<Velodyne::Point>),pattern_cloud(new pcl::PointCloud<Velodyne::Point>);

    clouds_proc_++;

    fromROSMsg(*laser_cloud, *velocloud);

    Velodyne::addRange(*velocloud); // For latter computation of edge detection

    // Plane segmentation
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

    pcl::SACSegmentation<Velodyne::Point> plane_segmentation;
    plane_segmentation.setModelType(pcl::SACMODEL_PARALLEL_PLANE);
    plane_segmentation.setDistanceThreshold(plane_distance_threshold_);
    plane_segmentation.setMethodType(pcl::SAC_RANSAC);
    plane_segmentation.setAxis(Eigen::Vector3f(0,0,1)); // Upright sensor and board required
    plane_segmentation.setEpsAngle(plane_max_angle_to_upright_);
    plane_segmentation.setOptimizeCoefficients(true);
    plane_segmentation.setMaxIterations(1000);

    pcl::ExtractIndices<Velodyne::Point> plane_extract;

    pcl::PointCloud<Velodyne::Point>::Ptr plane_cloud(new pcl::PointCloud<Velodyne::Point>);
    bool valid = false;

    pcl::PointCloud<Velodyne::Point>::Ptr tmp_cloud(new pcl::PointCloud<Velodyne::Point>);

    while (velocloud->size() > min_num_points_per_ring_ * min_num_rings_)
    {
        plane_segmentation.setInputCloud(velocloud);
        plane_segmentation.segment(*inliers, *coefficients);

        valid = true;

        if (inliers->indices.size() < min_num_points_per_ring_ * min_num_rings_)
        {
            valid = false;
            break;
        }

        if (DEBUG)
            ROS_INFO_STREAM("[Laser] Segmented plane with " << inliers->indices.size() << " inliers.");

        // Get inliers
        plane_extract.setInputCloud(velocloud);
        plane_extract.setIndices(inliers);
        plane_extract.setNegative(false);
        plane_extract.filter(*tmp_cloud);

        // Get rings with correct width
        vector<vector<Velodyne::Point*>> rings = Velodyne::getRings(*tmp_cloud, rings_count_);
        int num_rings = 0;
        double min_z = numeric_limits<double>::max();
        double max_z = numeric_limits<double>::min();
        for (vector<vector<Velodyne::Point*>>::iterator ring = rings.begin(); ring < rings.end(); ++ring)
        {
            if (ring->size() < min_num_points_per_ring_)
                continue;

            if (fabs(pcl::euclideanDistance(*(*(ring->begin())), *(*(ring->end() - 1))) - board_width_) > plane_max_width_deviation_)
                continue;

            if((*(ring->begin()))->z < min_z)
                min_z = (*(ring->begin()))->z;

            if((*(ring->begin()))->z > max_z)
                max_z = (*(ring->begin()))->z;

            for (vector<Velodyne::Point*>::iterator pt = ring->begin(); pt < ring->end() - 1; pt++)
            {
                plane_cloud->push_back(*(*pt));
            }

            num_rings++;
        }

        // Check number of rings
        if (DEBUG)
            ROS_INFO_STREAM("[Laser] Plane consists of " << num_rings << " rings.");
        if (num_rings < min_num_rings_)
        {
            if (DEBUG)
                ROS_WARN_STREAM("[Laser] Too few rings.");
            valid = false;
        }

        // Check height
        if(fabs((max_z - min_z) - board_height_) > plane_max_height_deviation_) {
            if (DEBUG)
              ROS_WARN_STREAM("[Laser] Height does not fit.");
            valid = false;
        }

        // Check normal angle
        Eigen::Vector3f plane_normal(coefficients->values[0], coefficients->values[1], coefficients->values[2]);
        double angle = std::atan2(plane_normal.cross(plane_normal_axis_).norm(), plane_normal.dot(plane_normal_axis_));
        if (DEBUG)
            ROS_INFO_STREAM("[Laser] Angle difference between plane normal and desired normal: " << angle);
        if (angle > plane_angle_threshold_to_normal_)
        {
            if (DEBUG)
                ROS_WARN_STREAM("[Laser] Angle does not fit.");
            valid = false;
        }

        if (valid)
        {
            break;
        }
        else
        {
            plane_cloud->clear();

            // Remove inliers
            plane_extract.setInputCloud(velocloud);
            plane_extract.setIndices(inliers);
            plane_extract.setNegative(true);
            plane_extract.filter(*tmp_cloud);
            velocloud->swap(*tmp_cloud);
        }
    }

    if (!valid)
    {
        ROS_ERROR_STREAM("[Laser] Could not estimate plane");
        return;
    }

    sensor_msgs::PointCloud2 plane_cloud_ros;
    pcl::toROSMsg(*plane_cloud, plane_cloud_ros);
    plane_cloud_ros.header = laser_cloud->header;
    plane_cloud_pub.publish(plane_cloud_ros);

    // Copy coefficients to proper object for further filtering
    Eigen::VectorXf coefficients_v(4);
    coefficients_v(0) = coefficients->values[0];
    coefficients_v(1) = coefficients->values[1];
    coefficients_v(2) = coefficients->values[2];
    coefficients_v(3) = coefficients->values[3];

    // Get edges points by range
    vector<vector<Velodyne::Point*>> rings = Velodyne::getRings(*plane_cloud, rings_count_);
    for (vector<vector<Velodyne::Point*>>::iterator ring = rings.begin(); ring < rings.end(); ++ring)
    {
        Velodyne::Point *prev, *succ;
        if (ring->empty())
            continue;

        (*ring->begin())->intensity = 0;
        (*(ring->end() - 1))->intensity = 0;
        for (vector<Velodyne::Point*>::iterator pt = ring->begin() + 1; pt < ring->end() - 1; pt++)
        {
            Velodyne::Point* prev = *(pt - 1);
            Velodyne::Point* succ = *(pt + 1);
            (*pt)->intensity = max(pcl::euclideanDistance(*prev,*(*pt)),pcl::euclideanDistance(*(*pt),*succ));
        }
    }

    pcl::PointCloud<Velodyne::Point>::Ptr edges_cloud(new pcl::PointCloud<Velodyne::Point>);
    for (pcl::PointCloud<Velodyne::Point>::iterator pt = plane_cloud->points.begin(); pt < plane_cloud->points.end(); ++pt)
    {
        if (pt->intensity > edge_threshold_)
        {
            edges_cloud->push_back(*pt);
        }
    }

    if (edges_cloud->points.size() == 0)
    {
        ROS_WARN("[Laser] Could not detect pattern edges.");
        return;
    }

    // Get points belonging to plane in pattern pointcloud
    pcl::SampleConsensusModelPlane<Velodyne::Point>::Ptr dit(
        new pcl::SampleConsensusModelPlane<Velodyne::Point>(edges_cloud));
    std::vector<int> inliers2;
    dit->selectWithinDistance(coefficients_v, max_plane_distance_, inliers2); // 0.05
    pcl::copyPointCloud<Velodyne::Point>(*edges_cloud, inliers2, *pattern_cloud);

    // Remove kps not belonging to circles by coords TODO CHECK IF CAN BE REMOVED
    pcl::PointCloud<pcl::PointXYZ>::Ptr circles_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    vector<vector<Velodyne::Point*>> rings2 = Velodyne::getRings(*pattern_cloud, rings_count_);

    for (vector<vector<Velodyne::Point*>>::iterator ring = rings2.begin(); ring < rings2.end(); ++ring)
    {
        for (vector<Velodyne::Point*>::iterator pt = ring->begin(); pt < ring->end(); ++pt)
        {
            // Velodyne specific info no longer needed for calibration
            // so standard point is used from now on
            pcl::PointXYZ point;
            point.x = (*pt)->x;
            point.y = (*pt)->y;
            point.z = (*pt)->z;
            circles_cloud->push_back(point);
        }
    }

    sensor_msgs::PointCloud2 velocloud_ros2;
    pcl::toROSMsg(*circles_cloud, velocloud_ros2);
    velocloud_ros2.header = laser_cloud->header;
    pattern_pub.publish(velocloud_ros2);

    // Rotate cloud to face pattern plane
    pcl::PointCloud<pcl::PointXYZ>::Ptr xy_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    Eigen::Vector3f xy_plane_normal_vector, floor_plane_normal_vector;
    xy_plane_normal_vector[0] = 0.0;
    xy_plane_normal_vector[1] = 0.0;
    xy_plane_normal_vector[2] = -1.0;

    floor_plane_normal_vector[0] = coefficients->values[0];
    floor_plane_normal_vector[1] = coefficients->values[1];
    floor_plane_normal_vector[2] = coefficients->values[2];

    Eigen::Affine3f rotation = getRotationMatrix(floor_plane_normal_vector, xy_plane_normal_vector);
    pcl::transformPointCloud(*circles_cloud, *xy_cloud, rotation);

    sensor_msgs::PointCloud2 ros_auxpoint;
    pcl::toROSMsg(*xy_cloud, ros_auxpoint);
    ros_auxpoint.header = laser_cloud->header;
    auxpoint_pub.publish(ros_auxpoint);

    pcl::PointCloud<pcl::PointXYZ>::Ptr aux_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointXYZ aux_point;
    aux_point.x = 0;
    aux_point.y = 0;
    aux_point.z = (-coefficients_v(3) / coefficients_v(2));
    aux_cloud->push_back(aux_point);

    pcl::PointCloud<pcl::PointXYZ>::Ptr auxrotated_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::transformPointCloud(*aux_cloud, *auxrotated_cloud, rotation);

    double zcoord_xyplane = auxrotated_cloud->at(0).z;

    // Extract circles
    pcl::ModelCoefficients::Ptr coefficients3(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers3(new pcl::PointIndices);

    // Ransac settings for circle detection
    pcl::SACSegmentation<pcl::PointXYZ> circle_segmentation;
    circle_segmentation.setModelType(pcl::SACMODEL_CIRCLE2D);
    circle_segmentation.setDistanceThreshold(0.02);
    circle_segmentation.setMethodType(pcl::SAC_RANSAC);
    circle_segmentation.setOptimizeCoefficients(true);
    circle_segmentation.setMaxIterations(1000);
    circle_segmentation.setRadiusLimits(circle_radius_ - 0.02, circle_radius_ + 0.02);

    pcl::PointCloud<pcl::PointXYZ>::Ptr copy_cloud(new pcl::PointCloud<pcl::PointXYZ>); // Used for removing inliers
    pcl::copyPointCloud<pcl::PointXYZ>(*xy_cloud, *copy_cloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr circle_cloud(new pcl::PointCloud<pcl::PointXYZ>);   // To store circle points
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_f(new pcl::PointCloud<pcl::PointXYZ>);        // Temp pc used for swaping

    // Force pattern points to belong to computed plane
    for (pcl::PointCloud<pcl::PointXYZ>::iterator pt = copy_cloud->points.begin(); pt < copy_cloud->points.end(); ++pt)
    {
        pt->z = zcoord_xyplane;
    }

    pcl::ExtractIndices<pcl::PointXYZ> extract;
    pcl::PointCloud<pcl::PointXYZ>::Ptr centroid_candidates(new pcl::PointCloud<pcl::PointXYZ>);

    if (DEBUG)
        ROS_INFO("[Laser] Searching for points in cloud of size %lu", copy_cloud->points.size());
    while (copy_cloud->points.size() > 3)
    {
        circle_segmentation.setInputCloud(copy_cloud);
        circle_segmentation.segment(*inliers3, *coefficients3);
        if (inliers3->indices.size() == 0)
        {
            break;
        }

        // Extract the inliers
        extract.setInputCloud(copy_cloud);
        extract.setIndices(inliers3);
        extract.setNegative(false);
        extract.filter(*circle_cloud);

        // Add center point to cloud
        pcl::PointXYZ center;
        center.x = *coefficients3->values.begin();
        center.y = *(coefficients3->values.begin() + 1);
        center.z = zcoord_xyplane;

        centroid_candidates->push_back(center);

        // Remove inliers from pattern cloud to find next circle
        extract.setNegative(true);
        extract.filter(*cloud_f);
        copy_cloud.swap(cloud_f);

        if (DEBUG)
            ROS_INFO("[Laser] Remaining points in cloud %lu", copy_cloud->points.size());
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr rotated_candidates(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::transformPointCloud(*centroid_candidates, *rotated_candidates, rotation.inverse());
    sensor_msgs::PointCloud2 range_ros2;
    pcl::toROSMsg(*rotated_candidates, range_ros2);
    range_ros2.header = laser_cloud->header;
    debug_pub.publish(range_ros2);

    if (DEBUG)
        ROS_INFO("[Laser] Found %lu center candidates", centroid_candidates->size());

    if (centroid_candidates->size() < 4)
    {
        ROS_WARN("[Laser] Not enough centers: %ld", centroid_candidates->size());
        //        return;
    }

    /**
      At this point, circles' center candidates have been computed (found_centers). Now we need to select
      the set of 4 candidates that best fit the calibration target geometry. For that, the following steps
      are followed:
        1) Create a cloud with 4 points representing the exact geometry of the calibration target
        2) For each possible set of 4 points:
          2.1) Rotate cloud to lay in XY plane
          2.2)
    **/
    if (centroid_candidates->size() >= 4)
    {
        std::vector<std::vector<int>> groups;
        comb(centroid_candidates->size(), min_centers_found_, groups);
        double groups_scores[groups.size()]; // -1: invalid; 0-1 normalized score
        for (int i = 0; i < groups.size(); ++i)
        {
            std::vector<pcl::PointXYZ> candidates;
            // Build candidates set
            for (int j = 0; j < groups[i].size(); ++j)
            {
                pcl::PointXYZ center;
                center.x = centroid_candidates->at(groups[i][j]).x;
                center.y = centroid_candidates->at(groups[i][j]).y;
                center.z = centroid_candidates->at(groups[i][j]).z;
                candidates.push_back(center);
            }

            // Compute candidates score
            Square square_candidate(candidates);
            groups_scores[i] = square_candidate.is_valid() ? 1.0 : -1; // -1 when it's not valid, 1 otherwise
        }

        int best_candidate_idx = -1;
        double best_candidate_score = -1;
        for (int i = 0; i < groups.size(); ++i)
        {
            if (best_candidate_score == 1 && groups_scores[i] == 1)
            {
                ROS_ERROR(
                    "More than one set of candidates fit target's geometry. Please, make sure your parameters are "
                    "well set. Exiting callback");
                return;
            }
            if (groups_scores[i] > best_candidate_score)
            {
                best_candidate_score = groups_scores[i];
                best_candidate_idx = i;
            }
        }

        if (best_candidate_idx == -1)
        {
            ROS_WARN("[Laser] Unable to find a candidate set that matches target's geometry");
            return;
        }

        // Build selected centers set TODO WARNING! This replaces center selection algorithm!!
        std::vector<pcl::PointXYZ> selected_centers;
        for (int j = 0; j < groups[best_candidate_idx].size(); ++j)
        {
            selected_centers.push_back(centroid_candidates->at(groups[best_candidate_idx][j]));

            pcl::PointXYZ center_rotated_back =
                pcl::transformPoint(centroid_candidates->at(groups[best_candidate_idx][j]), rotation.inverse());
            center_rotated_back.x = (-coefficients->values[1] * center_rotated_back.y -
                                     coefficients->values[2] * center_rotated_back.z - coefficients->values[3]) /
                                    coefficients->values[0];
            cumulative_cloud->push_back(center_rotated_back); // Build selected centers set TODO WARNING! This replaces
                                                              // center selection algorithm!!
        }

        nFrames++;
        clouds_used_ = nFrames;
    }

    sensor_msgs::PointCloud2 ros_pointcloud;
    pcl::toROSMsg(*cumulative_cloud, ros_pointcloud);
    ros_pointcloud.header = laser_cloud->header;
    cumulative_pub.publish(ros_pointcloud);

    circle_cloud.reset();
    copy_cloud.reset(); // Free memory
    cloud_f.reset();    // Free memory

    pcl_msgs::ModelCoefficients m_coeff;
    pcl_conversions::moveFromPCL(*coefficients, m_coeff);
    m_coeff.header = laser_cloud->header;
    coeff_pub.publish(m_coeff);

    ROS_INFO("[Laser] %d/%d frames: %ld pts in cloud", clouds_used_, clouds_proc_, cumulative_cloud->points.size());

    // Create cloud for publishing centers
    pcl::PointCloud<pcl::PointXYZ>::Ptr centers_cloud(new pcl::PointCloud<pcl::PointXYZ>);

    // Compute circles centers
    getCenterClusters(cumulative_cloud, centers_cloud, cluster_size_, nFrames / 7, nFrames);
    if (centers_cloud->points.size() > 4)
    {
        getCenterClusters(cumulative_cloud, centers_cloud, cluster_size_, 3.0 * nFrames / 4.0, nFrames);
    }

    if (centers_cloud->points.size() == 4)
    {

        sensor_msgs::PointCloud2 ros2_pointcloud;
        pcl::toROSMsg(*centers_cloud, ros2_pointcloud);
        ros2_pointcloud.header = laser_cloud->header;
        center_pc_pub.publish(ros2_pointcloud);

        velo2cam_calibration::ClusterCentroids to_send;
        to_send.header = laser_cloud->header;
        to_send.cluster_iterations = clouds_used_;
        to_send.total_iterations = clouds_proc_;
        to_send.cloud = ros2_pointcloud;

        centers_pub.publish(to_send);
        if (DEBUG)
            ROS_INFO("Pattern centers published");
    }
}

void param_callback(velo2cam_calibration::LaserConfig& config, uint32_t level)
{
    circle_radius_ = config.circle_radius;
    ROS_INFO("New pattern circle radius: %f", circle_radius_);
    plane_distance_threshold_ = config.plane_distance_threshold;
    ROS_INFO("New plane distance threshold: %f", plane_distance_threshold_);
    plane_normal_axis_[0] = config.normal_x;
    plane_normal_axis_[1] = config.normal_y;
    plane_normal_axis_[2] = config.normal_z;
    ROS_INFO("New normal axis for plane segmentation: %f, %f, %f",
             plane_normal_axis_[0],
             plane_normal_axis_[1],
             plane_normal_axis_[2]);
    plane_max_angle_to_upright_ = (config.max_angle_to_upright * M_PI) / 180.0;
    ROS_INFO("New angle threshold to parallel: %f",
             plane_max_angle_to_upright_);
    plane_angle_threshold_to_normal_ = (config.max_angle_to_normal * M_PI) / 180.0;
    ROS_INFO("New angle threshold to normal: %f", plane_angle_threshold_to_normal_);
    plane_max_width_deviation_ = config.max_width_deviation;
    ROS_INFO("New maximal width deviation: %f", plane_max_width_deviation_);
    plane_max_height_deviation_ = config.max_height_deviation;
    ROS_INFO("New maximal height deviation: %f", plane_max_height_deviation_);
    min_num_points_per_ring_ = config.min_num_points_per_ring;
    ROS_INFO("New minimum number of points per ring: %f", min_num_points_per_ring_);
    min_num_rings_ = config.min_num_rings;
    ROS_INFO("New minimum number of rings: %f", min_num_rings_);
    edge_threshold_ = config.edge_threshold;
    ROS_INFO("New edge threshold: %f", edge_threshold_);
    max_plane_distance_ = config.max_plane_distance;
    ROS_INFO("New max plane distance: %f", max_plane_distance_);
    accumulate_ = config.accumulate;
    ROS_INFO("New accumulation option: %f", accumulate_);

    if (config.reset)
    {
        nFrames = 0;
        cumulative_cloud->clear();
        config.reset = false;
        ROS_WARN("Reset cumulative cloud cloud");
    }
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "laser_pattern");
    ros::NodeHandle nh_("~"); // LOCAL
    ros::Subscriber sub = nh_.subscribe("cloud1", 1, callback);
    pcl::console::setVerbosityLevel(pcl::console::L_ALWAYS);

    plane_cloud_pub = nh_.advertise<PointCloud2>("plane_cloud", 1);
    pattern_pub = nh_.advertise<PointCloud2>("pattern_circles", 1);
    auxpoint_pub = nh_.advertise<PointCloud2>("rotated_pattern", 1);
    cumulative_pub = nh_.advertise<PointCloud2>("cumulative_cloud", 1);
    center_pc_pub = nh_.advertise<PointCloud2>("center_pc", 1);
    centers_pub = nh_.advertise<velo2cam_calibration::ClusterCentroids>("centers_cloud", 1);

    debug_pub = nh_.advertise<PointCloud2>("debug", 1);

    coeff_pub = nh_.advertise<pcl_msgs::ModelCoefficients>("plane_model", 1);

    nh_.param("cluster_size", cluster_size_, 0.1);
    nh_.param("min_centers_found", min_centers_found_, 4);
    nh_.param("rings_count", rings_count_, 64);
    nh_.param("circles_horizontal_distance", circles_horizontal_distance_, 0.3);
    nh_.param("circles_vertical_distance", circles_vertical_distance_, 0.3);
    nh_.param("circles_diagonal", circles_diagonal_, 0.4242);
    circles_perimeter_ = 2 * circles_horizontal_distance_ + 2 * circles_vertical_distance_;
    nh_.param("board_width", board_width_, 1.2);
    nh_.param("board_height", board_height_, 0.8);

    nFrames = 0;
    cumulative_cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);

    dynamic_reconfigure::Server<velo2cam_calibration::LaserConfig> server;
    dynamic_reconfigure::Server<velo2cam_calibration::LaserConfig>::CallbackType f;
    f = boost::bind(param_callback, _1, _2);
    server.setCallback(f);

    ros::spin();
    return 0;
}
