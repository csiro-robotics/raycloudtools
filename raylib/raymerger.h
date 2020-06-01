#if !defined RAY_MERGER
#define RAY_MERGER
#include "raycloud.h"
#include "rayutils.h"
namespace ray
{  
// Class for merging ray clouds. The merged result goes into 'result', and differences go into the 'differences' cloud.
class RAYLIB_EXPORT Merger
{
public: 
  // merging on a single cloud removes the transients and puts them in the 'differences' cloud.
  void mergeSingleCloud(const Cloud &cloud, const std::string &merge_type, double num_rays, bool colour_cloud);
  // merge a set of clouds into 'result'.
  void mergeMultipleClouds(std::vector<Cloud> &clouds, const std::string &merge_type, double num_rays);
  // 3-way merge of cloud1 and cloud2. Cloud1 and cloud2 are modified to be just the changed rays from base_cloud.
  void mergeThreeWay(const Cloud &base_cloud, Cloud &cloud1, Cloud &cloud2, const std::string &merge_type, double num_rays);

  Cloud &fixedCloud(){ return result_; }
  const Cloud &fixedCloud() const { return result_; }
  Cloud &differenceCloud(){ return differences_; }
  const Cloud &differenceCloud() const { return differences_; }
private:
  Cloud result_;
  Cloud differences_;
};
}
#endif