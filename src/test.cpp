#include <iostream>
#include <so3.h>
#include <sys/time.h>
#include <sophus/so3.hpp>

#include <opencv2/opencv.hpp>

template class Sophus::SO3<double, Eigen::AutoAlign>;
using SO3Type = Sophus::SO3<double>;
using Point = typename Sophus::SO3<double>::Point;

int main()
{
  // Timestamps
  timeval tvStart, tvEnd, tvDiff;
  
  int len = 10000;
  
  double p[3]={.5, .2, .1};
  rot3d::SO3<double> so3(p);
  
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
  {
    so3.rv(2) = 0.1*i;
    so3.vec2mat();
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv:\n" << so3.rv << std::endl;
  
  cv::Mat rv(3,1,CV_64FC1);
  rv.at<double>(0) = 5;
  rv.at<double>(1) = 0.2;
  rv.at<double>(2) = 0.1;
  
  cv::Mat dst;
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
  {
    rv.at<double>(2) = 0.1*i;
    cv::Rodrigues(rv, dst);
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  std::cout << "dst:\n" << dst << std::endl;
  
  rot3d::Matrix3d rot;
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
    rot = SO3Type::exp(Point(0.5, 0.2, 0.1*i)).matrix();
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  std::cout << "exp:\n" << rot << std::endl;
  
  return 0;
}