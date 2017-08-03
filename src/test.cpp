#include <iostream>
#include <so3.hpp>

#include <opencv2/opencv.hpp>

int main()
{
  float p[3]={.0, .2, .1};
  rot3d::SO3<float> so3(p);
  
  so3.vec2mat();
  
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv:\n" << so3.rv << std::endl;
  
  cv::Mat rv(3,1,CV_64FC1);
  rv.at<double>(0) = 0;
  rv.at<double>(1) = 0.2;
  rv.at<double>(2) = 0.1;
  
  cv::Mat dst;
  cv::Rodrigues(rv, dst);
  
  std::cout << "dst:\n" << dst << std::endl;
  
  return 0;
}