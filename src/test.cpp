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
  
  std::cout << "rv:\n" << so3.rv << std::endl;
  so3.cayley();
  so3.cayleyInv();
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv Cayley:\n" << so3.rv << std::endl;
  
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
  {
    so3.rv(2) = 0.1*i;
    so3.cayley();
    so3.cayleyInv();
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  std::cout << "r*r.transpose():\n" << so3.r*so3.r.transpose() << std::endl;
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv:\n" << so3.rv << std::endl;
  so3.cayley();
  std::cout << "rv:\n" << so3.rv << std::endl;
  so3.cayleyInv();
  std::cout << "r:\n" << so3.r << std::endl;
  
  std::cout << "0.5*(r-r.transpose()):\n" << 0.5*(so3.r-so3.r.transpose()) << std::endl;
  
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
  
  so3.r << -0.950146567583153, -6.41765854280073e-05, 0.311803617668748,
     -6.41765854277654e-05, -0.999999917385145, -0.000401386434914383,
      0.311803617668748, -0.000401386434914345, 0.950146484968298;
  so3.cayleyInv();
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "r.trace():\n" << so3.r.trace() << std::endl;
  std::cout << "rv Cayley:\n" << so3.rv << std::endl;
  
  return 0;
}