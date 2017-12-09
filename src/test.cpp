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
    so3.cayley();
//    so3.cayleyInv();
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("Cayley runtime# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
  {
    so3.rv(2) = 0.1*i;
    so3.rodrigues();
//    so3.rodriguesInv();
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("Rodrigues runtime# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
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
  printf("opencv Rodrigues runtime# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  rot3d::Matrix3d rot;
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
    rot = SO3Type::exp(Point(0.5, 0.2, 0.1*i)).matrix();
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("Sophus runtime# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  
  so3.r << -0.950146567583153, -6.41765854280073e-05, 0.311803617668748,
     -6.41765854277654e-05, -0.999999917385145, -0.000401386434914383,
      0.311803617668748, -0.000401386434914345, 0.950146484968298;
  so3.cayleyInv();
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv Cayley:\n" << so3.rv << std::endl;
  
  so3.cayley();
  std::cout << "r:\n" << so3.r << std::endl;
  so3.cayleyInv();
  std::cout << "rv Cayley:\n" << so3.rv << std::endl;
  
  double err = 0;
  rot3d::Matrix3d _r;
  srand(time(0));
  for (int i=0; i<len; i++)
  {
    so3.rv(0) = rand()%2000000-1000000;
    so3.rv(1) = rand()%2000000-1000000;
    so3.rv(2) = rand()%2000000-1000000;
//    std::cout << "rv Cayley:\n" << so3.rv << std::endl;
    so3.cayley();
    _r = so3.r;
//    std::cout << "r:\n" << so3.r << std::endl;
    so3.cayleyInv();
//    std::cout << "rv Cayley:\n" << so3.rv << std::endl << std::endl;
    so3.cayley();
    err += (_r-so3.r).norm();
  }
  std::cout << "rv Cayley absolute mean error:\n" << err/len << std::endl;
  
  err = 0;
  for (int i=0; i<len; i++)
  {
    so3.rv(0) = rand()%200-100;
    so3.rv(1) = rand()%200-100;
    so3.rv(2) = rand()%200-100;
    so3.rv /= so3.rv.norm();
    so3.rv *= M_PI;
//    std::cout << "rv Cayley:\n" << so3.rv << std::endl;
    so3.rodrigues();
//    std::cout << "r:\n" << so3.r << std::endl;
    _r = so3.r;
    so3.cayleyInv();
//    std::cout << "rv Cayley:\n" << so3.rv << std::endl;
    so3.cayley();
//    std::cout << "r:\n" << so3.r << std::endl << std::endl;
    err += (_r-so3.r).norm();
  }
  std::cout << "r Cayley absolute mean error:\n" << err/len << std::endl;

  std::cout << "r*r' Cayley:\n" << so3.r*so3.r.transpose() << std::endl;
  
  
  so3.rv << 1,2,3;
  rot3d::SO3<double> out;
  so3.cayleyMul(so3, out);
  std::cout << "so3.rv:\n" << so3.rv << std::endl;
  std::cout << "rv*rv Cayley:\n" << out.rv << std::endl;
  so3.cayley();
  out.cayley();
  rot3d::Matrix3<double> rr = so3.r*so3.r;
  std::cout << "out Cayley ratation:\n" << out.r << std::endl;
  std::cout << "so3.r*so3.r :\n" << rr << std::endl;
  
  
  return 0;
}