#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen; 
int main()
{
  // Timestamps
  timeval tvStart, tvEnd, tvDiff;
  
  int len = 1000000;
  
  Vector3d rv;
  rv << 0.2, 0.4, 0.6;
  std::cout << "rv=\n" << rv << std::endl; 
  Matrix3d A;
  A << 0, -rv(2), rv(1),
      rv(2), 0, -rv(0),
      -rv(1), rv(0),0;
      
  std::cout << "A=\n" << A << std::endl; 
  
  Matrix3d A2 = A*A;
  std::cout << "A2=\n" << A2 << std::endl; 
  
  Matrix3d R = 2*(A+A2)/(1-0.5*A2.trace());
  R(0) += 1;
  R(4) += 1;
  R(8) += 1;
  std::cout << "R=\n" << R << std::endl;
  std::cout << "R*R'=\n" << R*R.transpose() << std::endl;
  
  gettimeofday(&tvStart, NULL);
  for (int i=0;i<len;i++)
  {
    
  }
  gettimeofday(&tvEnd, NULL);
  timersub(&tvEnd, &tvStart, &tvDiff);
  printf("# %ld.%06lds\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  return 0;
}