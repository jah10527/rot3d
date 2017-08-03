#include <iostream>
#include <so3.hpp>

int main()
{
  rot3d::SO3<float> so3;
  
  
  
  std::cout << "r:\n" << so3.r << std::endl;
  std::cout << "rv:\n" << so3.rv << std::endl;
  
  return 0;
}