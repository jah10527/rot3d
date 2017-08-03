#ifndef _SO3_HPP
#define _SO3_HPP

//#include <Eigen/Dense>
#include <types.hpp>

namespace rot3d {
template <class Scalar_, int Options = 0>
class SO3{
  public:
  using Scalar = Scalar_;
  
  SO3()
  {
    rv.setZero();
    r.setZero();
  }
  ~SO3(){}
  
  void vec2mat()
  {
    
  }
  
  void mat2vec()
  {
    
  }
  
  Vector3<Scalar> rv;
  Matrix3<Scalar> r;
  
};

}





#endif