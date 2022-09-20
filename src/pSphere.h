#ifndef PSPHERE_H
#define PSPHERE_H

#include<memory>
#include "math.h"

/* Class pSphere
 * -------------------------------------------------------------------------------
 * c - vector  of  the pSphere  center coordinates , r - value of the sphere radius;
 * p - dimension in {2,...,20}
 * -------------------------------------------------------------------------------
 */
template <size_t p>

class pSphere {
public:
  double r;
  std::array<double,p> c;
  
  pSphere<p>(): r(0) {}
  
  pSphere<p>(std::array<double,p> center, double radius) : r{radius} {
    for (size_t  k = 0; k < p; k++) {
      c[k] = center[k];
    } 
  }
  
  double get_dist(std::array<double,p> pnt1, std::array<double,p> pnt2) {
    double res = 0;
    for (size_t  k = 0; k < p; k++){
      res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
    } 
    return sqrt(res);
  }
  
  bool isnotIntersection( pSphere<p>* &sphere) {
    if (r <= 0) return (true);
    double distance2 = 0.;
    double dif = 0.;
    for (size_t k = 0; k < p; k++) {
      dif = c[k] - sphere->c[k];
      distance2 = distance2 + dif*dif;
    }
    return((sphere->r + r) <= sqrt(distance2));
  }
  
  
  bool isInclusion( pSphere<p>  * &sphere) {
    if (r <= 0) return (true);
    double distance = 0.;
    double dif = 0.;
    for(size_t k = 0; k < p; k++) {
      dif = c[k] - sphere->c[k];
      distance = distance + dif*dif;
    }
    return ( (sqrt(distance) + r) <= sphere->r);
  }
  
};
#endif //PSPHERE_H
//------------------------------------------------------------------------------
