#ifndef PRECTANGLE_H
#define PRECTANGLE_H

#include "math.h"
#include "pSphere.h"
#include<memory>
#include <Rcpp.h>

/*Class pRectangle
 * -------------------------------------------------------------------------------
 * pointMin is  left lower corner
 * pointMax is upper right corner
 * p - dimension in {2,...,20}
 * -------------------------------------------------------------------------------
 */
using namespace Rcpp;
using namespace std;

template <size_t p>
class pRectangle {
public:
  unsigned int origine;
  std::array<double,p> pointMin;// left lower corner
  std::array<double,p> pointMax;	//right upper corner

  pRectangle<p>():origine(0){}

  pRectangle<p>(unsigned int tau) : origine{tau} {
    for (size_t k = 0; k < p; k++) {
      pointMin[k] = -INFINITY;
      pointMax[k] = INFINITY;
    }
  }

  pRectangle<p>(unsigned int tau, double * pointMin_, double * pointMax_) : origine{tau} {
    for (size_t k = 0; k < p; k++) {
      pointMin[k] = pointMin_[k];
      pointMax[k] = pointMax_[k];
    }
  }

  //tools-------------------------------------------------------------------------
  double min_ab(double a, double b) { if (a <= b) { return a; } else { return b; } }
  double max_ab(double a, double b) { if (a >= b) { return a; } else { return b; } }
  
  double get_dist(std::array<double,p>  pnt1, std::array<double,p>  pnt2) {
    double res = 0;
    for (size_t k = 0; k < p; k++) {
      res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
    }
    return sqrt(res);
  }
  
  void DoEmptyRect() { 
    pointMin[0] = pointMax[0]; 
  }
  
  bool IsEmptyRect() {
    bool res = true;
    size_t k = 0;
    while(res & (k < p)) {
      res = (pointMin[k] < pointMax[k]);
      k++;
    }
    return(!res);
  }
  
  // closest point of rectangle, relatives to point pnt (sphere center)
  void closestPoint(const std::array<double,p>  pnt, std::array<double,p>  &res) {
    for(size_t k = 0; k < p; k++) {
      if (pnt[k] < pointMin[k]) {
        res[k] = pointMin[k];
      } else if (pnt[k] > pointMax[k]) {
        res[k] = pointMax[k];
      } else {
        res[k] = pnt[k];
      }
    }
  }
  
  // farthest point of block, relatives to point pnt (sphere center)
  void farthestPoint(const std::array<double,p>  pnt, std::array<double,p>  &res) {
    for(size_t k = 0; k < p; k++) {
      if( pnt[k] < pointMin[k]) {
        res[k] = pointMax[k];
      } else if ( pnt[k] > pointMax[k]) {
        res[k] = pointMin[k];
      } else {
        if ((pnt[k] - pointMin[k]) < (pointMax[k] - pnt[k])) {
          res[k] = pointMax[k];
        } else {
          res[k] = pointMin[k];
        }
      }
    }
  }
  
  //smallest distance between rectangle and pnt
  double smallestDistanceToPoint(const std::array<double,p>  pnt) {
    std::array<double,p>  closestP;
    this->closestPoint(pnt, closestP);
    double dif = 0.;
    double res = 0.;
    for (size_t k = 0; k < p; k++) {
      dif = (closestP[k] - pnt[k]);
      res = res + dif*dif;
    }
    return(sqrt(res));
  }
  
  //biggest distance between block and Point (ball center)
  double biggestDistanceToPoint(std::array<double,p> pnt) {
    std::array<double,p>  farthestP;
    this->farthestPoint(pnt, farthestP);
    double dif = 0.;
    double res = 0.;
    for (size_t k = 0; k < p; k++) {
      dif = (farthestP[k] - pnt[k]);
      res = res + dif*dif;
    }
    return(sqrt(res));
  }
  
  bool isInside(pSphere<p>* &sphere) {
    return(sphere->r >= (this->biggestDistanceToPoint(sphere->c)));
  }
  
  bool isOutside(pSphere<p>* &sphere) {
    return(sphere->r <= (this->smallestDistanceToPoint(sphere->c)));
  }
 
  void IntersectionSphere(pSphere<p> *disk) {
    //Rcpp::Rcout<<"rI="<<disk->r<<endl;
    if (!(this->isInside(disk))) {// => possible restriction of the rectangle
      std::array<double,p> closestP;
      this->closestPoint(disk->c, closestP);
      double dif_k = 0.;
      double distance2 = 0.;
      for (size_t k = 0; k < p; k++) {
        dif_k = (closestP[k] - disk->c[k]);
        distance2 = distance2 + dif_k*dif_k;
      }
      if (disk->r <= sqrt(distance2)) {  //=> restriction of the EMPTY rectangle
        this->DoEmptyRect();
      } else {// possible restriction of the rectangle
        double r2 = disk->r * disk->r;
        dif_k = 0.;
        double possiblePmin, possiblePmax;//new boundaries
        double distance_mk;
        for (size_t k = 0; k < p; k++) {
          dif_k = disk->c[k] - closestP[k];
          distance_mk = sqrt(r2 - distance2 + dif_k*dif_k);
          possiblePmin = disk->c[k] - distance_mk;//solution of equations
          possiblePmax = disk->c[k] + distance_mk;
          if (possiblePmin > pointMin[k]){
            pointMin[k] = possiblePmin;
          } 
          if(possiblePmax < pointMax[k]) {
            pointMax[k] = possiblePmax;
          }
        }
      }
    }
  }
  
  void ExclusionSphere(pSphere<p> *&disk) {
   // Rcpp::Rcout<<"rE="<<disk->r<<endl;
    std::array<double,p>  farthestP;
    this->farthestPoint(disk->c, farthestP);
    double distance2 = 0.;
    double dif_k = 0.;
    for(size_t k = 0; k < p; k++) {
      dif_k = (farthestP[k] - disk->c[k]);
      distance2 = distance2 + dif_k*dif_k;
    }
    if( sqrt(distance2) <= disk->r) {  // everything to exclude
      this->DoEmptyRect();
    } else {
      double distance_mk;
      double possiblePmin, possiblePmax;//new boundaries of block
      double r2 = disk->r * disk->r;
      for (size_t k = 0; k < p; k++) {//discriminant
        dif_k = farthestP[k] - disk->c[k];
        distance_mk = sqrt(r2 - distance2 + dif_k*dif_k);
        if (distance_mk > 0) {
          possiblePmax = disk->c[k] + distance_mk;
          possiblePmin = disk->c[k] - distance_mk;
          if ((farthestP[k] == pointMin[k]) && (pointMax[k] <=  possiblePmax)) {
            if (pointMax[k] > possiblePmin) {
              pointMax[k] = possiblePmin;
            }
          }
          if ((farthestP[k] == pointMax[k]) && (pointMin[k] >=  possiblePmin)) {
            if (pointMin[k] < possiblePmax) {
              pointMin[k] = possiblePmax;
            }
          }
        }
      }
    }
  }
};

#endif //PRECTANGLE_H




