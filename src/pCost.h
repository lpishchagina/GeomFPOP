#ifndef PCOST_H
#define PCOST_H
#include <memory>
#include <array>
#include "pSphere.h"
/* Class pCost
 *-------------------------------------------------------------------------------
 * The Gaussian cost for the interval (i,t) in p-dimension. q^i_t = m_{i-1}+\beta + \sum_{j = i}^{t} ||y_j - \mu||^2 = A2\mu^2 + A1\mu + A0
 * A2 = (t-i+1);    A1[k] = -2sum_{j=i}^{t}y^{k}_{j} , k = 0,..,p-1;     A0 = sum_{j=i}^{t} sum_{k=0}^{p-1} (y^k_j)^2
 * mean for the interval (i,t): mean_[k] = (-0.5)*A1[k]/A2;
 * min_value(q^i_t) = A0 + sum_{k=0}^{p-1} ( mean_[k]*(A1[k] + A2 * mean_[k])  );
 * sphere for cost q^i_t : center: c[k] = mean[k] ,  radius : -A0/A2 + sum_{k=0}^{p-1} (sphere->mean[k])^2
 * p - dimension in {2,...,20}
 *-------------------------------------------------------------------------------
 */
template <size_t p>

class pCost {
public:
  double A0; 
  std::array<double, p> A1;
  double A2;

  pCost<p>() {
    A0 = 0.;
    for(size_t k = 0; k < p; k++) {
      A1[k] = 0.;
    } 
    A2 = 0.;
  }
  
  void convertAiToZero() {
    A0 = 0.;
    A2 = 0.;
    for (size_t k = 0; k < p; k++) {
      A1[k] = 0.;
    }
  }
  
  void getMean(std::array<double,p> &mean) {
    for (size_t k = 0; k < p; k++) {
      mean[k] = (-0.5)*A1[k]/A2;
    }
  }
  
  double costInPnt(std::array<double,p> &pnt) {
    double res = A0;
    for(size_t k = 0; k < p; k++) {
      res = res + pnt[k] * (A1[k] + A2 * pnt[k]);
    }
    return(res);
  }
  
  double getMin() {
    std::array<double,p> mean;
    this->getMean(mean);
    double res = this->costInPnt(mean);
  return(res);
  }
  
  void getSphere(pSphere<p> * &sphere) {
    for (size_t k = 0; k < p ; k++) {
      sphere->c[k] = -0.5*A1[k]/A2;
    }
    sphere->r = -A0/A2;//E(Y^2)
    for (size_t k = 0; k < p; k++) {
      sphere->r = sphere->r + (sphere->c[k]) * (sphere->c[k]);
    }
    if (sphere->r > 0.) {
      sphere->r = sqrt(sphere->r);
    } else {
      sphere->r = 0.;
    }
  }
  
  void addCost(pCost * &cost) {
    A0 = A0 + (cost->A0);
    A2 = A2 + (cost->A2);
    for (size_t k = 0; k < p; k++) {
      A1[k] = A1[k]+ (cost->A1[k]);
    }
  }
  
  void addPnt(std::array<double,p>  &pnt) {
    A2 = A2+1;
    for (size_t k = 0; k < p; k++) {
      A0 = A0 + pnt[k]*pnt[k];
      A1[k] = A1[k] -2*pnt[k];
    }
  }
  
  //delete
  void initialize(std::array<double,p> & y){
    this->convertAiToZero();
    for(size_t k = 0; k < p; k++) {
      A0 = A0 + y[k]*y[k];
      A1[k] = A1[k] -2*y[k];
    }
    A2 = A2+1;
  }
  
  void initialize(double a0){
    this->convertAiToZero();
    A0 = a0;
    for (size_t i = 0; i < p; i++) A1[i] = 0.;
    A2 = 0.;
  }
  void initialize(pCost * &cost1, pCost * &cost2){
    A0 = (cost1->A0) - (cost2->A0);
    A2 = (cost1->A2) - (cost2->A2);
    for (size_t k = 0; k < p; k++) {
      A1[k] = (cost1->A1[k]) - (cost2->A1[k]);
    }
  }
  
 };
 
#endif // PCOST
