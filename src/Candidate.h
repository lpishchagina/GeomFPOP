#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "pCost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>
#include<memory>

/*Class candidate
 * -------------------------------------------------------------------------------
 * tau - label of change-point candidate
 * cost_tau_t - value of q^i_t
 * isPruning - if true => remove this candidate
 * p - dimension in {2,...,20}
 * -------------------------------------------------------------------------------
 */

template <size_t p>

class candidate {
public:
  unsigned int tau;
  
  bool isPruning;
  pCost<p>* cost_tau_t;//value of cost at time t

  candidate():tau(0),  isPruning(false) {
    cost_tau_t = new pCost<p>();
    
  }
  
  candidate(unsigned int t): tau(t), isPruning(false) {
    cost_tau_t = new pCost<p>();
    
  }
  
};
  //-
#endif //CANDIDATE_H
