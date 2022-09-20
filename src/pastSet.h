#ifndef PASTSET_H
#define PASTSET_H

#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>
#include "pCost.h"

/*Class pastset
 * -------------------------------------------------------------------------------
 * tau -changepoint candidate
 * pastCost - past cost functions q_{tau-1}^{j}  for change-point candidate tau (j<tau)
 * p - dimension in {2,...,20}
 * -------------------------------------------------------------------------------
 */

template <size_t p>
class pastset {
public:
  unsigned int tau;
  std::list<pSphere<p>*> pastSpheres;
  pastset<p>(unsigned int t): tau(t)  {
    pastSpheres.clear();
  }
};
#endif //PASTSET_H
