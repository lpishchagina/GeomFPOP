#ifndef ALGOS_H
#define ALGOS_H

#include "Candidate.h"
#include "pCost.h"
#include "pRectangle.h"
#include "pastSet.h"

#include <Rcpp.h>
#include "math.h"

/*Class Algos
 * -------------------------------------------------------------------------------
 * OP, PELT, GeomFPOP
 * GeomFPOP(rectangle): 
 * intersection : all, last, (last+random) 
 * exclusion : empty, random, all (past spheres of potential candidates), all* (all past spheres)
 * p - dimension in {2,...,20}
 *  GeomFPOP(sphere): 
 *  check 1 : inclusion of Last sphere to Past spheres
 *  check 2 : intersection of Last sphere to Future spheres
 * -------------------------------------------------------------------------------
 */

using namespace Rcpp;
using namespace std;

template <size_t p>
  
class Algos {
public:
  unsigned int N;
  double Penalty;
  double* VectOfCosts;          //UnpenalizedCost = VectOfCosts[n] - Changes.size()*Penalty
  unsigned int* LastChpt;       //vector of the best last changepoints
  std::vector <unsigned int> Changes;
  std::vector <std::vector <double>> SegmentMeans;
  std::vector <unsigned int> NbOfCandidats;
  double UnpenalizedCost;
  std::array<double,p>* Data;
  
  Algos<p>(Rcpp::NumericMatrix data, double penalty) {
    N = (unsigned int)data.ncol();
    Penalty = penalty;
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    Data = new std::array<double, p>[N];
    for (size_t  i = 0; i < N; i++) {
      for(size_t k = 0; k < p; k++) { Data[i][k] = data(k,i);}
    }
  }
  
  ~Algos<p>() {
    delete [] VectOfCosts;
    delete [] LastChpt;
    VectOfCosts = NULL;
    LastChpt = NULL;
    delete [] Data;
    Data = NULL;
  }
  
  Algos<p>(const Algos &cand) {
    N = cand.N;
    Penalty = cand.Penalty;
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    Data = new std::array<double, p>[N];
    for (size_t i = 0; i < N + 1; i++) {VectOfCosts[i] = cand.VectOfCosts[i];}
    for (size_t i = 0; i < N; i++) {
      LastChpt[i] = cand.LastChpt[i];
      for (size_t k = 0; k < p; i++) {Data[i][k] = cand.Data[i][k];}
    }
    Changes = cand.Changes;
    SegmentMeans = cand.SegmentMeans;
    NbOfCandidats = cand.NbOfCandidats;
    UnpenalizedCost = cand.UnpenalizedCost;
  }
  
  //PRUNING-------------------------------------------------------------------//
  void pruning_t(std::list<unsigned int> &changes_t, candidate<p> ** &candidates, unsigned int nbCds) {
    std::list<unsigned int>::iterator iter = changes_t.begin();
    for (size_t i = 0; i < nbCds; i++) {
      if (candidates[*iter]->isPruning) {iter = changes_t.erase(iter);} 
      else { iter++;}
    }
  }
  //RETURN OPTIMAL SEGMENTATION-----------------------------------------------//
  void backtracking() {
    unsigned int chp = N;
    while (chp > 0) {
      Changes.push_back(chp);
      chp = LastChpt[chp-1];
    }
    Changes.push_back(0);
    unsigned int j = 1;
    chp = N - 1;
    while (chp > 0) {
      chp = Changes[j];
      j++;
    }
    reverse(Changes.begin(), Changes.end());//{0,...,N}
    std::vector<double> MeanOneSegment;
    for (unsigned int s = 0; s < Changes.size()-1; s++){
      double coef = Changes[s+1] - Changes[s];
      MeanOneSegment.clear();
      for (unsigned int k = 0; k < p; k++){
        double mean = 0.;
        for (unsigned int count = Changes[s]; count < Changes[s+1]; count++) {mean = mean + Data[count][k];}
        mean = mean/coef;
        MeanOneSegment.push_back(mean);
      }
      SegmentMeans.push_back(MeanOneSegment);
    }
    Changes.pop_back();//remove N
    reverse(Changes.begin(), Changes.end());
    Changes.pop_back();//remove 0
    reverse(Changes.begin(), Changes.end());
    UnpenalizedCost = VectOfCosts[N] - Penalty * (Changes.size());
  }
  //SHOW RESULT---------------------------------------------------------------//
  List ResAlgo(bool showNbCands) {
    List res;
    res["changes"] = Changes;
    res["means"] =  SegmentMeans;
    res["UnpenalizedCost"] = UnpenalizedCost;
    if (showNbCands) {res["NumberOfCandidats"] = NbOfCandidats;}
    return res;
  }
  //RETURN TYPE GEOMFPOP(R-type)----------------------------------------------//
  // for intersection: -1-last, 1-last+random, 10-all
  //for exclusion: 0-empty, 1-random, 10-all
  void nbInterExclGeomFPOP(unsigned int type_algo, int &nbInter, int &nbExcl) {
    if (type_algo == 11) {return;}
    if ((type_algo >= 2) && (type_algo <= 4)) {nbExcl = 10;}
    if ((type_algo >= 5) && (type_algo <= 7)) {nbExcl = 1;}
    if ((type_algo == 2) || (type_algo == 5)|| (type_algo == 8) || (type_algo == 12)) {nbInter = 10;}
    if ((type_algo == 3) || (type_algo == 7) || (type_algo == 9) || (type_algo == 13)) {nbInter = -1;}
    if ((type_algo == 4) || (type_algo == 6) || (type_algo == 10) || (type_algo == 14)) {nbInter = 1;}
    if ((type_algo == 8) || (type_algo == 9) || (type_algo == 10)) {nbExcl = 0;}
  }
  
  //ALGORITHM-----------------------------------------------------------------//
  List algosOP(unsigned int type_algo, bool showNbCands, int nbRandInter, int nbRandExcl){
    srand(time(0));
    //INITIALIZATION------------------------------------//
    candidate<p> ** candidates;
    pRectangle<p> ** blocks;
    pastset<p> ** pastSets;
    pCost<p> *costDif = new pCost<p>();
    VectOfCosts[0] = 0;
    // Initialize all objects
    candidates = new candidate<p> * [N+1];// for all methods
    if (type_algo >= 2 && type_algo <= 14) {//for GeomFPOP
      if (type_algo != 11) { blocks = new pRectangle<p> * [N+1];} //for GeomFPOP(R-type)GeomFPOP(R-type*)
      if (type_algo >= 11 && type_algo <= 14) {pastSets = new pastset<p>* [N+1];}//for GeomFPOP(R-type*), Geom-FPOP(S-type)
    }
    for (size_t i = 0; i < (N+1); i++) {
      candidates[i] = new candidate<p>(i);
      if (type_algo >= 2 && type_algo <= 14) {
        if (type_algo != 11) { blocks[i] = new pRectangle<p>(i);} 
        if (type_algo >= 11 && type_algo <= 14) {pastSets[i] = new pastset<p>(i);}
      }
    } 
    //Initialize list with potential candidates at t
    std::list<unsigned int> changes_t;
    changes_t.push_front(0);
    //number of intersection and exclusion for GeomFPOP
    int nbInter;
    int nbExcl;
    if (type_algo >= 2 && type_algo <=14) {nbInterExclGeomFPOP(type_algo, nbInter, nbExcl);}
    //INITIALIZATION END//
    
    //STEP 1. OPTIMAL PARTITIONING: MIN and ARGMIN//
    pCost<p> * funCtt = new pCost<p>();
    double q_bestTau_t;
    double q_Tau_t;
    //ALGO
    for (size_t t = 0; t < N; t++) {
      // add new data-point in cost functions 
      funCtt->initialize(Data[t]);
      for (std::list<unsigned int>::iterator iter = changes_t.begin(); iter != changes_t.end(); iter++) {
        (candidates[*iter]->cost_tau_t)->addCost(funCtt);
      }
      // min
      q_bestTau_t = (candidates[changes_t.front()]->cost_tau_t)->getMin();
      LastChpt[t] = 0;//best solution of last change at time t
      for (std::list<unsigned int>::iterator iter = changes_t.begin(); iter != changes_t.end(); iter++) {
        q_Tau_t = (candidates[*iter]->cost_tau_t)->getMin();
        if(q_Tau_t <= q_bestTau_t) {
          q_bestTau_t = q_Tau_t;
          LastChpt[t] = candidates[*iter]->tau;
        }
      }
      // store min
      VectOfCosts[t+1] = q_bestTau_t + Penalty;//m_t+\beta
      NbOfCandidats.push_back(changes_t.size());
      (candidates[t+1]->cost_tau_t)->initialize(VectOfCosts[t+1]); //add m_t+\beta in cost at time t+1
      
    //save pastSpheres for GeomFPOP(S-type) and GeomFPOP(R-type*) 
    if (type_algo >= 11 && type_algo <= 14) {
      for (std::list<unsigned int>::iterator iter = changes_t.begin(); iter != changes_t.end(); iter++) {
        costDif->initialize(candidates[*iter]->cost_tau_t, candidates[t+1]->cost_tau_t);
        pSphere<p> * pSph = new pSphere<p>();
        costDif->getSphere(pSph);
        (pastSets[t+1]->pastSpheres).push_back(pSph);
      }
    }
    
    //STEP 1: END//
    // STEP 2: PRUNING//
    //PELT
    if (type_algo == 1) {oneIterPruningPELT(VectOfCosts[t+1], t, changes_t, candidates);}
    //GEOMFPOP
    if (type_algo >= 2 && type_algo <=14) {
      if (type_algo == 11) {oneIterPruningStypeGeomFPOP(t, changes_t, candidates, pastSets);} 
      else if  ( type_algo < 11) {oneIterPruningGeomFPOP(nbInter, nbExcl,t, changes_t,candidates, blocks,nbRandIner,nbRandExcl);} 
      else {oneIterPruningGeomFPOPallPast(nbInter, t, changes_t, candidates, blocks, pastSets,nbRandIner);}
    }
      changes_t.push_back(t+1);//ADD NEW CHANGE
    //STEP 2: END//
    }
    //RETURN OPTIMAL SEGMENTATION//
    backtracking();
    return ResAlgo(showNbCands);
  }
  
  //PRUNING-------------------------------------------------------------------//
  //PELT----------------------------------------------------------------------//
  void oneIterPruningPELT(double boundary, unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates) {
    
    std::list<unsigned int>::iterator iter = changes_t.begin();
    for (unsigned int i = 0; i < NbOfCandidats[t]; i++) {
      if ( (candidates[*iter]->cost_tau_t)->getMin() >= boundary) {
        iter = changes_t.erase(iter);
      } else { 
        iter++;
      }
    }
  }
  //GEOMFPOP------------------------------------------------------------------//
  void RintersectionGeomFPOP (unsigned int k, bool &isPruning, int nbInter, unsigned int nbCds, unsigned int t, unsigned int * &orderedchanges, candidate<p> ** &candidates, pRectangle<p> ** &blocks, pCost<p> *&costDif, pSphere<p> * &sphere, int nbRandIner) {
    //last intersection
    costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[t+1]->cost_tau_t);
    costDif->getSphere(sphere);
    blocks[orderedchanges[k]]->IntersectionSphere(sphere);
    isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
    
    //others intersections (+random or all)
    if(!isPruning && (nbInter > 0) && (k < (nbCds-1))) {
      if (nbRandIner >= (nbCds - k - 1)) {
        nbInter = 10;
      }
      unsigned int j = k;
      unsigned int countRand = 0;
      while (!isPruning && (j < (nbCds-1)) && (countRand != nbRandIner)) {//without last sphere
        if (nbInter == 1) {
          j = k + (std::rand() % (nbCds - k - 1));//{k,..,t-1}
        }
        costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[orderedchanges[j+1]]->cost_tau_t);
        costDif->getSphere(sphere);
        blocks[orderedchanges[k]]->IntersectionSphere(sphere);
        isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
        if(nbInter == 1) {
          countRand++;
        } else {
          j++;
        }
      }
    }
  }
  
  //GEOMFPOP(R-type)----------------------------------------------------------//
  void oneIterPruningGeomFPOP(int nbInter, int nbExcl, unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates, pRectangle<p> ** &blocks, int nbRandIner, int nbRandExcl) {
    //INITIALIZATION
    unsigned int nbCds = NbOfCandidats[t];
    pCost<p> *costDif = new pCost<p>();
    pSphere<p> * sphere = new pSphere<p>();
    bool isPruning;
    std::list<unsigned int>::iterator iter = changes_t.begin();
    //copy potential candidates at time t
    unsigned int * orderedchanges = new unsigned int[nbCds]; 
    for (size_t i = 0; i < NbOfCandidats[t]; i++) {
      orderedchanges[i] = *iter;
      iter++;
    }
    //Update candidates
    for (size_t k = 0; k < nbCds; k++) {
      isPruning = false;//candidate is not empty
      RintersectionGeomFPOP (k, isPruning, nbInter, nbCds,  t,orderedchanges, candidates, blocks, costDif, sphere,nbRandIner);
      //exclusion(random or all)
      if(!isPruning && (nbExcl > 0) && (k != 0)) {
        unsigned int j = 0;
        if (nbRandExcl > k) {
          nbExcl = 10;
        }
        unsigned int countRand = 0;
        while (!isPruning && (j < k) && (countRand != nbRandExcl)) {
          if (nbExcl == 1) {
           // srand(time(0));
            j = std::rand() % k; //{0,...,k-1}
          }
          costDif->initialize( candidates[orderedchanges[j]]->cost_tau_t, candidates[orderedchanges[k]]->cost_tau_t);
          costDif->getSphere(sphere);
           blocks[orderedchanges[k]]->ExclusionSphere(sphere);
          isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
          if(nbExcl == 1) {countRand++;} 
          else {j++;}
        }
      }
      if (isPruning) {candidates[orderedchanges[k]]->isPruning = true;}
    }
    //candidate pruning
    pruning_t(changes_t,candidates, NbOfCandidats[t]);
    delete [] orderedchanges;
    orderedchanges = NULL;
  }
  
  //GEOMFPOP allPast spheres(R-type)----------------------------------------------------------//
  void oneIterPruningGeomFPOPallPast(int nbInter,unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates, pRectangle<p> ** &blocks, pastset<p>** &pastSets, int nbRandIner) {
    //INITIALIZATION
    unsigned int nbCds = NbOfCandidats[t];
    pCost<p> *costDif = new pCost<p>();
    pSphere<p> * sphere = new pSphere<p>();
    bool isPruning;
    std::list<unsigned int>::iterator iter = changes_t.begin();
    //copy potential candidates at time t
    unsigned int * orderedchanges = new unsigned int[nbCds]; 
    for (size_t i = 0; i < NbOfCandidats[t]; i++) {
      orderedchanges[i] = *iter;
      iter++;
    }
    //Update candidates
    for (size_t k = 0; k < nbCds; k++) {
      isPruning = false;//candidate is not empty
      RintersectionGeomFPOP (k, isPruning, nbInter, nbCds,  t,orderedchanges, candidates, blocks, costDif, sphere,nbRandIner);
      //exclusion(random or all)
      if (((pastSets[orderedchanges[k]]->pastSpheres).size() > 0) && !isPruning) {
        typename std::list<pSphere<p>*>::iterator iterP = (pastSets[orderedchanges[k]]->pastSpheres).begin();
        while ( iterP != (pastSets[orderedchanges[k]]->pastSpheres).end() && !isPruning) {
          if (blocks[orderedchanges[k]]->isOutside( *iterP)) {
            iterP = (pastSets[orderedchanges[k]]->pastSpheres).erase(iterP);
          }
         else {
            blocks[orderedchanges[k]]->ExclusionSphere (*iterP);
            isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
            ++iterP;
          }
        }
      }
      if (isPruning) { candidates[orderedchanges[k]]->isPruning = true;}
    }
    //candidate pruning
    pruning_t(changes_t,candidates, NbOfCandidats[t]);
    delete [] orderedchanges;
    orderedchanges = NULL;
  }
  
  
  //GEOMFPOP(S-type)----------------------------------------------------------//
  void oneIterPruningStypeGeomFPOP( unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates, pastset<p>** &pastSets) {
    unsigned int nbCds = NbOfCandidats[t];
    std::list<unsigned int>::iterator iter = changes_t.begin();
    bool isPruning;
    unsigned int * orderedchanges = new unsigned int[nbCds]; 
    for (size_t i = 0; i < NbOfCandidats[t]; i++) {
      orderedchanges[i] = *iter;
      iter++;
    }
    pSphere<p> * sphereL = new pSphere<p>();
    pSphere<p> * sphere = new pSphere<p>();
    pCost<p> *costDif = new pCost<p>();
    
    for (size_t k = 0; k < nbCds; k++) {
      isPruning = false;//candidate is not empty
      //last sphere
      costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[t+1]->cost_tau_t);
      costDif->getSphere(sphereL);
      //exclusion to past
      if (((pastSets[orderedchanges[k]]->pastSpheres).size() > 0) && !isPruning) {
        typename std::list<pSphere<p>*>::iterator iterP = (pastSets[orderedchanges[k]]->pastSpheres).begin();
        while ( iterP != (pastSets[orderedchanges[k]]->pastSpheres).end() && !isPruning) {
          isPruning = sphereL->isInclusion(*iterP);
          ++iterP;
        }
      }
      //intersection
      if(!isPruning && (k < (nbCds-1))) {
        unsigned int j = k;
        while (!isPruning && (j < (nbCds-1))) {//without last sphere
          costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[orderedchanges[j+1]]->cost_tau_t);
          costDif->getSphere(sphere);
          isPruning = sphereL->isnotIntersection(sphere);
          j++;
        }
      }
      if (isPruning) {candidates[orderedchanges[k]]->isPruning = true;}
    }
    //pruning
    pruning_t(changes_t,candidates, NbOfCandidats[t]);
    delete [] orderedchanges;
    orderedchanges = NULL;  
  }
};

#endif //ALGOS_H

 
