#include <iostream>
#include "math.h"
#include <Rcpp.h>
#include "Algos.h"

using namespace Rcpp;
using namespace std;
/* OP: method = "OP" 
 * PELT: method = "PELT" 
 * GeomFPOP(R-type)
 * method = "GeomFPOP"
 * type = "R"
 * type = "S"
 * intersection = "all" intersection = "last" intersection = "random"
 *  exclusion = "all" exclusion = "all*" exclusion = "random" exclusion = "empty" 
 *  nbRandInter, nbRandExcl - min number of intersections an exclusions, respectively
 */

//converting parameters("method","type", "intersection", "exclusion") to a numeric value.
unsigned int typeAlgo(std::string method, std::string type, std::string intersection,  std::string exclusion) {
  unsigned int type_algo = INFINITY;
  if (method == "OP") { type_algo = 0; }
  if (method == "PELT") { type_algo = 1; }
  if (method == "GeomFPOP") { 
    if (type == "R") {
      if (exclusion == "all") { 
        if (intersection == "all") { type_algo = 2;}
        if (intersection == "last") {type_algo = 3;}
        if (intersection == "random") {type_algo = 4;}
      }
      if (exclusion == "random") {
        if (intersection == "all") {type_algo = 5;}
        if (intersection == "random") {type_algo = 6;}
        if (intersection == "last") {type_algo = 7;}
      }
      if (exclusion == "empty") {
        if (intersection == "all") {type_algo = 8;}
        if (intersection == "last") {type_algo = 9;}
        if (intersection == "random") {type_algo = 10;}
      }
      if (exclusion == "all*") {
        if (intersection == "all") {type_algo = 12;}
        if (intersection == "last") {type_algo = 13;}
        if (intersection == "random") {type_algo = 14;}
      }
      
    }
    if (type == "S") { type_algo = 11;}
  }

  return type_algo;
}

//' @title getChangePoints
//'
//' @description Multiple Changepoint Detection using methods: OP, PELT, GeomFPOP.
//' @param data is a matrix of data (p-rows x n-columns, where p - dimension, n - length).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param method is the algorithm: 'OP', 'PELT', 'GeomFPOP'.
//' @param type is the approximation type for GeomFPOP: 'R' -rectangle or 'S'- last sphere ( by default, 'R').
//' @param intersection is the type of intersection for GeomFPOP(method = 'R'):  'all', 'last' or 'random' (by default, 'last').
//' @param exclusion is the type of exclusion for GeomFPOP(method = 'R'): 'empty', 'all', 'all*' or 'random'(by default, 'all').
//' @param showNbCands is the logical parameter (if "true", than to show the number of candidates at each iteration).
//' @param nbRandInter is the minimum number of intersections at each iteration for 'intersection = random' (by default, nbRandInter = 1)
//' @param nbRandInter is the minimum number of exclusions at each iteration for 'exclusion = random'(by default, nbRandExcl = 1)
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
//' }
//'
//' @examples
//' N <- 100000
//' Chpt <-5000
//' Means <-  matrix(c(0,1,1,10), nrow = 2) 
//' Noise <- 1
//' Dim <- 2
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' Dim <- 3
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = 3, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 3), noise = 1)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'OP',showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'PELT', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'S', showNbCands = FALSE)#GeomFPOP(S-type)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'last', exclusion = 'all', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'last', exclusion = 'all*', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'last', exclusion = 'random', nbRandExcl = 1, showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'last', exclusion = 'empty', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'all', exclusion = 'all', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'all', exclusion = 'all*', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'all', exclusion = 'random', nbRandExcl = 1, showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'all', exclusion = 'empty', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'random',nbRandInter = 1, exclusion = 'all', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'random',nbRandInter = 1, exclusion = 'all*', showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'random', nbRandInter = 1, exclusion = 'random', nbRandExcl = 1, showNbCands = FALSE)
//' getChangePoints(data = time_series, penalty = Penalty, method = 'GeomFPOP', type = 'R', intersection = 'random', nbRandInter = 1, exclusion = 'empty', showNbCands = FALSE)

// [[Rcpp::export]]
List getChangePoints(Rcpp::NumericMatrix data, double penalty, std::string method = "GeomFPOP", std::string type = "R", std::string intersection = "all",  std::string exclusion = "all", bool showNbCands = false, int nbRandInter = 1, int nbRandExcl = 1) {
  unsigned int type_algo = typeAlgo(method, type, intersection, exclusion);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_algo == INFINITY){throw std::range_error("This combination of parameters is not available.");}
//  if( ((unsigned int)data.nrow() < 1) ||  ((unsigned int)data.nrow() > 20)) {throw std::range_error("The dimension of time series can not exceed 20.");}
  //----------------------------------------------------------------------------
  unsigned int p = (unsigned int)data.nrow();
  if (p == 1){
    Algos<1> X = Algos<1>(data, penalty);
    return X.algosOP(type_algo, showNbCands,nbRandInter, nbRandExcl);
  } else 
  if (p == 2){
    Algos<2> X = Algos<2>(data, penalty);
    return X.algosOP(type_algo, showNbCands,nbRandInter, nbRandExcl);
  } else if (p == 3){
    Algos<3> X = Algos<3>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 4){
    Algos<4> X = Algos<4>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 5){
    Algos<5> X = Algos<5>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 6){
    Algos<6> X = Algos<6>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 7){
    Algos<7> X = Algos<7>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 8){
    Algos<8> X = Algos<8>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 9){
    Algos<9> X = Algos<9>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 10){
    Algos<10> X = Algos<10>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 11){
    Algos<11> X = Algos<11>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 12){
    Algos<12> X = Algos<12>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 13){
    Algos<13> X = Algos<13>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 14){
    Algos<14> X = Algos<14>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 15){
    Algos<15> X = Algos<15>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 16){
    Algos<16> X = Algos<16>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 17){
    Algos<17> X = Algos<17>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 18){
    Algos<18> X = Algos<18>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 19){
    Algos<19> X = Algos<19>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  } else if (p == 20){
    Algos<20> X = Algos<20>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  }
//p=100
  else if (p == 100){
    Algos<100> X = Algos<100>(data, penalty);
    return X.algosOP(type_algo, showNbCands, nbRandInter, nbRandExcl);
  }
  //
  return NULL;
}
