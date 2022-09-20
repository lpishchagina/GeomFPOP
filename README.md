<a id="top"></a>
#  GeomFPOP Vignette
### Liudmila Pishchagina
### September 20, 2022

## Quick Start

` GeomFPOP ` is an R package written in Rcpp/C++ and developed to detect changes in using the methods: Optimal Partitioning (OP), Pruned Exact Linear Time (PELT) and Geometric Functional Pruning Optimal Partitioning (GeomFPOP) in `p`-variate time series (`p`-variate Gaussian Model) of length `n`, where  `p` in `{2,.., 20}`. 


## The package installation

We install the package from Github:

```r
#devtools::install_github("lpishchagina/GeomFPOP")
library(GeomFPOP)
```

## Generation of time series (`p`-variate Gaussian Model)

We use the following function for data generation:

###function chpt_rnorm

The `rnormChanges` is the generation of data (`p`-variate Gaussian Model) of dimension `p` with a given values of means and changes.

`n`  is the time series length.

`p`  is the time series dimension.

`changes` is the changepoint vector that gives the last index of each segment.

The last element of `changes` is always less than to the length of time series.

By default, `changes = NULL` (for the data without changes). 

`means` is the matrix of successive means for the `p`-variate time series.

By default, `means = matrix(0, ncol = 1, nrow = p)` (for the data without changes). 

The length of each matrix row is equal to the length of `changes` plus one.

`noise` is a variance of the time series. By default, `noise = 1`.

```r
set.seed(21)
N <- 1000
Chpt <-500
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
Penality <- 2*Dim*log(N)

##the time series with one change
tsOneChange <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)

tsNoChange <- rnormChanges(p = Dim, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = Dim), noise = Noise)
```


## Multiple Changepoint Detection

We consider the multiple change-point detection problem in large data. Dynamic programming algorithms exist to solve this problem exactly. Currently, there is an increasing need for methods that can not only accurately, but also fast  detect change-points in `p`-variate large time-series. This package implements three methods of dynamic programming: two basic methods, Optimal Partitioning (OP) and Pruned Exact Linear Time (PELT),  and Geometric Functional Pruning Optimal Partitioning (GeomFPOP).

OP solves our problem at quadratic time.  PELT uses an inequality test, inequality-based pruning, to discard some indices in the minimization operator. Its computational cost, on average, is linear in the number of data points `n` when the number of changes is linear in `n` yet quadratic if the number of changes is smaller. It is a problem for many real-world applications. To get around this problem for uni-variate time series, another method, called Functional Pruning Optimal Partitioning (FPOP), has recently been developed. In practice, even in the absence of changes, the obtained time complexity is `O(nlog(n))`. FPOP is fast even in the absence of true changes but it only works for uni-dimensional data. 

The extension of FPOP to the multivariate setting is not straightforward as one needs to follow non-convex sets obtained as the intersection and difference of sets in `R^p`. In this package we propose an extension of FPOP to the `p`-dimensional setting (`p`-variate Gaussian Model), GeomFPOP. Our idea is to getChangesimate those sets using simpler geometric shapes that are easy to update (hyper-spheres and hyper-rectangles,`type = S` and `type = R`, respectively).
We consider GeomFPOP with different number of intersections and exclusions at each iterations. In this package we implemented the following modifications:

`intersection = all` - we use for intersection  all possible sets.

`intersection = last` - we use for intersection  only last set.

`intersection = random` - we use for intersection   last set and one another random set from the possible sets.

`exclusion = empty` - we do not make exclusions.

`exclusion = random` - we make exclusion with one random set from the possible sets.

`intersection = all` - we use for exclusions  all possible sets (their changes exist at time 't') .

`intersection = all*` - we use for exclusions  all possible sets.

Empirical estimation of pruning and run-time efficiency of GeomFPOP demonstrates its fast work for long `2`-dimensional and `3`-dimensional signals.

### The function getChangePoints

The ` getChangePoints ` function returns the result of the segmentation.

` data ` is the `p`-variate time series (matrix of real numbers with `p`-rows (dimension) and `n`-columns(length)).

` penalty ` is the value of penalty (a non-negative real number  equals to a classic `2*p*(noise^2)*log(n)`. 

` method ` is the type of algorithm : `'OP'`, `'PELT'`, `'GeomFPOP'`  (by default, `'GeomFPOP'`).

` type ` is the getChangesimation type for GeomFPOP : `'R'`, `'S'` (by default, `'R'`).

` intersection ` is the type of intersection : `'all'`, `'last'` or  `'random'`.

` exclusion ` is the type of exclusion : `'empty'`, `'all'` , `'all*'` or `'random' `.

` showNbCands `is the logical parameter (` showNbCands = TRUE ` - to print the number of change-point candidates for each iteration).

```r


getChanges <- list()

getChanges[[1]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'OP', showNbCands = FALSE)

getChanges[[2]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'PELT', showNbCands = FALSE)

getChanges[[3]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'S', showNbCands = FALSE)

getChanges[[4]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'all', exclusion = 'all', showNbCands = FALSE)

getChanges[[5]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'last', exclusion = 'all', showNbCands = FALSE)

getChanges[[6]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'random', exclusion = 'all', showNbCands = FALSE)


getChanges[[7]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'all', exclusion = 'all*', showNbCands = FALSE)

getChanges[8] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'last', exclusion = 'all*', showNbCands = FALSE)

getChanges[[9]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'random', exclusion = 'all*', showNbCands = FALSE)


getChanges[[10]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'all', exclusion = 'random', showNbCands = FALSE)

getChanges[11] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'last', exclusion = 'random', showNbCands = FALSE)

getChanges[[12]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'random', exclusion = 'random', showNbCands = FALSE)


getChanges[[13]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'all', exclusion = 'empty', showNbCands = FALSE)

getChanges[14] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'last', exclusion = 'empty', showNbCands = FALSE)

getChanges[[15]] <- getChangePoints(data = tsNoChange, penalty = Penality, method = 'GeomFPOP',  type = 'R',  intersection = 'random', exclusion = 'empty', showNbCands = FALSE)

```

```r
#OP
getChanges[[1]] 
#PELT
getChanges[[2]] 
#GeomFPO(S-type)
getChanges[[3]] 
#GeomFPO(R-type, all/all)
getChanges[[4]] 

#$changes
#numeric(0)

#$means
#$means[[1]]
#[1]  0.03372487 -0.01246343


#$UnpenalizedCost
#[1] 1975.554

```

`changes` is the change-point vector that gives the last index of each segment.

The last element of `changes` always equals to the length of time series.

`means` is the list of successive means for the p-variate time series.

`globalCost` is the overall Gaussian cost of the segmented data. 

[Back to Top](#top)
