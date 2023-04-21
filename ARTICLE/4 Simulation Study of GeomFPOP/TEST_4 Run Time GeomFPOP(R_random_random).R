
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# TEST 4:                                                                     |
# Time Complexity of GeomFPOP(R-type: random/random)                          |
#for p = 2,..,10,100 (data without changes)                                   |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We consider ts ~ N_p(0,1) with n = 2^10:2^25 (without changes)          |                                  
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |  
# nbSimus_ = 100 is the number of simulations.                            |   
# The test results are saved in files of the following type:              | 
# 'TC_Res_p_1_GeomFPOP_R_random_random_it_100_.txt',                      |
# with p = 2,..,10,100 and Time limit:TimeLimit = 200 s                   |
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#library(devtools)
#devtools::install_github("lpishchagina/GeomFPOP", force = TRUE)
library(GeomFPOP)
library(iterators)
library(stats)
library(parallel)

#function----------------------------------------------------------------|
#returns time complexity for ts ~ N_p(0,1) (without changes) 
RuntimeFpopOneSimu <- function(cnst_, 
                               n_, 
                               p_, 
                               pentype_ = 1, 
                               method_, 
                               type_, 
                               intersection_, 
                               exclusion_, 
                               func = "getChangePoints", 
                               noise = 1) {
  pen_ <- ifelse( pentype_ == 0, 2*log(n_), 2*p_*log(n_))
  set.seed(cnst_ + 10)
  ts_ <- rnormChanges(p_, 
                      n_, 
                      changes = NULL, 
                      means = matrix(0, ncol = 1, nrow = p_), 
                      noise)
  res <- system.time(getChangePoints(ts_, 
                                     pen_, 
                                     method_, 
                                     type_, 
                                     intersection_, 
                                     exclusion_))[[1]]
  return (res);
}

#method parameters--------------------------------------------------------|
method_ = 'GeomFPOP'
type_ = R
intersection_ = 'random'
exclusion_ = 'random'

#parameters---------------------------------------------------------------|
degree <- 10:25                                                              
n_ <-2^(degree)
dim_ <- c(2, 3, 4, 5, 6, 7, 8, 9, 10,100)
#By default, penalty = 2*p*log(n)
pentype_ = 1
#number of simulations
nbSimus_ <- 100
#Time limit
TimeLimit <- 200

#calculations-------------------------------------------------------------|
for (p in 1 : length(dim_)) {
  simulation_res <- matrix(NA, length(n_),  nbSimus_)
  index <- 1
  Time <- 0
  while ((index <= length(n_)) && (Time <= TimeLimit)) {
    N <- n_[index]
    simulation_res[index , ]  <- do.call(cbind, parallel::mclapply(
      1:nbSimus_, function(i) RuntimeFpopOneSimu(cnst_ =i, 
                                                 n_ = N, 
                                                 p_ = dim_[p], 
                                                 pentype_ = 1, 
                                                 method_, 
                                                 type_, 
                                                 intersection_, 
                                                 exclusion_, 
                                                 func = "getChangePoints", 
                                                 noise = 1), mc.cores = nbSimus_))
    Time <- max(simulation_res[index , ])
    index <- index + 1
  }
  #save result : average runtime
  file <- paste('TC', 
                      dim_[p], 
                      '1_GeomFPOP_R_random_random_it', 
                      nbSimus_,
                      '.txt', 
                      sep = '_')
  write.table(data.frame (n = n_, rowMeans(simulation_res, na.rm = TRUE)), file, row.names = TRUE, col.names = FALSE)
}
################################################################################
########################### END ################################################
################################################################################


