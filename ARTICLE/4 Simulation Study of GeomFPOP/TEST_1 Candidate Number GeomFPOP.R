################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.1. The Number of Potential Change-Point Candidates Stored Over Time      ##
################################################################################

################################################################################
#-----------------------------------------------------------------------------|
# TEST_1:                                                                     |
# The number of candidates of change stored over time by GeomFPOP (R and S)   |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We consider ts ~ N_p(0,1) with n = 10^4 (without changes),  1 < p < 11  |                                  
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |   
# nbSimus_ = 100 is the number of simulations.                            |   
# The test results are saved in files of the following type:              | 
# 'NbCds_Res_10000_p_1_GeomFPOP_type_all_all*_it_100_.txt',             |
# with 1 < p < 11  and type = R or S                                      | 
# penalty = 2*p*log(n)                                                    |
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
#library(devtools)
#devtools::install_github("lpishchagina/GeomFPOP", force = TRUE)
library(GeomFPOP)
library(iterators)
library(stats)
library(parallel)

#function-----------------------------------------------------------------|
#returns number of candidates for ts ~ N_p(0,1) (without changes) 
CandsNbOneSimu <- function(cnst_,
                           n_, 
                           dim_, 
                           pentype_ = 1, 
                           method_ ,
                           type_, 
                           intersection_, 
                           exclusion_,
                           noise_ = 1) {
  pen_  <- ifelse(pentype_ == 0, 2*log(n_), 2*dim_*log(n_))
  set.seed(cnst_ + 10)
  
  ts_ <- rnormChanges(dim_, 
                      n_, 
                      changes = NULL, 
                      means = matrix(0, ncol = 1, nrow = dim_), 
                      1)
  
  res <-getChangePoints(ts_, 
                        pen_, 
                        method_, 
                        type_, 
                        intersection_, 
                        exclusion_, 
                        showNbCands = TRUE)$NumberOfCandidats
  return (res);
}

#method parameters--------------------------------------------------------|
method_ = 'GeomFPOP'
type_ = c('R','S')
intersection_ = 'all'
exclusion_ = 'all*'

#parameters---------------------------------------------------------------|
#number of data points
n_ <-10^4
#dimensions
dim_ <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
#type of penalty
#By default, penalty = 2*p*log(n)
pentype_ = 1
#number of simulations
nbSimus_ <- 100

#calculations-------------------------------------------------------------|
for (p in 1 : length(dim_)) {
  for ( t in 1 : length(type_)) {
    simulation_res <- matrix(0, n_, nbSimus_)
    simulation_res <- do.call(cbind, parallel::mclapply(
      1:nbSimus_,
      function(i) CandsNbOneSimu(cnst_ = i,
                                 n_, 
                                 dim_[p], 
                                 pentype_, 
                                 method_, 
                                 type_[t], 
                                 intersection_, 
                                 exclusion_, 
                                  1),
              mc.cores = nbSimus_
    ))
    #save result: average number of candidates
    txtFileRes <- paste('NbCds_Res', 
                        n_, 
                        dim_[p], 
                        '1_GeomFPOP',
                        type_[t], 
                        'all_all*_it',
                        nbSimus_, '.txt', 
                        sep = '_')
    write.table(data.frame (time = c(1 : n_), 
                rowMeans(simulation_res)), 
                txtFileRes, 
                row.names = FALSE, col.names = FALSE)
  }
}
################################################################################
########################### END ################################################
################################################################################

dim_

