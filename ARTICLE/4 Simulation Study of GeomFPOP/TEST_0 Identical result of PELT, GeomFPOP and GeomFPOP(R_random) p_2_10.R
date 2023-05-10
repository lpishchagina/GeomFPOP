################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
################################################################################

#-----------------------------------------------------------------------------|
# TEST_0:                                                                     |
# The results of PELT, GeomFPOP and GeomFPOP(r/r).                            |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We consider ts with n =10^4 data points                                 |     
# number of segments = c(100, 50, 10, 5, 1)                               |
# The mean for even segments was equal to 1, for odd segments - 0.        |
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |  
# nbSimus_ = 100 is the number of simulations.                            |   
# The test results are saved in files of the following type:              | 
# 'ResIsIdentical_N_Ints_p_1_it_100_.txt',                                |
# with p = 2,.., 10                                                       |
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#install.packages("devtools")
#library(devtools)
#devtools::install_github("lpishchagina/GeomFPOP", force = TRUE)
library(GeomFPOP)
#library(iterators)
library(stats)
library(parallel)

#function----------------------------------------------------------------|
#function returns true if the results are identical-----------------------------
res_One_simulation <- function(nb_ints, 
                               cnst_, 
                               n_, p_, 
                               pentype_ = 1, 
                               func = "getChangePoints", 
                               noise = 1) {
  pen_ <- ifelse( pentype_ == 0, 2*log(n_), 2*p_*log(n_))
  set.seed(cnst_ + 10)
  deltaInt <- n_/nb_ints
  changes_<-NULL
  means_ <-matrix(0, nrow = p_, ncol = nb_ints,byrow = TRUE)
  if(nb_ints != 1){
    for (i in 1:(nb_ints-1)) {
      changes_ <- c(changes_, i*deltaInt)
    }
    i=1:(nb_ints/2)
    means_[,2*i] <- 1 
  }  
  ts_ <- rnormChanges(p_, n_, changes_, means_, noise)
  resPELT <- getChangePoints(ts_, pen_, 'PELT', 'R', 'all', 'all*')$UnpenalizedCost
  resGeomFPOP <- getChangePoints(ts_, pen_, 'GeomFPOP', 'R', 'all', 'all*')$UnpenalizedCost
  resGeomFPOP_RR <- getChangePoints(ts_, pen_, 'GeomFPOP', 'R', 'random', 'random')$UnpenalizedCost
  res <- (resPELT == resGeomFPOP)&&(resGeomFPOP == resGeomFPOP_RR)
  return (res);
}
#parameters--------------------------------------------------------------|
n <- 10^4
p <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
#segment number
nb_ints <- c(100, 50, 10, 5, 1)
#penalty
pentype_ <- 1 # 0 # if '0' then 'penalty = 2log(n)'; if '1' then 'penalty = 2dim*log(n)'
#simulation number
nbSimus_ <- 100

#calculations------------------------------------------------------------|
for (t_dim in 1:length(p)) {
  test_results <- matrix(NA, length(nb_ints),  nbSimus_)
  index <- 1
  ints <- NULL
  while (index <= length(nb_ints)) {
    ints <- nb_ints [index]
    test_results[index , ]  <- do.call(cbind, parallel::mclapply(
      1:nbSimus_,
      function(i) res_One_simulation(nb_ints = ints, 
                                             cnst_ = i, 
                                             n_ = n, 
                                             p_ = p[t_dim], 
                                             pentype_ = 1, 
                                             func = "getChangePoints", 
                                             noise = 1), mc.cores = nbSimus_
    ))
    index <- index + 1
  }
  #save result--------------------------------------------------------------------
  res_id <- (rowSums(test_results, na.rm = TRUE) == rep(nbSimus_, length(nb_ints)))
  file <- paste('ResIsIdentical_',n, 
                p[t_dim], 
                pentype_, 
                'it', 
                nbSimus_,
                '.txt', 
                sep = '_')
  write.table(data.frame (ints = nb_ints, identify = res_id ), file, row.names = TRUE, col.names = FALSE)
}
################################################################################
########################### END ################################################
################################################################################
p

