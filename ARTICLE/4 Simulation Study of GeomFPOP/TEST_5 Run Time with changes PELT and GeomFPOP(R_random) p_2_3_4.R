################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# TEST_5:                                                                     |
# The time complexity dependence on the number of segments.                   |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We consider ts with n =10^6 data points                                 |     
# number of segments = (1,2,5)x10^i,10^4 , i = 0,..,3                     |
# The mean for even segments was equal to 1, for odd segments - 0.        |
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |  
# nbSimus_ = 100 is the number of simulations.                            |   
# The test results are saved in files of the following type:              | 
# 'DependenceTC_Ints_p_1_GeomFPOP_R_random_random_it_100_.txt',           |
# with p = 2, 3, 4                                                        |
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#library(devtools)
#devtools::install_github("lpishchagina/GeomFPOP", force = TRUE)
library(GeomFPOP)
library(iterators)
library(stats)
library(parallel)

#function----------------------------------------------------------------|
#returns time complexity
RuntimeFpopOneSimuChangesk <- function(nb_ints, 
                                       cnst_, 
                                       n_, p_, 
                                       pentype_ = 1, 
                                       method_, 
                                       type_, 
                                       intersection_, 
                                       exclusion_, 
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
  res <- system.time(getChangePoints(ts_, 
                                     pen_, 
                                     method_, 
                                     type_, 
                                     intersection_, 
                                     exclusion_))[[1]]
  return (res);
}

#method parameters-------------------------------------------------------|
method_ = c('GeomFPOP', 'PELT')
type_ = 'R'
intersection_ = c('random','all')
exclusion_ = c('random','all*')

#parameters--------------------------------------------------------------|
n <- 10^6
p <- c(2, 3, 4)
#segment number
nb_ints <- c(1, 2, 5, 
             10, 20, 50, 
             100, 200, 500,
             1000, 2000, 5000, 
             10000)
#penalty
pentype_ <- 1 # 0 # if '0' then 'penalty = 2log(n)'; if '1' then 'penalty = 2dim*log(n)'
#simulation number
nbSimus_ <- 100

#calculations------------------------------------------------------------|
for (t_dim in 1:length(p)) {
  for (t_a in 1:length(method_)) {
    test_results <- matrix(NA, length(nb_ints),  nbSimus_)
    index <- 1
    ints <- NULL
    while (index <= length(nb_ints)) {
      ints <- nb_ints [index]
      test_results[index , ]  <- do.call(cbind, parallel::mclapply(
        1:nbSimus_,
        function(i) RuntimeFpopOneSimuChangesk(nb_ints = ints, 
                                               cnst_ = i, 
                                               n_ = n, 
                                               p_ = p[t_dim], 
                                               pentype_ = 1, 
                                               method_[t_a], 
                                               type_, 
                                               intersection_[t_a], 
                                               exclusion_[t_a], 
                                               func = "getChangePoints", 
                                               noise = 1), mc.cores = nbSimus_
      ))
      #file1 <- paste('DependenceTC_Ints', ints, p[t_dim], pentype_, method_[t_a], type_, intersection_[t_a], exclusion_[t_a], 'it', nbSimus_,'.txt', sep = '_')
      #write.table(data.frame (ints = ints, mean(test_results[index,])), file1, row.names = TRUE, col.names = FALSE)
      index <- index + 1
    }
    #save result--------------------------------------------------------------------
    file <- paste('DependenceTC_Ints', 
                        p[t_dim], 
                        pentype_, 
                        method_[t_a], 
                        type_, 
                        intersection_[t_a], 
                        exclusion_[t_a], 
                        'it', 
                        nbSimus_,
                        '.txt', 
                        sep = '_')
    write.table(data.frame (ints = nb_ints, rowMeans(test_results, na.rm = TRUE)), file, row.names = TRUE, col.names = FALSE)
    
  }
}
################################################################################
########################### END ################################################
################################################################################


