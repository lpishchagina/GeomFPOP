
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 7:                                                                   |
# Run-time of the (random/random) approache of GeomFPOP (R-type) and PELT     |
# using multivariate time-series without change-points.                       |
# The maximum run-time of the algorithms is 3 minutes.                        |        
# Averaged over 100 data sets, p = 2,..,10, 100                               |
#-----------------------------------------------------------------------------|
# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_3 and TEST_4                        |                                  
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(base)
library(rstream)
library(tidyverse)
library(stats)
library(RColorBrewer)
library(ggpubr)

#parameters---------------------------------------------------------------|
dim_ <- c (2, 3, 4, 5, 6, 7, 8, 9, 10, 100)
Length <-2^(10:25)
namedim_ <- c ("p = 2",
               "p = 3", 
               "p = 4", 
               "p = 5", 
               "p = 6", 
               "p = 7", 
               "p = 8", 
               "p = 9", 
               "p = 10", 
               "p = 100")
typeA <-c("GeomFPOP (R-type: random / random)", "PELT")

#read data frame----------------------------------------------------------|
nbSimus_ <- 100
typeAlgo <-c("GeomFPOP_R_random_random_it", "PELT_R_all_all*_it")
pentype_ <- 1
NbDims <- length(dim_)
read_file <- NULL
test_results <- list()
for (j in  1 : length(typeA)) {
  test_results[[j]] <- data.frame(n = Length, matrix(0, length(Length), NbDims))
  colnames(test_results [[j]]) <- c('n', namedim_)
  for (k  in 1 :NbDims) {
    read_file <- paste ('TC', 
                        dim_[k], 
                        pentype_,
                        typeAlgo[j], 
                        NbSimus,
                        '.txt', 
                        sep='_')
    test_results[[j]] [, k+1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}

#modification of data frame-----------------------------------------------|
res <- do.call(rbind, lapply(
  1:length(test_results), function(i) { 
    data.frame(n = test_results[[i]]$n,
               y = c(test_results[[i]][, 2], 
                     test_results[[i]][, 3], 
                     test_results[[i]][, 4], 
                     test_results[[i]][, 5], 
                     test_results[[i]][, 6], 
                     test_results[[i]][, 7], 
                     test_results[[i]][, 8], 
                     test_results[[i]][, 9], 
                     test_results[[i]][, 10], 
                     test_results[[i]][, 11]),
               dimension = rep(colnames(test_results[[i]][-1]), each = nrow(test_results[[i]])),
               Method = typeA[i]
    )
  }))

#plot---------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = n, y = y, color = Method)) +
  facet_wrap(~dimension, scales = "fixed", ncol = 5) +
  geom_line() + 
  ylab("Seconds") +
  xlab("Number of data points of time series") +
  geom_hline(yintercept = 180, linetype = "dashed") + 
  scale_x_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  theme_bw() +
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.text=element_text(size = 14), 
        legend.position = "bottom") +
  scale_color_brewer(palette = 'Set1')

#save---------------------------------------------------------------------|
pdf(file = "Figure7.pdf",  width = 8, height = 4)
print(Plot) 
dev.off()
################################################################################
########################### END ################################################
################################################################################
