
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.1. The Number of Potential Change-Point Candidates Stored Over Time      ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 5:                                                                   |
# Percentage of candidate change-points stored over time by GeomFPOP          |
#  with R (left) or S (right) type pruning for dimension 1 < p < 11           |
# We simulated 100 i.i.d Gaussian data N_p (0, I) and report the average.     |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_1                                   |
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(base)
library(rstream)
library(tidyverse)
library(stats)
library(RColorBrewer)
library(ggpubr)

#parameters (by default)--------------------------------------------------|
n_ <- 10^4
dim_ <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
namedim_ <- c ("p = 2",
               "p = 3", 
               "p = 4", 
               "p = 5", 
               "p = 6", 
               "p = 7", 
               "p = 8", 
               "p = 9", 
               "p = 10")
algo_ <- c("GeomFPOP (S-type)", 
           "GeomFPOP (R-type)")

nbSimus_ <- 100
mtd_ <- c("GeomFPOP_S_all_all*", 
          "GeomFPOP_R_all_all*")
pentype_ <- 1
nbDim_ <- length(dim_)

#read results-------------------------------------------------------------|
read_file <- NULL
test_results <- list()
for (j in  1 : length(algo_)) {
  test_results[[j]] <- data.frame(time = c(1 : n_), matrix(0, n_, nbDim_))
  colnames(test_results [[j]]) <- c('Time', namedim_)
  for (k  in 1 : nbDim_) {
    read_file <- paste('NbCds_Res', 
                       n_, 
                       dim_[k], 
                       pentype_, 
                       mtd_[j],
                       'it',
                       nbSimus_, 
                       '.txt', 
                       sep = '_')
    test_results[[j]] [, k+1] <- read.table(file = read_file, row.names = 1)
  }
}
#Ratio--------------------------------------------------------------------|
Ratio <- list()
for (t in 1:2) {
  Ratio[[t]] <- test_results [[t]]
  Ratio[[t]][, -1] <- test_results[[t]][, -1]/test_results[[t]][, 1]
}

# modification of data frame----------------------------------------------|
res <- do.call(rbind, lapply(
  1:length(Ratio), function(i) { 
    data.frame(Time = Ratio[[i]]$Time,
               y = c(Ratio[[i]][,2], 
                     Ratio[[i]][,3], 
                     Ratio[[i]][,4], 
                     Ratio[[i]][,5], 
                     Ratio[[i]][,6], 
                     Ratio[[i]][,7], 
                     Ratio[[i]][,8], 
                     Ratio[[i]][,9], 
                     Ratio[[i]][,10]),
               dimension = rep(colnames(Ratio[[i]][-1]), each = nrow(Ratio[[i]])),
               talgo_ = algo_[i])
  }))

#plot--------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = Time, y = y, color = dimension)) +
  facet_grid(cols = vars(talgo_), scales = "fixed") +
  geom_line() +
  scale_x_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  scale_y_continuous(trans = "log10",
                     labels = scales::math_format(10^.x, format = log10)) +
  ylab("Ratio value") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')

#save---------------------------------------------------------------------|
pdf(file = "Figure5.pdf",  width = 8, height = 4)
print(Plot)
dev.off()
################################################################################
########################### END ################################################
################################################################################

