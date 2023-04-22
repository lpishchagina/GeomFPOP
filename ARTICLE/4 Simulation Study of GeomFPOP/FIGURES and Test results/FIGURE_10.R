
################################################################################
##  APPENDIX E. Optimization Strategies for GeomFPOP(R-type)                  ##
##  The Number of Potential Change-Point Candidates Stored Over Time          ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 10:                                                                  |
# Percentage of candidate change-points stored over time                      |
# by  different optimization approaches of GeomFPOP                           |
# p = 2, 3, 4                                                                 |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_1 and TEST_6                        |                                  
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(base)
library(rstream)
library(tidyverse)
library(iterators)
library(stats)
library(RColorBrewer)
library(ggpubr)

#parameters (by default)--------------------------------------------------|
dim_ <- c (2, 3, 4)
nbSimus_ <- 100
n_ <- 10^4
pentype_ <- 1
namedim_ <- c ("p = 2",
               "p = 3", 
               "p = 4")

algo_ <- c("(all / all)", 
           "(all / empty)",
           "(all / random)",
           "(last / all)",
           "(last / empty)",
           "(last / random)",
           "(random / all)",
           "(random / empty)",
           "(random / random)")

mtd_ <- c("GeomFPOP_R_all_all*",
          "GeomFPOP_R_all_empty",
          "GeomFPOP_R_all_random",
          "GeomFPOP_R_last_all*",
          "GeomFPOP_R_last_empty",
          "GeomFPOP_R_last_random",
          "GeomFPOP_R_random_all*",
          "GeomFPOP_R_random_empty",
          "GeomFPOP_R_random_random")

#read results-------------------------------------------------------------|
nbMethods_ <- length(algo_)
readfile_ <- NULL
test_results <- list()
for (j in  1 : length(dim_)) {
  test_results[[j]] <- data.frame(time = c(1 : n_), matrix(0, n_, nbMethods_))
  colnames(test_results[[j]]) <- c('Time', algo_)
  for (k  in 1 : nbMethods_) {
    readfile_ <- paste('NbCds_Res', 
                       n_, dim_[j], 
                       pentype_, 
                       mtd_[k],
                       'it',
                       nbSimus_, 
                       '.txt', 
                       sep = '_')
    test_results[[j]] [, k+1] <- read.table(file = readfile_, row.names = 1)
  }
}

#Ratio--------------------------------------------------------------------|
Ratio <- list()
for (t in 1:3) {
  Ratio[[t]] <- test_results [[t]]
  Ratio[[t]][, -1] <- test_results [[t]][, -1]/test_results [[t]][, 1]
}

# modification of data frame----------------------------------------------|
res <- do.call(rbind, lapply(
  1:length(Ratio), function(i) { 
    data.frame(Time = Ratio[[i]]$Time,
               y = c(Ratio[[i]][, 2],
                     Ratio[[i]][, 3], 
                     Ratio[[i]][, 4], 
                     Ratio[[i]][, 5], 
                     Ratio[[i]][, 6], 
                     Ratio[[i]][, 7], 
                     Ratio[[i]][, 8], 
                     Ratio[[i]][, 9], 
                     Ratio[[i]][, 10]),
               Approach = rep(colnames(Ratio[[i]][-1]), each = nrow(Ratio[[i]])),
               dimension = namedim_[i])
  }))

#plot--------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = Time, y = y, color = Approach)) +
  facet_grid(cols = vars(dimension), scales = "fixed") +
  geom_line() +
  scale_x_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  scale_y_continuous(trans = "log10",
                     labels = scales::math_format(10^.x, format = log10)) +
  ylab("Ratio value") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom") +
  scale_color_brewer(palette = 'Paired')

#save---------------------------------------------------------------------|
pdf(file = "Figure10.pdf",  width = 10, height = 4)
print(Plot)
dev.off()
################################################################################
########################### END ################################################
################################################################################
