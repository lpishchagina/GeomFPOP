
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 9:                                                                   |
# Time complexity dependence of (random/random) approach of GeomFPOP (R-type) |
# on the number of segments in time series with 106 data points.              |
#-----------------------------------------------------------------------------|
# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_5                                   |                                  
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(base)
library(rstream)
library(tidyverse)
library(iterators)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)
#parameters---------------------------------------------------------------|
Dim <- c (2, 3, 4)
nameDim <- c("p = 2", "p = 3", "p = 4")
nameMethod <- c("GeomFPOP (\U0052-type : random / random)", "PELT")
Algo <- c("GeomFPOP_R_random_random","PELT_R_all_all*")
TP <- 1 
nbChanges <- c(1, 2, 5, 
               10, 20, 50, 
               100, 200, 500, 
               1000, 2000, 5000, 
               10000)
NbSimus <- 100
#read data frame----------------------------------------------------------|
NbMethods <- length(Algo)
file <- NULL
test_results <- list()
for (j in  1 : length(Dim)) {
  test_results[[j]] <- data.frame(nb_changes = nbChanges, 
                                     matrix(0, length(nbChanges), 
                                     NbMethods))
  colnames(test_results [[j]]) <- c('nb_changes', nameMethod)
  for (k  in 1 :NbMethods) {
    file <- paste ('DependenceTC_Ints', 
                               Dim[j], 
                               TP, 
                               Algo[k],
                               "it", 
                               NbSimus,
                               '.txt', 
                               sep='_')
    test_results[[j]] [, k+1] <- read.table(file = file, row.names = 1)[,2]
  }
}

#transformation of data frame---------------------------------------------|
res <- do.call(rbind, lapply(1:length(test_results), function(i) {
  data.frame(nb_changes = test_results[[i]]$nb_changes,
             y = c(test_results[[i]][, 2], 
                   test_results[[i]][, 3]),
             Method = rep(colnames(test_results[[i]][-1]), each = nrow(test_results[[i]])),
             Dimension = paste("p =",Dim[i])
  )
}
))
#plot---------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = nb_changes, y = y, color = Method)) +
  facet_grid(cols = vars(Dimension), scales ="fixed") +
  geom_point() +
  geom_line(linetype = "dashed") + 
  ylab("Seconds") +
  xlab("Number of segments into a time series with 10\U2076 data points") +
  scale_x_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  theme_bw() +
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom") +
  scale_color_brewer(palette = 'Set1')

#save---------------------------------------------------------------------|
pdf(file = "Figure9.pdf",  width = 8, height = 4)
print(Plot) 
dev.off()

################################################################################
########################### END ################################################
################################################################################


