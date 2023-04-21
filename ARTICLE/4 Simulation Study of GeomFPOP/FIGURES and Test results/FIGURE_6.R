
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 6:                                                                   |
# Run-time of GeomFROP (S and R types) and PELT                               |
# using multivariate time-series without change-points.                       |
# The maximum run-time of the algorithms is 3 minutes.                        |        
# Averaged over 100 data sets.                                                |
#-----------------------------------------------------------------------------|
# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_2 and TEST_3                        |                                  
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(base)
library(rstream)
library(tidyverse)
library(stats)
library(RColorBrewer)
library(ggpubr)

#parameters---------------------------------------------------------------|
Dim <- c (2, 3, 4)
nameDim <- c("p = 2", 
             "p = 3", 
             "p = 4")

Methods <- c("GeomFPOP (\U0053-type)", 
             "GeomFPOP (\U0052-type)", 
             "PELT")

Algo <- c("GeomFPOP_S_all_all*", 
          "GeomFPOP_R_all_all*",
          "PELT_R_all_all*")

TP <- 1 
Length <-2^(10:25)
NbSimus <- 100

#read results-------------------------------------------------------------|
NbMethods <- length(Algo)
read_file <- NULL
test_results <- list()
for (j in  1 : length(Dim)) {
  test_results[[j]] <- data.frame(n = Length, matrix(0, length(Length), NbMethods))
  colnames(test_results [[j]]) <- c('n', Methods)
  for (k  in 1 :NbMethods) {
    read_file <- paste ('TC', 
                        Dim[j], 
                        TP, 
                        Algo[k],
                        "it", 
                        NbSimus,
                        '.txt', 
                        sep = '_')
    test_results[[j]] [, k+1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}

#transformation of data frame---------------------------------------------|
res <- do.call(rbind, lapply(1:length(test_results), function(i) {
  data.frame(n = test_results[[i]]$n,
             y = c(test_results[[i]][, 2], 
                   test_results[[i]][, 3],
                   test_results[[i]][, 4]),
             Method = rep(colnames(test_results[[i]][-1]), each = nrow(test_results[[i]])),
             tDim = paste("p =",Dim[i])
  )
}
))

#plot---------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = n, y = y, color = Method)) +
  facet_grid(cols = vars(tDim), scales = "fixed") +
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
        legend.text = element_text(size = 12), 
        legend.position = "bottom") +
  scale_color_brewer(palette = 'Set1')
#save---------------------------------------------------------------------|
pdf(file = "Figure6.pdf",  width = 8, height = 4)
print(Plot) 
dev.off()
################################################################################
########################### END ################################################
################################################################################
