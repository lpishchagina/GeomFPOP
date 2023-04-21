
################################################################################
##  APPENDIX E. Optimization Strategies for GeomFPOP(R-type)                  ##
##  The Number of Potential Change-Point Candidates Stored Over Time          ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 11:                                                                  |
# Run-time of the different optimization approaches of GeomFPOP               |
# p = 2, 3, 4                                                                 |
#-----------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_3 and TEST_7                        |                                  
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
TP <- 1 
Length <-2^(10:25)
NbSimus <- 100
Dim <- c (2, 3, 4)
nameDim <- c("p = 2", 
             "p = 3", 
             "p = 4")
nameMethod <- c("(all / all)", 
                "(all / empty)",
                "(all / random)",
                "(last / all)",
                "(last / empty)",
                "(last / random)",
                "(random / all)",
                "(random / empty)",
                "(random / random)")
Algo <- c("GeomFPOP_R_all_all*",
          "GeomFPOP_R_all_empty",
          "GeomFPOP_R_all_random",
          "GeomFPOP_R_last_all*",
          "GeomFPOP_R_last_empty",
          "GeomFPOP_R_last_random",
          "GeomFPOP_R_random_all*",
          "GeomFPOP_R_random_empty",
          "GeomFPOP_R_random_random")

#read results-------------------------------------------------------------|

NbMethods <- length(Algo)
read_file <- NULL
test_results <- list()
for (j in  1 : length(Dim)) {
  test_results[[j]] <- data.frame(n = Length, 
                                  matrix(0, length(Length), NbMethods))
  colnames(test_results [[j]]) <- c('n', nameMethod)
  for (k  in 1 : NbMethods) {
    read_file <- paste ('TC', 
                        Dim[j], 
                        TP, 
                        Algo[k],
                        'it', 
                        NbSimus,'.txt', 
                        sep='_')
    test_results[[j]] [, k+1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}

# modification of data frame----------------------------------------------|
res <- do.call(rbind, lapply(1:length(test_results), function(i) {
  data.frame(n = test_results[[i]]$n,
             y = c(test_results[[i]][, 2], 
                   test_results[[i]][, 3],
                   test_results[[i]][, 4],
                   test_results[[i]][, 5],
                   test_results[[i]][, 6],
                   test_results[[i]][, 7],
                   test_results[[i]][, 8],
                   test_results[[i]][, 9],
                   test_results[[i]][, 10]),
             Method = rep(colnames(test_results[[i]][-1]), each=nrow(test_results[[i]])),
             tDim = paste("p =",Dim[i])
  )
}
))

#plot--------------------------------------------------------------------|
Plot <- ggplot(res, aes(x = n, y = y, color = Method)) +
  facet_grid(cols = vars(tDim), scales = "fixed") +
  geom_line() + 
  ylab("Seconds") +
  xlab("Number of data points of time series") +
  geom_hline(yintercept = 180, linetype="dashed") + 
  scale_x_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::math_format(10^.x, format = log10)) +
  theme_bw() +
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom") +
  scale_color_brewer(palette = 'Paired')

#save---------------------------------------------------------------------|
pdf(file = "Figure11.pdf",  width = 10, height = 4)
print(Plot)
dev.off()
################################################################################
########################### END ################################################
################################################################################

