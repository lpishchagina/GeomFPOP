
################################################################################
##  4. Simulation Study of GeomFPOP                                           ##
## 4.2. Empirical Time Complexity of GeomFPOP                                 ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 8:                                                                   |
# Time complexity dependence of (random/random) approach                      |
# of GeomFPOP (R-type) on dimension p.                                        |
#-----------------------------------------------------------------------------|
# Remark------------------------------------------------------------------|
# We use the results obtained in TEST_4                                   |                                  
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
NbSimus <- 100
typeA <-c("GeomFPOP (R-type: random / random)")
typeAlgo <-c("GeomFPOP_R_random_random_it")
pentype_ <- 1

#data frame---------------------------------------------------------------|
NbDims <- length(dim_)
read_file <- NULL
test_results <- list()
for (j in  1 : length(typeA)) {
  test_results[[j]] <- data.frame(n = Length, 
                                  matrix(0, length(Length), NbDims))
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

#data=> log10(data)-------------------------------------------------------|
for (i in 1:2){
  test_results[[i]]<-log10(test_results[[i]])
}

#regression---------------------------------------------------------------|
GeomFPOPregres<-list()
alpha<-list()
for (i in 1:length(dim_)){
  GeomFPOPregres[[i]]<- lm(test_results[[1]][, i+1] ~ test_results[[1]][, 1])
  alpha[[i]]<-GeomFPOPregres[[i]]$coefficients[2]
}

#convert list to vector---------------------------------------------------|
alpha <- unlist(alpha)
alpha <- as.vector(alpha,'numeric')
#plot---------------------------------------------------------------------|
table_alpha <- data.frame(dim_, alpha)
colnames(table_alpha) <- c('p', '\u03b1')

Plot<- ggplot(table_alpha, aes(x = table_alpha[, 1], y = table_alpha[, 2])) +
  geom_point( color = '#CB2027') +
  geom_line( linetype = "dashed", color = '#CB2027') +
  ylab("Slope, \u03b1") +
  xlab("Dimension, p") + 
  geom_hline(yintercept = 2, linetype = "dotted", color = '#265DAB') + 
  geom_hline(yintercept = 1, linetype = "dotted",  color = '#265DAB') + 
  theme_bw() + 
  xlim(c(1,100)) +
  ylim(c(1,2.1))
#Plot = Plot + scale_x_continuous(breaks = pretty_breaks(n = 10))
#Plot = Plot + scale_y_continuous(breaks = pretty_breaks(n = 10))

#save---------------------------------------------------------------------|
pdf(file = "Figure8.pdf",  width = 8, height = 4)
print(Plot) 
dev.off()
################################################################################
########################### END ################################################
################################################################################