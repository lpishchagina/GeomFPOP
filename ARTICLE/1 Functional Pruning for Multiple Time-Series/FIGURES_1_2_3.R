
################################################################################
##  1. Functional Pruning for Multiple Time-Series                            ##
################################################################################


#packages-----------------------------------------------------------------|
#library(devtools)
#devtools::install_github("vrunge/opac", force = TRUE)
#devtools::install_github("lpishchagina/GeomFPOP", force = TRUE)
library(opac)
library(GeomFPOP)
library(tidyverse)
library(plot3D)
library(gridExtra)
library(ggforce)

################################################################################
## 1.2. Functional Pruning Dynamic Programming Algorithm                      ##
################################################################################
#-----------------------------------------------------------------------------|
# Figure 1:                                                                   |
# The sets Z_t^i over time for bi-variate independent Gaussian                |
# time-series without change                                                  |
# y={{0,29;1,86;0,9;−1,26;1,22},{1,93;−0,02;−2,51;0,91;1,11}}.                |
# Index1 associated to quadratics 2 is pruned at iteration 3.                 |
# Each time-sequence of Z_t^i with i fixed is a nested sequence of sets.      |
#-----------------------------------------------------------------------------|
set.seed(617)
# Remark------------------------------------------------------------------|
# Simulation by the package "opac" (V.Runge)                              |                                     
# time series: ts ~ N_2(0, I); dimension: p = 2; number of data: chpts = 5|                                          
#-------------------------------------------------------------------------|

#data generation----------------------------------------------------------|
ts <- dataGenerator_2D(chpts = 5, means1 = 0, means2 = 0)
#ts
#y1         y2
#1  0.2939169  1.9338383
#2  1.8691489 -0.0209778
#3  0.9017728 -2.5184432
#4 -1.2685418  0.9993188
#5  1.2261388  1.1100762

#grid--------------------------------------------------------------------|
mu <- theta_grid_square_2D(ts, bound1 = c(0,3), bound2 = c(0,3), margin = rep(1, 2))
#plot--------------------------------------------------------------------|
pdf(file = "Figure1.pdf",  width = 20, height = 4)
par(mfrow = c(1,5))
plot_2D(ts, mu, color = "ramp")
dev.off()

################################################################################
## 1.3. Geometric Formulation of Functional Pruning                           ##
################################################################################

#-----------------------------------------------------------------------------|
# Figure 2:                                                                   |
# Three examples of the level curves of function s_{ij}                       |
#for bi-variate time-series {x,y}.                                            |
# We use the following simulations for univariate time series :               |
#  (a) x ∼ N (0, 1), y ∼ N (0, 1),                                            |
#  (b) x ∼ P(1), y ∼ P(3),                                                    |
#  (c) x ∼ NB(0.5,1), y ∼ NB(0.8,1).                                          |
#  i = 1, j = 10                                                              |
#-----------------------------------------------------------------------------|

#parameter definition-----------------------------------------------------|
i <- 1
j <- 10
#-----------------------------------------------------|
#(a) x ∼ N (0, 1), y ∼ N (0, 1)
#data generation 
set.seed(21)
ts1_gauss <- rnorm(j)
ts2_gauss <- rnorm(j)
#mean parameter
mean1_gauss <- mean(ts1_gauss)
mean2_gauss <- mean(ts2_gauss)
# mean1_gauss
# [1] 0.1976666
# mean2_gauss 
# [1] 0.07501655
#-----------------------------------------------------| 
#(b) x ∼ P(1), y ∼ P(3),   
#data generation
set.seed(21)
ts1_poisson <- rpois(j, 1)
ts2_poisson <- rpois(j, 2)
#mean parameter
mean1_poisson <- mean(ts1_poisson)
mean2_poisson <- mean(ts2_poisson)
const <- 0
#cost parameters
A <- j - i + 1
B1 <- sum(ts1_poisson)
B2 <- sum(ts2_poisson)
C <- 0
for (u in i : j) { 
  C <- C + 
       log(factorial(ts1_poisson[u])) +
       log(factorial(ts2_poisson[u]))
}
# A
# [1] 10
# B1
# [1] 1.4
# B2
# [1] 1.6
# C
# [1] 13.69304
#-----------------------------------------------------|
# (c) x ∼ NB(0.5,1), y ∼ NB(0.8,1). 
# data generation
set.seed(21)
ts1_neg_bin <- rnbinom(j, prob = 0.5, size = 1)
ts2_neg_bin <- rnbinom(j, prob = 0.8, size = 1)
#cost parameters
NA1 <- sum(ts1_neg_bin)
NA2 <- sum(ts2_neg_bin)
NB <- j - i + 1
NC <- 0
for (u in i : j) {
  NC <- NC + 
        log(factorial(ts1_neg_bin[u] + 9)/(factorial(ts1_neg_bin[u])*factorial(9))) + 
        log(factorial(ts2_neg_bin[u]+9)/(factorial(ts2_neg_bin[u])*factorial(9))) 
}
# NA1 
# [1] 11
# NA2
# [1] 3
# NB
# [1] 10
# NC
#[1] 28.92855

#function definition------------------------------------------------------|
#(a) x ∼ N (0, 1), y ∼ N (0, 1)
x_gauss <- y_gauss <- seq(-2, 2, by = .025)
fun_gauss <- function(x,y) { z <- (x - mean1_gauss)^2 +
                                (y - mean2_gauss)^2 }
z_gauss <- outer(x_gauss, y_gauss, fun_gauss)
r_gauss <- 1:nrow(z_gauss)
p_gauss <- 1:ncol(z_gauss)
#-----------------------------------------------------|
#  (b) x ∼ P(1), y ∼ P(3),
x_poisson <- y_poisson <- seq(0.001,10, by = .025)
fun_poisson <- function(x,y) { z <- A*x + 
                                A*y -  
                                B1*log(x) - 
                                B2*log(y) - 
                                C }
z_poisson <- outer(x_poisson,y_poisson,fun_poisson)
r_poisson <- 1:nrow(z_poisson)
p_poisson <- 1:ncol(z_poisson)
#-----------------------------------------------------|
# (c) x ∼ NB(0.5,1), y ∼ NB(0.8,1). 
x_neg_bin <- y_neg_bin <- seq(0.001,0.99999, by = .025)
fun_neg_bin <- function(x,y) { z <- -NC - 
                                  NB * log(1 - x) - 
                                  NB * log(1 - y) - 
                                  NA1 * log(x) - 
                                  NA2 * log(y) }
z_neg_bin <- outer(x_neg_bin, y_neg_bin, fun_neg_bin)
r_neg_bin <- 1:nrow(z_neg_bin)
p_neg_bin <- 1:ncol(z_neg_bin)

#plot#-----------------------------------------------------|
pdf(file = "Figure2.pdf",  width = 8, height = 6)
par(mfrow=c(2,3))
contour3D(x = x_gauss, 
          y = y_gauss, 
          z = z_gauss,  
          colvar = z_gauss, 
          bty = "b2", dDepth = 1, theta = 60, nlevels = 100)

contour3D(x = x_poisson,
          y = y_poisson,
          z = z_poisson, 
          colvar = z_poisson, 
          bty = "b2",dDepth = 1, theta = 60, nlevels = 100)

contour3D(x = x_neg_bin,
          y = y_neg_bin,
          z = z_neg_bin, 
          colvar = z_neg_bin, 
          bty = "b2",dDepth = 1, theta = 60, nlevels = 100)

contour3D(x = r_gauss,
          y = p_gauss,
          z = 0, 
          colvar = z_gauss, 
          bty = "b2", dDepth = 1, theta = 60, nlevels = 100, colkey = FALSE)

contour3D(x = r_poisson,
          y = p_poisson,
          z = 0, 
          colvar = z_poisson, 
          bty = "b2", dDepth = 1, theta = 60, nlevels = 100, colkey = FALSE)

contour3D(x = r_neg_bin,
          y = p_neg_bin,
          z = 0, 
          colvar = z_neg_bin, 
          bty="b2", dDepth = 1, theta = 60, nlevels = 100, colkey = FALSE)
dev.off()

#-----------------------------------------------------------------------------|
# Figure 3:                                                                   |
# Examples of building a set Z_t^i with |Pi| = |Fi(t)| = 3                    |
# for ts ~ N_2(0, I).                                                         |
# The green disks are S-type sets of the past set P^i.                        |
# The blue disks are S-type sets of the future set F^i(t).                    |
#-----------------------------------------------------------------------------|

#functions-----------------------------------------------------------------|
#S-type set in 2D for Gaussian distribution is the ball (see Remark 2)
getDisk <- function(i, j, x, Cost_t) { ## j > i
  center <- rowMeans(x[, (i+1):j, drop = F])
  radius2 <- (Cost_t[j] - Cost_t[i] - sum(x[, (i+1):j]^2) + (j-i) * sum(center^2)) / (j-i)
  return (list(radius2 = radius2, center = center))
}
#data generation 
set.seed(21)
#p=2, n=10 N(0,I)
ts <- rnormChanges(2, 10, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), 1)
## Cost of one data point
Cost_t <- rowSums(apply(ts, 1, function(y){
  cumsum(y^2) - cumsum(y)^2/(1:length(y))
}))

#generation of 3 past and 2 future disks 
pastDisks <- list()
futureDisks <- list()

past <- 4
tau <- past + 1

stepPS <- c(0, 1, 2)
stepFS <- c(1:3)

for(i in 1:length(stepPS)) pastDisks[[i]] <- getDisk(past - stepPS[i], tau, ts, Cost_t);
for(i in 1:length(stepFS)) futureDisks[[i]]   <- getDisk(tau, tau+stepFS[i]+1, ts, Cost_t);
#table
#future= darkblue , past=darkgreen
r2 <- as.double(unlist(c(futureDisks[[1]][1], futureDisks[[2]][1], futureDisks[[3]][1],
                         pastDisks[[1]][1],pastDisks[[2]][1], pastDisks[[3]][1])))
r <- sqrt(r2)
y <- as.double(unlist(c(futureDisks[[1]][2], futureDisks[[2]][2],futureDisks[[3]][2],
                        pastDisks[[1]][2],pastDisks[[2]][2], pastDisks[[3]][2])))
i = c(1, 3, 5, 7, 9, 11)

circles <- data.frame(
  r = r,
  y1 = y[i],
  y2 = y[i+1],
  alpha = rep(0.01,6),
  colour  = c(rep("darkblue",3),rep("darkgreen",3))
)
#plot
#boundaries
xlim = c(min(circles[,2]-circles[,1]),max(circles[,2]+circles[,1]))
ylim = c(min(circles[,3]-circles[,1]),max(circles[,3]+circles[,1]))

#plot
Plot<-ggplot(circles) +
  geom_circle(aes(x0=y1, y0=y2, r=r, fill=colour, alpha=alpha, color=colour), show.legend=FALSE) + 
  scale_colour_identity() + 
  scale_fill_identity() +
  xlab("y\U00B9") +
  ylab("y\U00B2") +
  xlim(xlim) +
  ylim(ylim) +
  theme_bw()+
  theme(text = element_text(size = 15))

print(Plot)
#pdf(file = "Figure3.pdf",  width = 4, height = 4)
print(Plot)
#dev.off

################################################################################
########################### END ################################################
################################################################################



