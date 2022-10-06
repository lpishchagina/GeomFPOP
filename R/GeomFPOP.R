#' @title rnormChanges
#'
#' @description Generation of data (normal distribution) of dimension p with a given values of means and changes
#'
#' @param n number of data point.
#' @param p parameter of dimension.
#' @param changes a vector of increasing changepoint indices (last element is always less then n).
#' By default, 'changes = NULL' (the data without changepoints).
#' @param means matrix of successive means for data, by default 'means = matrix(0, ncol = 1, nrow = p)'.
#' @param noise standard deviation of an additional normal noise, by default 'noise = 1'.
#'
#' @return matrix of data of dimension p x n with a given values of means by the segmentation.
#'
#' @examples
#'N <- 10
#' Chpt <-5
#' Means <-  matrix(c(0,1,1,10), nrow = 2)
#' Noise <- 1
#' Dim <- 2
#' Penality <- 2*Dim*log(N)
#'time_series1 <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
#' time_series2 <- rnormChanges(p = Dim, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = Dim), noise = Noise)

rnormChanges <- function(p, n, changes = NULL, means = matrix(0, ncol = 1, nrow = p), noise = 1) {
  #---stop---#
  if (!is.null(changes) && n <= changes[length(changes)]) {
    stop('The last element of changepoints is always less than n!')
  }
  if(!is.null(changes) && !is.numeric(changes)) {
    stop('changes are not all numeric!')
  }
  if(is.unsorted(changes)) {
    stop('changes should be an increasing vector!')
  }
  if(!is.numeric(means)) {
    stop('means are not all numeric!')
  }
  if ((length(changes) + 1) !=  length(means[1,])) {
    stop('The length of the means[,] is always equal to the number of changes plus one!')
  }
  if (!is.double(noise)) {
    stop('noise is not a double!')
  }
  if (noise < 0) {stop('noise must be non-negative!')
  }
  #---function---#
  res <- matrix(0,p,n)
  InttT <- diff(c(0, changes, n))
  for (i in 1:p) {
    res[i,] <- rep(means[i,], InttT) + rnorm(n, 0, noise)}
  return(res)
}
