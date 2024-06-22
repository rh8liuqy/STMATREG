#' The 'STMATREG' package.
#'
#' @description This R package provides a linear regression based on the matrix variant skew-t distribution with the asynchronous distributed algorithm.
#' @name STMATREG-package
#' @importFrom MASS mvrnorm
#' @import tidyverse
#' @import numDeriv
#' @importFrom parallel makeCluster mclapply
#' @import async
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats pnorm rbinom rexp rgamma rnorm rpois sd runif
#' @importFrom rootSolve uniroot.all
#' @importFrom methods is
#' @useDynLib STMATREG, .registration = TRUE
NULL
#> NULL
