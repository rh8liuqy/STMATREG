library(MASS)
library(tidyverse)
library(numDeriv)
library(parallel)
library(async)
library(rootSolve)
# create list_rbind function as tidyverse package in the cluster of out of date
# x: a list of data frame.
list_rbind <- function(x) {
  if (! is(x, "list")) {
    stop("x must be a list.")
  }
  output <- do.call(rbind,x)
  return(output)
}

# Random number generation part -------------------------------------------

# random number generator of the matrix variate normal distribution
# size: the number of random numbers generated
# M: location matrix. n by p matrix.
# Sigma: A scale matrix. n by n matrix.
# Psi: A scale matrix. p by p matrix.
rmatN <- function(size,
                  M,
                  Sigma,
                  Psi) {
  if ((nrow(M) != nrow(Sigma)) | (ncol(M) != ncol(Psi))) {
    stop("Check dimension of M, Sigma and Psi.")
  }
  
  output <- vector(mode = "list",
                   length = size)
  # the location parameter of the multivariate normal distribution
  M_vec <- as.vector(M)
  # the covariance matrix of the multivariate normal distribution
  cov_mat <- kronecker(X = Psi,
                       Y = Sigma)
  # calculate size of matrix
  n <- nrow(Sigma)
  p <- nrow(Psi)
  
  for (i in 1:size) {
    value_vec <- mvrnorm(n = 1,
                         mu = M_vec,
                         Sigma = cov_mat)
    value <- matrix(data = value_vec,
                    nrow = n,
                    ncol = p,
                    byrow = FALSE)
    output[[i]] <- value
  }
  return(output)
}

# random number generator of the matrix variate skewt distribution
# size: the number of random numbers generated
# M: location matrix. n by p matrix.
# A: skewness matrix. n by p matrix.
# Sigma: A scale matrix. n by n matrix.
# Psi: A scale matrix. p by p matrix.
# nu: the degree of freedom.
rmatskewT <- function(size,
                      M,
                      A,
                      Sigma,
                      Psi,
                      nu) {
  output <- vector(mode = "list",
                   length = size)
  for (i in 1:size) {
    W <- 1/rgamma(n = 1,shape = 0.5*nu, rate = 0.5*nu)
    V_mat <- rmatN(size = 1,
                   M = matrix(data = 0,
                              nrow = nrow(A),
                              ncol = ncol(A)),
                   Sigma = Sigma,
                   Psi = Psi)
    V_mat <- V_mat[[1]]
    value <- M + W*A + sqrt(W)*V_mat
    output[[i]] <- value
  }
  return(output)
}

# Create the correlation matrix using the DEC structure
# \operatorname{corr}\left(Y_{i j}, Y_{i k}\right)= \rho^{\left|s_{i k}-s_{i j}\right|^{\theta}}
# x: a numeric vector containing time lapse
DEC_correlation <- function(x,
                            rho,
                            theta) {
  if (rho < 0 | rho > 1) {
    stop("0 <= rho < 1 must hold!")
  }
  if (theta < 0) {
    stop("theta must be postitive")
  }
  dist_mat <- outer(X = x,
                    Y = x, 
                    FUN = function(x,y){return(abs(x-y))})
  dist_mat_theta <- dist_mat^theta
  output <- rho^(dist_mat_theta)
  return(output)
}

# create the Sigma matrix
# ti: a numeric vector.
# DEC_vec: a numeric vector of length 2.
DEC_Sigma_mat_func <- function(ti,
                               DEC_vec) {
  output <- DEC_correlation(x = ti,
                            rho = DEC_vec[1],
                            theta = DEC_vec[2])
  return(output)
  
}

diagonal_Sigma_mat_func <- function(ti) {
  ni <- length(ti)
  output <- diag(x = 1, nrow = ni, ncol = ni)
  return(output)
}


# loglikelihood part ------------------------------------------------------

# log-likelihood from MST package for i-th observation
# Y: n by p matrix. observed values.
# M: n by p matrix. location parameter.
# A: n by p matrix. skewness parameter.
# Sigma: n by n matrix. A covariance matrix.
# Psi: p by p matrix. A covariance matrix.
# n: the number of rows in Y.
# p: the number of columns in Y.
# log_pdf_type: the type of loglikelihood
MV_logdensity <- function(Y, 
                          M, 
                          A, 
                          Sigma, 
                          Psi, 
                          nu, 
                          n, 
                          p, 
                          log_pdf_type) {
  if (! log_pdf_type %in% c("skewT", "skewN", "N")) {
    stop("log_pdf_type must be one of 'skewT', 'skewN' and 'N'.")
  }
  if (log_pdf_type != "N") {
    lambda <- -(nu + n * p)/2
    Delta <- sum(diag((solve(Sigma) %*% (Y - M) %*% solve(Psi) %*% t(Y - M)))) + nu
    rho <- sum(diag(solve(Sigma) %*% A %*% solve(Psi) %*% t(A)))
    BessX <- besselK(sqrt(Delta*rho),lambda)
    # for numerical stability
    if (identical(BessX,0)){
      BessX <- 1e-100
    }
    
    dens <-log(2)+(nu/2)*log(nu/2)-log(2*pi)*(n*p/2)-log(det(Sigma))*(p/2)-log(det(Psi))*(n/2)-log(gamma(nu/2))+(lambda/2)*(log(Delta)-log(rho))+log(BessX)+sum(diag(solve(Sigma)%*%(Y-M)%*%solve(Psi)%*%t(A)))
    
    return(dens)
  } else {
    YMinusM <- Y - M
    Psi_Sigma <- kronecker(Psi,Sigma)
    logDens <- mvtnorm::dmvnorm(x = as.numeric(YMinusM),
                                sigma = Psi_Sigma,
                                log = TRUE)
    return(logDens)
  }
}

# Y: a list of length N.
# M: a list of length N.
# A: a list of length N.
# Sigma: a list of length N.
# Psi: a p by p matrix.
# nu: a number.
# log_pdf_type: the type of loglikelihood.
# log-likelihood for all observations
MV_logpdf <- function(Y, 
                      M, 
                      A, 
                      Sigma, 
                      Psi, 
                      nu,
                      log_pdf_type) {
  if (! log_pdf_type %in% c("skewT", "skewN", "N")) {
    stop("log_pdf_type must be one of 'skewT', 'skewN' and 'N'.")
  }
  N <- length(Y)
  output <- numeric(N)
  for (i in 1:N) {
    Yi <- Y[[i]]
    Mi <- M[[i]]
    Ai <- A[[i]]
    Sigmai <- Sigma[[i]]
    ni <- nrow(Yi)
    p_i <- ncol(Yi)
    output[i] <- MV_logdensity(Y = Yi, 
                               M = Mi,
                               A = Ai,
                               Sigma = Sigmai, 
                               Psi = Psi,
                               nu = nu,
                               n = ni, 
                               p = p_i,
                               log_pdf_type)
  }
  output <- sum(output)
  return(output)
}

create_symmetric_matrix <- function(upper_triangular_vector) {
  n <- floor((-1 + sqrt(1 + 8 * length(upper_triangular_vector))) / 2)
  
  symmetric_matrix <- matrix(0, nrow = n, ncol = n)
  
  k <- 1
  for (i in 1:n) {
    for (j in i:n) {
      symmetric_matrix[i, j] <- upper_triangular_vector[k]
      symmetric_matrix[j, i] <- upper_triangular_vector[k]
      k <- k + 1
    }
  }
  
  return(symmetric_matrix)
}

# Used for calculating the hessian matrix at the end of EM
# Y: a list of length N.
# X: a list of length N.
# ti: a list of length N.
# log_pdf_type: the type of loglikelihood.
MV_logpdf_hess <- function(parameter,
                           Y,
                           X,
                           ti,
                           Sigma_type,
                           log_pdf_type) {
  if (! (Sigma_type %in% c("DEC","diagonal"))) {
    stop('Sigma_type %in% c("DEC","diagonal") must be TRUE.')
  }
  if (! log_pdf_type %in% c("skewT", "skewN", "N")) {
    stop("log_pdf_type must be one of 'skewT', 'skewN' and 'N'.")
  }
  N <- length(Y)
  p <- ncol(Y[[1]])
  q <- ncol(X[[1]])
  beta <- parameter[1:(q*p)]
  beta_mat <- matrix(data = beta, nrow = q, ncol = p, byrow = FALSE)
  if (log_pdf_type == "skewT") {
    A_vec <- parameter[(q*p+1):(q*p+p)]
    if (Sigma_type == "DEC") {
      DEC_vec <- parameter[(q*p+p+1):(q*p+p+2)]
      Psi_vec <- parameter[(q*p+p+3):(q*p+p+2+p*(p+1)/2)]
      nu <- parameter[(q*p+p+2+p*(p+1)/2)+1]
    } else if (Sigma_type == "diagonal") {
      Psi_vec <- parameter[(q*p+p+1):(q*p+p+p*(p+1)/2)]
      nu <- parameter[(q*p+p+p*(p+1)/2) + 1]
    }
    
    # create Psi matrix
    Psi_mat_value <- create_symmetric_matrix(upper_triangular_vector = Psi_vec)
    
    # initialize other matrix list
    Sigma_mat_value <- vector(mode = "list",
                              length = N)
    A_mat_value <- vector(mode = "list",
                          length = N)
    M_mat_value <- vector(mode = "list",
                          length = N)
    
    for (i in 1:N) {
      # create M matrix
      Xi <- X[[i]]
      M_mat_value[[i]] <- Xi %*% beta_mat
      
      # create A matrix
      Yi <- Y[[i]]
      ni <- nrow(Yi)
      A_mat_value[[i]] <- matrix(data = rep(A_vec,ni),
                                 nrow = ni,
                                 byrow = TRUE)
      
      # create Sigma matrix
      tii <- ti[[i]]
      if (Sigma_type == "DEC") {
        Sigma_mat_value[[i]] <- DEC_Sigma_mat_func(ti = tii,
                                                   DEC_vec = DEC_vec)
      } else {
        Sigma_mat_value[[i]] <- diagonal_Sigma_mat_func(ti = tii)
      }
    }
    
    output <- MV_logpdf(Y = Y, 
                        M = M_mat_value, 
                        A = A_mat_value, 
                        Sigma = Sigma_mat_value, 
                        Psi = Psi_mat_value, 
                        nu = nu,
                        log_pdf_type = log_pdf_type)
    return(output)
  } else if (log_pdf_type == "skewN") {
    A_vec <- parameter[(q*p+1):(q*p+p)]
    if (Sigma_type == "DEC") {
      DEC_vec <- parameter[(q*p+p+1):(q*p+p+2)]
      Psi_vec <- parameter[(q*p+p+3):(q*p+p+2+p*(p+1)/2)]
      nu <- 200
    } else if (Sigma_type == "diagonal") {
      Psi_vec <- parameter[(q*p+p+1):(q*p+p+p*(p+1)/2)]
      nu <- 200
    }
    
    # create Psi matrix
    Psi_mat_value <- create_symmetric_matrix(upper_triangular_vector = Psi_vec)
    
    # initialize other matrix list
    Sigma_mat_value <- vector(mode = "list",
                              length = N)
    A_mat_value <- vector(mode = "list",
                          length = N)
    M_mat_value <- vector(mode = "list",
                          length = N)
    
    for (i in 1:N) {
      # create M matrix
      Xi <- X[[i]]
      M_mat_value[[i]] <- Xi %*% beta_mat
      
      # create A matrix
      Yi <- Y[[i]]
      ni <- nrow(Yi)
      A_mat_value[[i]] <- matrix(data = rep(A_vec,ni),
                                 nrow = ni,
                                 byrow = TRUE)
      
      # create Sigma matrix
      tii <- ti[[i]]
      if (Sigma_type == "DEC") {
        Sigma_mat_value[[i]] <- DEC_Sigma_mat_func(ti = tii,
                                                   DEC_vec = DEC_vec)
      }
      else {
        Sigma_mat_value[[i]] <- diagonal_Sigma_mat_func(ti = tii)
      }
    }
    
    output <- MV_logpdf(Y = Y, 
                        M = M_mat_value, 
                        A = A_mat_value, 
                        Sigma = Sigma_mat_value, 
                        Psi = Psi_mat_value, 
                        nu = nu,
                        log_pdf_type = log_pdf_type)
    return(output)
  } else if (log_pdf_type == "N") {
    A_vec <- rep(0,p)
    if (Sigma_type == "DEC") {
      DEC_vec <- parameter[(q*p+1):(q*p+2)]
      Psi_vec <- parameter[(q*p+3):(q*p+2+p*(p+1)/2)]
      nu <- 200
    } else if (Sigma_type == "diagonal") {
      Psi_vec <- parameter[(q*p+1):(q*p+p*(p+1)/2)]
      nu <- 200
    }
    
    # create Psi matrix
    Psi_mat_value <- create_symmetric_matrix(upper_triangular_vector = Psi_vec)
    
    # initialize other matrix list
    Sigma_mat_value <- vector(mode = "list",
                              length = N)
    A_mat_value <- vector(mode = "list",
                          length = N)
    M_mat_value <- vector(mode = "list",
                          length = N)
    
    for (i in 1:N) {
      # create M matrix
      Xi <- X[[i]]
      M_mat_value[[i]] <- Xi %*% beta_mat
      
      # create A matrix
      Yi <- Y[[i]]
      ni <- nrow(Yi)
      A_mat_value[[i]] <- matrix(data = rep(A_vec,ni),
                                 nrow = ni,
                                 byrow = TRUE)
      
      # create Sigma matrix
      tii <- ti[[i]]
      if (Sigma_type == "DEC") {
        Sigma_mat_value[[i]] <- DEC_Sigma_mat_func(ti = tii,
                                                   DEC_vec = DEC_vec)
      } else {
        Sigma_mat_value[[i]] <- diagonal_Sigma_mat_func(ti = tii)
      }
    }
    
    output <- MV_logpdf(Y = Y, 
                        M = M_mat_value, 
                        A = A_mat_value, 
                        Sigma = Sigma_mat_value, 
                        Psi = Psi_mat_value, 
                        nu = nu,
                        log_pdf_type = log_pdf_type)
    return(output)
  }
}

# Regression data simulation part -----------------------------------------

#' Simulation program for the regression model
#'
#' @param N The number of subject.
#' @param beta_mat The predictor coefficients (q by p matrix). q is the number of predictors. p is the number of columns of Yi.
#' @param A_vec A numeric vector whose length is p.
#' @param DEC_vec A numeric vector whose length is 2.
#' @param Psi_mat A p by p matrix.
#' @param nu The degree of freedom.
#' @param Sigma_type The type of the correlation matrix. Must be "DEC" or "diagonal".
#' @param log_pdf_type Must be "skewT" or "N".
#'
#' @return The list of simulated data.
#' @export
#'
#' @examples
#' N <- 2500
#' beta_mat <- matrix(c(0.5,1.5,-0.5,0.5,1.5,-0.5),nrow = 3)
#' A_vec <- c(2.0,-2.0)
#' DEC_vec <- c(0.9,0.8)
#' Psi_mat <- matrix(c(1.0,-0.5,-0.5,1.0), nrow = 2)
#' nu <- 5.0
#' Sigma_type <- "DEC"
#' log_pdf_type <- "skewT"
#' simulated_data <- reg_simulation(N,
#'                                  beta_mat,
#'                                  A_vec,
#'                                  DEC_vec,
#'                                  Psi_mat,
#'                                  nu,
#'                                  Sigma_type,
#'                                  log_pdf_type)
#' print(str(simulated_data))
reg_simulation <- function(N,
                           beta_mat,
                           A_vec,
                           DEC_vec,
                           Psi_mat,
                           nu,
                           Sigma_type,
                           log_pdf_type) {
  if (! (Sigma_type %in% c("DEC","diagonal"))) {
    stop('Sigma_type %in% c("DEC","diagonal") must be TRUE.')
  }
  
  if (! (log_pdf_type %in% c("skewT","N"))) {
    stop('log_pdf_type %in% c("skewT","N") must be TRUE.')
  }
  
  output_Y <- vector(mode = "list", length = N)
  output_X <- vector(mode = "list", length = N)
  output_ti <- vector(mode = "list", length = N)
  output_Tvec <- vector(mode = "list", length = N)
  output_ni <- vector(mode = "list", length = N)
  q <- nrow(beta_mat)
  if (q != 3) {
    stop("The number of row of beta_mat must be 3 !")
  }
  p <- ncol(beta_mat)
  
  for (i in 1:N) {
    ni <- rpois(n = 1, lambda = 8) + 2
    ti <- rnorm(n = ni)
    Xi <- data.frame(x1 = rexp(n = ni),
                     x2 = rnorm(n = ni),
                     x3 = rbinom(n = ni, size = 1, prob = 2*pnorm(abs(ti)) - 1))
    Xi <- as.matrix(Xi)
    output_X[[i]] <- Xi
    output_ti[[i]] <- ti
    output_ni[[i]] <- ni
  }
  
  # calculate the location and scaling coefficients
  
  X_mat <- do.call(rbind,output_X)
  loc1 <- mean(X_mat[,1])
  scale1 <- sd(X_mat[,1])
  loc2 <- mean(X_mat[,2])
  scale2 <- sd(X_mat[,2])
  loc3 <- mean(X_mat[,3])
  scale3 <- sd(X_mat[,3])
  
  for (i in 1:N) {
    ni <- output_ni[[i]]
    ti <- output_ti[[i]]
    Xi <- output_X[[i]]
    
    Xi[,1] <- (Xi[,1] - loc1)/scale1
    Xi[,2] <- (Xi[,2] - loc2)/scale2
    Xi[,3] <- (Xi[,3] - loc3)/scale3
    
    etai <- Xi%*%beta_mat
    
    if (log_pdf_type == "skewT") {
      if (Sigma_type == "DEC") {
        Ui <- rmatskewT(size = 1,
                        M = matrix(data = 0, nrow = ni, ncol = p),
                        A = matrix(data = rep(A_vec,ni), nrow = ni, ncol = p , byrow = TRUE),
                        Sigma = DEC_correlation(x = ti,
                                                rho = DEC_vec[1],
                                                theta = DEC_vec[2]),
                        Psi = Psi_mat,
                        nu = nu)
      } else {
        Ui <- rmatskewT(size = 1,
                        M = matrix(data = 0, nrow = ni, ncol = p),
                        A = matrix(data = rep(A_vec,ni), nrow = ni, ncol = p , byrow = TRUE),
                        Sigma = diagonal_Sigma_mat_func(ti = ti),
                        Psi = Psi_mat,
                        nu = nu)
      }
    } else if (log_pdf_type == "N") {
      
      if (Sigma_type == "DEC") {
        Ui <- rmatN(size = 1,
                    M = matrix(data = 0, nrow = ni, ncol = p),
                    Sigma = DEC_correlation(x = ti,
                                            rho = DEC_vec[1],
                                            theta = DEC_vec[2]),
                    Psi = Psi_mat)
      } else {
        Ui <- rmatN(size = 1,
                    M = matrix(data = 0, nrow = ni, ncol = p),
                    Sigma = diagonal_Sigma_mat_func(ti = ti),
                    Psi = Psi_mat)
      }
      
    }
    Ui <- Ui[[1]]
    Yi <- etai + Ui 
    output_Y[[i]] <- Yi
    output_X[[i]] <- Xi
    
  }
  output <- list(Y = output_Y,
                 X = output_X,
                 ti = output_ti,
                 ni = output_ni)
  return(output)
}

# Distribued EM algorithm part --------------------------------------------


# N: the number of subjects
# M: the number of subsets
# output of this function: a list of indices of M subsets
subsets_indices <- function(N, M) {
  # Initialize subset sizes and indices
  subset_sizes <- rep(0, M)
  subset_indices <- vector("list", M)
  
  # Generate random indices for the data
  data_indices <- sample(1:N, N)
  
  for (i in 1:N) {
    subset_idx <- i %% M + 1
    subset_indices[[subset_idx]] <- c(subset_indices[[subset_idx]], data_indices[i])
    subset_sizes[subset_idx] <- subset_sizes[subset_idx] + 1
  }
  
  return(subset_indices)
}

#' The function create subset based on the whole data set.
#'
#' @param Y The list of response matrices.
#' @param X The list of design matrices.
#' @param ti The list of measurement times.
#' @param N_subset The number of subsets.
#'
#' @return A list of subsets.
#' @noRd
create_subset <- function(Y,
                          X,
                          ti,
                          N_subset) {
  N <- length(Y)
  output <- vector(mode = "list", length = N_subset)
  indices_list <- subsets_indices(N = N, M = N_subset)
  for (i in 1:N_subset) {
    indecei <- indices_list[[i]]
    output[[i]] <- list(Y = Y[indecei],
                        X = X[indecei],
                        ti = ti[indecei])
  }
  return(output)
}

# delta function in the E step
delta_func <- function(y_mat,
                       M_mat,
                       Sigma_inv,
                       Psi_inv) {
  value <- Sigma_inv%*%(y_mat-M_mat)%*%Psi_inv%*%t((y_mat-M_mat))
  output <- sum(diag(value))
  return(output)
}

# rho function in the E step
rho_func <- function(A_mat,
                     Sigma_inv,
                     Psi_inv) {
  value <- Sigma_inv%*%A_mat%*%Psi_inv%*%t(A_mat)
  output <- sum(diag(value))
  return(output)
}

# ai function in the E step
ai_func <- function(delta,
                    rho,
                    nu,
                    K_lambda_1,
                    K_lambda_0) {
  p1 <- sqrt((delta + nu)/rho)
  p2 <- K_lambda_1/K_lambda_0
  output <- p1*p2
  return(output)
}

# bi function in the E step
bi_func <- function(delta,
                    rho,
                    nu,
                    K_lambda_1,
                    K_lambda_0,
                    n,
                    p) {
  p1 <- sqrt(rho/(delta + nu))
  p2 <- K_lambda_1/K_lambda_0
  p3 <- nu + n*p
  p4 <- delta + nu
  output <- p1*p2 + p3/p4
  return(output)
}

# ci function in the E step
ci_func <- function(delta,
                    rho,
                    nu,
                    K_lambda_0,
                    derivative_K) {
  p1 <- log(sqrt((delta + nu)/rho))
  p2 <- 1/K_lambda_0
  p3 <- derivative_K
  output <- p1+p2*p3
  return(output)
}

# E step

Estep <- function(Y_list,
                  X_list,
                  ti_list,
                  beta_mat,
                  A_mat_list,
                  Sigma_inv_list,
                  Psi_inv,
                  nu,
                  log_pdf_type,
                  EM_type,
                  Sigma_type,
                  DEC2_lb,
                  DEC2_ub,
                  DEC_value) {
  N <- length(Y_list)
  output <- vector(mode = "list", length = 3)
  names(output) <- c("a","b","c")
  output$a <- numeric(N)
  output$b <- numeric(N)
  output$c <- numeric(N)
  
  # output for updating beta
  output$beta_p1 <- vector(mode = "list", length = N)
  output$beta_p2 <- vector(mode = "list", length = N)
  
  # output for updating nu
  output$nu_statistics <- numeric(N)
  output$nu_N <- N
  
  # output for updating A
  output$tJ_Sigma_inv <- vector(mode = "list", length = N)
  output$SA1 <- vector(mode = "list", length = N)
  output$SA2 <- vector(mode = "list", length = N)
  
  # output for updating Psi
  if (EM_type != "ADECM"){
    output$Sigma_inv <- vector(mode = "list", length = N) 
  }
  output$ni <- numeric(N)
  output$SPsi <- vector(mode = "list", length = N)
  
  # Psi_value for the evaluation of log-likelihood function
  Psi_value <- solve(Psi_inv)
  
  if (Sigma_type == "DEC") {
    # number of grid for approximation
    n.grid <- 11
    # output for updating DEC parameters
    output$DEC1_mat <- matrix(data = NA, nrow = n.grid, ncol = N) 
    output$DEC2_mat <- matrix(data = NA, nrow = n.grid, ncol = N) 
    DEC1_grid <- seq(0,1,length.out = n.grid)
    # for numerical stability
    DEC1_grid[1] <- DEC1_grid[1] + 1e-5
    DEC1_grid[n.grid] <- DEC1_grid[n.grid] - 1e-5
    DEC2_grid <- seq(DEC2_lb,DEC2_ub,length.out = n.grid)
    # for numerical stability
    DEC2_grid[1] <- DEC2_grid[1] + 1e-5
    DEC2_grid[n.grid] <- DEC2_grid[n.grid] - 1e-5
    output$DEC1_objective_grid <- numeric(n.grid) # In M step, one can use this output to pick the DEC1 value that maximizes the objective function
    output$DEC2_objective_grid <- numeric(n.grid) # In M step, one can use this output to pick the DEC2 value that maximizes the objective function
    output$DEC1_grid <- DEC1_grid
    output$DEC2_grid <- DEC2_grid
  }
  
  for (i in 1:N) {
    yi_mat <- Y_list[[i]]
    Xi_mat <- X_list[[i]]
    Mi_mat <- Xi_mat%*%beta_mat
    yi_minus_Mi <- yi_mat - Mi_mat
    Sigmai_inv <- Sigma_inv_list[[i]]
    Ai_mat <- A_mat_list[[i]]
    ni <- nrow(yi_mat)
    p_i <- ncol(yi_mat)
    lambdai <- -0.5*(nu + ni*p_i)
    deltai_trace <- delta_func(y_mat = yi_mat,
                               M_mat = Mi_mat,
                               Sigma_inv = Sigmai_inv,
                               Psi_inv = Psi_inv)
    rhoi_trace <- rho_func(A_mat = Ai_mat,
                           Sigma_inv = Sigmai_inv,
                           Psi_inv = Psi_inv)
    kai <- sqrt(rhoi_trace*(deltai_trace + nu))
    
    K_lambdai_1 <- besselK(x = kai, 
                           nu = lambdai + 1)
    K_lambdai_0 <- besselK(x = kai, 
                           nu = lambdai + 0)
    
    if (is.infinite(K_lambdai_1)) {
      K_lambdai_1 <- 1e5
    }
    if (is.infinite(K_lambdai_0)) {
      K_lambdai_0 <- 1e5
    }
    
    # the derivative of the modified Bessel function of the third kind
    mod_Bessel <- function(lambda) {
      output <- besselK(x = kai,
                        nu = lambda)
    }
    
    # in the matrix-variant normal case, derivative_K = 1.
    derivative_K <- tryCatch({grad(func = mod_Bessel,
                                   x = lambdai + 0,
                                   method="Richardson")},
                             error = function(e) {
                               1
                             })
    
    # for numerical stability
    if (identical(K_lambdai_1,0)) {
      K_lambdai_1 <- 1e-11
    }
    if (identical(K_lambdai_0,0)) {
      K_lambdai_0 <- 1e-11
    }
    
    output$a[i] <- ai_func(delta = deltai_trace,
                           rho = rhoi_trace,
                           nu = nu,
                           K_lambda_1 = K_lambdai_1,
                           K_lambda_0 = K_lambdai_0)
    if (is.infinite(output$a[i])) {
      output$a[i] <- 1e5
    }
    
    output$b[i] <- bi_func(delta = deltai_trace,
                           rho = rhoi_trace,
                           nu = nu,
                           K_lambda_1 = K_lambdai_1,
                           K_lambda_0 = K_lambdai_0,
                           n = ni,
                           p = p_i)
    
    output$c[i] <- ci_func(delta = deltai_trace,
                           rho = rhoi_trace,
                           nu = nu,
                           K_lambda_0 = K_lambdai_0,
                           derivative_K = derivative_K)
    
    if (is.infinite(output$c[i])) {
      output$c[i] <- 1e5
    }
    
    # Sufficient Statistics for beta
    bi <- output$b[i]
    tXSigmai <- t(Xi_mat) %*% Sigmai_inv
    output$beta_p1[[i]] <- bi * (tXSigmai %*% Xi_mat)
    output$beta_p2[[i]] <- -tXSigmai %*% Ai_mat + bi*(tXSigmai %*% yi_mat)
    
    # Sufficient Statistics for nu
    ci <- output$c[i]
    output$nu_statistics[i] <- bi + ci
    
    # Sufficient Statistics for A
    Ji <- matrix(1,nrow = ni,ncol = 1)
    tJi <- t(Ji)
    output$tJ_Sigma_inv[[i]] <- tJi %*% Sigmai_inv
    output$SA1[[i]] <- output$tJ_Sigma_inv[[i]] %*% yi_minus_Mi
    output$SA2[[i]] <- output$a[i] * sum(Sigmai_inv)
    
    # Sufficient Statistics for Psi
    ai <- output$a[i]
    if (EM_type != "ADECM"){
      output$Sigma_inv[[i]]  <- Sigmai_inv 
    }
    output$ni[i] <- ni
    tyi_minus_Mi <- t(yi_minus_Mi)
    temp1 <- bi* (tyi_minus_Mi %*% Sigmai_inv %*% yi_minus_Mi)
    tAi_mat <- t(Ai_mat)
    temp2 <- -tAi_mat %*%  Sigmai_inv %*% yi_minus_Mi
    temp3 <- -tyi_minus_Mi %*% Sigmai_inv %*% Ai_mat
    temp4 <- ai * (tAi_mat %*% Sigmai_inv %*% Ai_mat)
    output$SPsi[[i]] <- temp1 + temp2 + temp3 + temp4
    
    # Sufficient Statistics for DEC(1) parameters
    if (EM_type == "ADECM") {
      if (Sigma_type == "DEC") {
        ti <- ti_list[[i]]
        Amati <- A_mat_list[[i]]
        tAmati <- t(Amati)
        for (j in 1:n.grid) {
          DEC_value_grid <- c(DEC1_grid[j],DEC_value[2])
          Sigmaj <- Sigma_mat_func(ti = ti,
                                   DEC_vec = DEC_value_grid) 
          # for numerical stability

          output$DEC1_mat[j,i] <- tryCatch({MV_logdensity(Y = yi_mat,
                                                          M = Mi_mat,
                                                          A = Ai_mat, 
                                                          Sigma = Sigmaj, 
                                                          Psi = Psi_value, 
                                                          nu = nu, 
                                                          n = nrow(yi_mat), 
                                                          p = ncol(yi_mat), 
                                                          log_pdf_type = log_pdf_type)},
                                           error = function(x) {
                                             return(-Inf)
                                           })
          
        }
      } 
    }
    
    # Sufficient Statistics for DEC(2) parameters
    if (EM_type == "ADECM") {
      if (Sigma_type == "DEC") {
        ti <- ti_list[[i]]
        Amati <- A_mat_list[[i]]
        tAmati <- t(Amati)
        for (j in 1:n.grid) {
          DEC_value_grid <- c(DEC_value[1],DEC2_grid[j])
          Sigmaj <- Sigma_mat_func(ti = ti,
                                   DEC_vec = DEC_value_grid) 
          # for numerical stability
          
          output$DEC2_mat[j,i] <- tryCatch({MV_logdensity(Y = yi_mat,
                                                          M = Mi_mat,
                                                          A = Ai_mat, 
                                                          Sigma = Sigmaj, 
                                                          Psi = Psi_value, 
                                                          nu = nu, 
                                                          n = nrow(yi_mat), 
                                                          p = ncol(yi_mat), 
                                                          log_pdf_type = log_pdf_type)},
                                           error = function(x) {
                                             return(-Inf)
                                           })

        }
      } 
    }
    
  }
  # Summation of Sufficient Statistics for beta
  output$beta_p1 <- Reduce("+",output$beta_p1)
  output$beta_p2 <- Reduce("+",output$beta_p2)
  # Summation of Sufficient Statistics for nu
  output$nu_statistics <- sum(output$nu_statistics)
  # Summation of Sufficient Statistics for A
  output$SA1 <- Reduce("+",output$SA1)
  output$SA2 <- Reduce("+",output$SA2)
  # Summation of Sufficient Statistcs for Psi
  output$SPsi <- Reduce("+",output$SPsi)
  output$SPsi_n <- sum(output$ni)
  # Summation of Sufficient Statistics for DEC
  if (Sigma_type == "DEC") {
    output$DEC1_objective_grid <- rowSums(output$DEC1_mat)
    output$DEC2_objective_grid <- rowSums(output$DEC2_mat)
  }
  
  if (EM_type == "ADECM") {
    if (Sigma_type != "DEC") {
      output <- output[c("beta_p1","beta_p2","nu_statistics","nu_N","SA1","SA2","SPsi","SPsi_n")] 
    } else if (Sigma_type == "DEC") {
      output <- output[c("beta_p1","beta_p2","nu_statistics","nu_N","SA1","SA2",
                         "SPsi","SPsi_n",
                         "DEC1_objective_grid",
                         "DEC2_objective_grid",
                         "DEC1_grid",
                         "DEC2_grid")] 
    }
  }
  return(output)
}

# M steps

# update beta
update_beta <- function(Estep_output) {
  N_subset <- length(Estep_output)
  q <- ncol(Estep_output[[1]]$beta_p1)
  p <- ncol(Estep_output[[1]]$beta_p2)
  p1 <- matrix(0, nrow = q, ncol = q)
  p2 <- matrix(0, nrow = q, ncol = p)
  for (i in 1:N_subset) {
    p1 <- p1 + Estep_output[[i]]$beta_p1
    p2 <- p2 + Estep_output[[i]]$beta_p2
  }
  output <- solve(a = p1, b = p2)
  return(output)
}

# update nu
update_nu <- function(Estep_output) {
  N_subset <- length(Estep_output)
  nu_statistics <- 0
  nu_N <- 0
  for (i in 1:N_subset) {
    nu_statistics <- nu_statistics + Estep_output[[i]]$nu_statistics
    nu_N <- nu_N + Estep_output[[i]]$nu_N
  }
  obj_func <- function(nu) {
    p1 <- log(nu*0.5)
    p2 <- 1
    p3 <- -digamma(0.5*nu)
    p4 <- -(1/nu_N)*nu_statistics
    output <- p1 + p2 + p3 + p4
    return(output)
  }
  # the maximum value of nu is 100
  output <- uniroot.all(f = obj_func,
                        interval = c(0.01,1000))
  
  if (length(output) == 0) {
    print("Fail to update the degrees of freedom.")
    print("Set nu to be 2.")
    output <- 2.0
  } else if (length(output) > 1) { # in just there is more than one root
    print("More than one root for nu has been found. Keep the smallest value.")
    output <- min(output)
  }
  return(output)
}

# update A - regular updating mechanism for regular ECM and regular parallel ECM only
update_A_regular <- function(Estep_output, beta_mat, subset_list) {
  N_subset <- length(Estep_output)
  p <- ncol(Estep_output[[1]]$beta_p2)
  A_numerator <- rep(0,p)
  A_denominator <- 0
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    # the number of subject in subset i
    Ni <- length(Estep_output_subseti$tJ_Sigma_inv)
    for (j in 1:Ni) {
      Yj <- subset_list[[i]]$Y[[j]]
      Xj <- subset_list[[i]]$X[[j]]
      tJ_Sigma_inv <- Estep_output_subseti$tJ_Sigma_inv[[j]]
      yjstar <- Yj - Xj %*% beta_mat
      aj <- Estep_output_subseti$a[j]
      A_numerator <- A_numerator + as.numeric(tJ_Sigma_inv %*% yjstar)
      A_denominator <- A_denominator + aj * sum(tJ_Sigma_inv)
    }
  }
  output <- A_numerator/A_denominator
  return(output)
}

# update A in ADECM
update_A_ADECM <- function(Estep_output) {
  N_subset <- length(Estep_output)
  A_dim <- length(as.numeric(Estep_output[[1]]$SA1))
  SA1 <- numeric(A_dim)
  SA2 <- 0
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    SA1 <- SA1 +  as.numeric(Estep_output_subseti$SA1)
    SA2 <- SA2 + Estep_output_subseti$SA2
  }
  output <- SA1/SA2
  return(output)
}

# update Psi - regular updating mechanism for regular ECM and regular parallel ECM only
update_Psi_regular <- function(Estep_output, beta_mat, A_value, subset_list) {
  N_subset <- length(Estep_output)
  p <- ncol(Estep_output[[1]]$beta_p2)
  p1 <- matrix(data = 0, nrow = p, ncol = p)
  p2 <- matrix(data = 0, nrow = p, ncol = p)
  p3 <- matrix(data = 0, nrow = p, ncol = p)
  p4 <- matrix(data = 0, nrow = p, ncol = p)
  p5 <- 0
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    # the number of subject in subset i
    Ni <- length(Estep_output_subseti$Sigma_inv)
    for (j in 1:Ni) {
      Yj <- subset_list[[i]]$Y[[j]]
      Xj <- subset_list[[i]]$X[[j]]
      Mj <- Xj %*% beta_mat
      yjstar <- Yj - Mj
      tyjstar <- t(yjstar)
      nj <- Estep_output_subseti$ni[j]
      Amatj <- matrix(data = rep(A_value,nj), nrow = nj, byrow = TRUE)
      tAmatj <- t(Amatj)
      aj <- Estep_output_subseti$a[j]
      bj <- Estep_output_subseti$b[j]
      Sigma_invj <- Estep_output_subseti$Sigma_inv[[j]]
      tyjstar_Sigma_invj <- tyjstar %*% Sigma_invj
      tAmatj_Sigma_invj <- tAmatj %*% Sigma_invj
      p1 <- p1 + bj*(tyjstar_Sigma_invj %*% yjstar)
      p2 <- p2 - tAmatj_Sigma_invj %*% yjstar
      p3 <- p3 - tyjstar_Sigma_invj %*% Amatj
      p4 <- p4 + aj*(tAmatj_Sigma_invj %*% Amatj)
      p5 <- p5 + nj
    }
  }
  output <- p1 + p2 + p3 + p4
  output <- output/p5
  return(output)
}

update_Psi_ADECM <- function(Estep_output) {
  N_subset <- length(Estep_output)
  Psi_nrow <- nrow(Estep_output[[1]]$SPsi)
  Psi_ncol <- ncol(Estep_output[[1]]$SPsi)
  SPsi1 <- matrix(0,nrow = Psi_nrow,ncol = Psi_ncol)
  SPsi2 <- 0
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    SPsi1 <- SPsi1 + Estep_output[[i]]$SPsi
    SPsi2 <- SPsi2 + Estep_output[[i]]$SPsi_n
  }
  output <- SPsi1/SPsi2
  return(output)
}

# create the Sigma matrix
# ti: a numeric vector.
# DEC_vec: a numeric vector of length 2.
Sigma_mat_func <- function(ti,
                           DEC_vec) {
  output <- DEC_correlation(x = ti,
                            rho = DEC_vec[1],
                            theta = DEC_vec[2])
  return(output)
  
}

# DEC1: a number.
# DEC2: a number.
# beta_mat: a q by p matrix. 
# A_value: a numeric vector.
# Psi_inv_value: a p by p matrix.
# nu_value: a number.
# subset: One subset from the whole data set.
# log_pdf_type: Either skewT or N.
DEC_obj_fun <- function(DEC1,
                        DEC2,
                        beta_mat, 
                        A_value, 
                        Psi_inv_value,
                        nu_value,
                        subset,
                        log_pdf_type) {
  N <- length(subset$Y)
  Psi_value <- solve(Psi_inv_value)
  output <- numeric(N)
  for (i in 1:N) {
    yi_mat <- subset$Y[[i]]
    Xi_mat <- subset$X[[i]]
    Mi_mat <- Xi_mat %*% beta_mat
    ti <- subset$ti[[i]]
    Sigmai <- Sigma_mat_func(ti = ti,
                             DEC_vec = c(DEC1,DEC2))
    ni <- nrow(yi_mat)
    Ai_mat <- matrix(data = rep(A_value,ni), nrow = ni, byrow = TRUE)
    output[i] <- tryCatch({MV_logdensity(Y = yi_mat,
                                         M = Mi_mat,
                                         A = Ai_mat, 
                                         Sigma = Sigmai, 
                                         Psi = Psi_value, 
                                         nu = nu_value, 
                                         n = nrow(yi_mat), 
                                         p = ncol(yi_mat), 
                                         log_pdf_type = log_pdf_type)},
                          error = function(x) {
                            return(-Inf)
                          })
  }
  output <- sum(output)
  return(output)
}

# update the first DEC parameter
update_DEC1_regular <- function(beta_mat, 
                                A_value, 
                                Psi_inv_value,
                                nu_value,
                                DEC2,
                                cl,
                                subset_list,
                                log_pdf_type,
                                N_cores) {
  N_subset <- length(subset_list)
  
  output <- mclapply(X = 1:N_subset,
                     FUN = function(x) {
                       DEC1_values <- seq(0,1,length.out = 11)
                       DEC1_values[1] <- DEC1_values[1] + 1e-5
                       DEC1_values[11] <- DEC1_values[11] - 1e-5
                       output <- numeric(length(DEC1_values))
                       for (i in 1:length(DEC1_values)) {
                         output[i] <- DEC_obj_fun(DEC1 = DEC1_values[i],
                                                  DEC2 = DEC2,
                                                  beta_mat = beta_mat,
                                                  A_value = A_value,
                                                  Psi_inv_value = Psi_inv_value,
                                                  nu_value = nu_value,
                                                  subset = subset_list[[x]],
                                                  log_pdf_type = log_pdf_type)
                       }
                       return(output)
                     },
                     mc.cores = N_cores)
  
  output <- Reduce("+",output)
  DEC1_values <- seq(0,1,length.out = 11)
  DEC1_values[1] <- DEC1_values[1] + 1e-5
  DEC1_values[11] <- DEC1_values[11] - 1e-5
  
  df_DEC <- data.frame(DEC = DEC1_values,
                       value = output)

  df_DEC <- df_DEC[! is.infinite(df_DEC$value),]
  
  output <- df_DEC[which.max(df_DEC$value),"DEC"]
  
  return(output)
}

# update the second DEC parameter
update_DEC2_regular <- function(beta_mat, 
                                A_value, 
                                Psi_inv_value,
                                nu_value,
                                DEC1,
                                DEC2_lb,
                                DEC2_ub,
                                cl,
                                subset_list,
                                log_pdf_type,
                                N_cores) {
  N_subset <- length(subset_list)
  
  output <- mclapply(X = 1:N_subset,
                     FUN = function(x) {
                       DEC2_values <- seq(DEC2_lb,DEC2_ub,length.out = 11)
                       DEC2_values[1] <- DEC2_values[1] + 1e-5
                       DEC2_values[11] <- DEC2_values[11] - 1e-5
                       output <- numeric(length(DEC2_values))
                       for (i in 1:length(DEC2_values)) {
                         output[i] <- DEC_obj_fun(DEC1 = DEC1,
                                                  DEC2 = DEC2_values[i],
                                                  beta_mat = beta_mat,
                                                  A_value = A_value,
                                                  Psi_inv_value = Psi_inv_value,
                                                  nu_value = nu_value,
                                                  subset = subset_list[[x]],
                                                  log_pdf_type = log_pdf_type)
                       }
                       return(output)
                     },
                     mc.cores = N_cores)
  
  output <- Reduce("+",output)
  
  DEC2_values <- seq(DEC2_lb,DEC2_ub,length.out = 11)
  DEC2_values[1] <- DEC2_values[1] + 1e-5
  DEC2_values[11] <- DEC2_values[11] - 1e-5
  
  df_DEC <- data.frame(DEC = DEC2_values,
                       value = unlist(output))
  
  df_DEC <- df_DEC[! is.infinite(df_DEC$value),]
  
  output <- df_DEC[which.max(df_DEC$value),"DEC"]
  
  return(output)
}

update_DEC1_ADECM <- function(Estep_output) {
  N_subset <- length(Estep_output)
  DEC1_obj <- numeric(length = length(Estep_output[[1]]$DEC1_objective_grid))
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    DEC1_obj <- DEC1_obj + Estep_output_subseti$DEC1_objective_grid
  }

  df_DEC <- data.frame(DEC1_obj = DEC1_obj, # value of objective function
                       DEC1_grid = Estep_output[[1]]$DEC1_grid) # value of DEC(1)/first parameter
  
  df_DEC <- df_DEC[! is.infinite(df_DEC$DEC1_obj),]
  
  # find the DEC(1) value that maximize the objective function
  index_max <- which.max(df_DEC$DEC1_obj)
  output <- df_DEC[index_max,"DEC1_grid"]
  return(output)
}

update_DEC2_ADECM <- function(Estep_output) {
  N_subset <- length(Estep_output)
  DEC2_obj <- numeric(length = length(Estep_output[[1]]$DEC2_objective_grid))
  for (i in 1:N_subset) {
    Estep_output_subseti <- Estep_output[[i]]
    DEC2_obj <- DEC2_obj + Estep_output_subseti$DEC2_objective_grid
  }

  df_DEC <- data.frame(DEC2_obj = DEC2_obj, # value of objective function
                       DEC2_grid = Estep_output[[1]]$DEC2_grid) # value of DEC(2)/second parameter
  
  df_DEC <- df_DEC[! is.infinite(df_DEC$DEC2_obj),]
  
  # find the DEC(1) value that maximize the objective function
  index_max <- which.max(df_DEC$DEC2_obj)
  output <- df_DEC[index_max,"DEC2_grid"]
  return(output)
}
