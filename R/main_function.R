#' The ECME function for the linear regression model based on the matrix-variant skewT distribution.
#'
#' @param EM_type Must be one of "ADECME", "RPECME" and "RNECME".
#' @param N_cores  The number of cores in the computer.0
#' @param N_wait  Wait for the number of cores to complete the E steps. Must be less than num_cores.
#' @param Sigma_type A character. Must be either "DEC" or "diagonal". 
#' @param log_pdf_type A character. Must be either "skewT" or "N". 
#' @param Y The list of response matrices.
#' @param X The list of design matrices.
#' @param ti The list of measurement times.
#'
#' @return NULL.
#' @export
ECME <- function(EM_type,N_cores,N_wait,Sigma_type,log_pdf_type,Y,X,ti) {
  
  # ensure EM type is one of "ADECME", "RPECME", "RNECME"
  if (! (EM_type %in% c("ADECME", "RPECME", "RNECME"))) {
    stop("EM type must be either ADECME, RPECME or RNECME")
  }
  # for regular ECME algorithm, N_wait = N_cores
  if (EM_type == "RPECME") {
    N_wait <- N_cores
  }
  # for non-parallel ECME algorithm, N_wait = N_cores = 1
  if (EM_type == "RNECME") {
    N_cores <- 1
    N_wait <- 1
  }
  # save setting in the global environment
  assign("EM_type",EM_type, envir = .GlobalEnv)
  assign("N_cores",N_cores, envir = .GlobalEnv)
  assign("prob_all_subsets",0.0, envir = .GlobalEnv) #ensure using asynchronous parallel computer after the first iteration
  assign("N_wait",N_wait, envir = .GlobalEnv)
  assign("DEC2_lb",0.0, envir = .GlobalEnv)
  assign("DEC2_ub",1.0, envir = .GlobalEnv)
  assign("asy_ci",FALSE, envir = .GlobalEnv) #do not calculate the asymptotic confidence interval
  assign("tol",1e-7, envir = .GlobalEnv)
  assign("maxit",1000, envir = .GlobalEnv)
  assign("update_iter",1, envir = .GlobalEnv)
  
  # save the data in the global environment
  assign("Y",Y, envir = .GlobalEnv)
  assign("X",X, envir = .GlobalEnv)
  assign("ti",ti, envir = .GlobalEnv)
  
  p <- ncol(Y[[1]])
  q <- ncol(X[[1]])
  
  # initial values for the ECM function
  assign("beta_value",matrix(rnorm(n = q*p), nrow = q), envir = .GlobalEnv)
  assign("A_value",rnorm(n = p), envir = .GlobalEnv)
  assign("DEC_value",sample(c(0.5,0.6,0.7,0.8,0.9), size = 2, replace = TRUE), envir = .GlobalEnv)
  assign("Psi_dim",p, envir = .GlobalEnv)
  assign("Psi_half",matrix(runif(Psi_dim^2)*2-1, ncol=Psi_dim), envir = .GlobalEnv)
  assign("Psi_value",t(Psi_half) %*% Psi_half, envir = .GlobalEnv)
  assign("nu_value",4 + rgamma(n = 1, shape = 2, rate = 1), envir = .GlobalEnv)
  script_path1 <- system.file("/source_code/AD_ECM.R",package = "STMATREG")
  script_path2 <- system.file("/source_code/AD_ECM_global.R",package = "STMATREG")
  source(script_path1)
  source(script_path2)
  return(output)
}
