#' The ECME function for the linear regression model based on the matrix-variant skewT distribution.
#'
#' @param EM_type must be one of "ADECME", "RPECME" and "RNECME".
#' @param N_cores  the number of cores in the computer.0
#' @param N_wait  wait for the number of cores to complete the E steps. Must be less than num_cores.
#' @param Sigma_type a character. Must be either "DEC" or "diagonal". 
#' @param log_pdf_type a character. Must be either "skewT" or "N". 
#' @param Y the list of response matrices.
#' @param X the list of design matrices.
#' @param ti the list of measurement times.
#' @param beta_value the initial value of beta.
#' @param A_value the initial value of A.
#' @param DEC_value the initial value of DEC.
#' @param Psi_value the initial value of Psi.
#' @param nu_value the initial value of nu.
#'
#' @return NULL.
#' @export
ECME <- function(EM_type,
                 N_cores,
                 N_wait,
                 Sigma_type,
                 log_pdf_type,
                 Y,
                 X,
                 ti,
                 beta_value,
                 A_value,
                 DEC_value,
                 Psi_value,
                 nu_value) {
  
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
  
  if (EM_type == "ADECME") {
    EM_type <- "ADECM"
  } else if (EM_type == "RPECME") {
    EM_type <- "RPECM"
  } else {
    EM_type <- "RNECM"
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
  assign("log_pdf_type",log_pdf_type,envir = .GlobalEnv)
  assign("Sigma_type",Sigma_type,envir = .GlobalEnv)
  
  # save the data in the global environment
  assign("Y",Y, envir = .GlobalEnv)
  assign("X",X, envir = .GlobalEnv)
  assign("ti",ti, envir = .GlobalEnv)
  
  p <- ncol(Y[[1]])
  q <- ncol(X[[1]])
  
  if (nrow(beta_value) != q | ncol(beta_value) != p) {
    stop("The dimension of beta_value is incorrect.")
  }
  
  if (length(A_value) != p) {
    stop("The dimension of A_value is incorrect.")
  }
  
  if (! all(DEC_value %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) | length(DEC_value) != 2 ) {
    stop("DEC_value must be two of c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9).")
  }
  
  if (! is.matrix(Psi_value) | ! isSymmetric(Psi_value)) {
    stop("The initial value of Psi is not correct.")
  }
  
  if (nu_value <= 0.0) {
    stop("nu_value must be positive")
  }
  
  # initial values for the ECM function
  assign("beta_value",beta_value, envir = .GlobalEnv)
  assign("A_value",A_value, envir = .GlobalEnv)
  assign("DEC_value",DEC_value, envir = .GlobalEnv)
  assign("Psi_value",Psi_value, envir = .GlobalEnv)
  assign("nu_value",nu_value, envir = .GlobalEnv)
  script_path1 <- system.file("/source_code/AD_ECM.R",package = "STMATREG")
  script_path2 <- system.file("/source_code/AD_ECM_global.R",package = "STMATREG")
  source(script_path1)
  source(script_path2)
  return(output)
}
