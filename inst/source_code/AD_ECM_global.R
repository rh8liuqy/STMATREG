# function argument validation --------------------------------------------

if (!log_pdf_type %in% c("skewT","skewN","N")) {
  stop("log_pdf_type must be one of 'skewT', 'skewN', or 'N'.")
}

if (!Sigma_type %in% c("DEC","diagonal")) {
  stop("Sigma_type must be one of 'DEC' or 'diagonal'.")
}

# manage initial values ---------------------------------------------------

Psi_inv_value <- solve(Psi_value)
N <- length(Y)
diff <- 1e3

if (log_pdf_type != "skewT") {
  nu_value <- 200
}

if (log_pdf_type == "N") {
  A_value <- numeric(length = length(A_value))
}

E_step_time <- numeric(1)
update_beta_time <- numeric(1)
update_A_time <- numeric(1)
update_DEC_time <- numeric(1)
update_Psi_time <- numeric(1)
update_nu_time <- numeric(1)

# subsets creation ------------------------------------------------------

subsets_usage_output <- list()

subset_list <- create_subset(Y = Y,
                             X = X,
                             ti = ti,
                             N_subset = N_cores)


# create cluster ----------------------------------------------------------

cl <- makeCluster(N_cores)

# get the pids of the worker processes
procIds <- clusterCall(cl, Sys.getpid)

# export procIds into the cluster
clusterExport(cl = cl, 
              varlist = c("procIds"),
              envir = environment())

# define a new async cluster function for efficient communication
clusterAsyncCall <- function (cl = NULL, k, pids, fun, ...) {
  cl <- defaultCluster(cl)
  stopifnot(k >= 1 && k <= length(cl))
  stopifnot(all(sapply(cl, function(x) x$host) == "localhost"))
  afun <- function(...)
    tryCatch(allowInterrupts(fun(...)), interrupt = function(e) NULL)
  
  val <- vector("list", length(cl))
  done <- logical(length(cl))
  for (i in seq_along(cl))
    sendCall(cl[[i]], afun, list(...), tag = i)
  for (i in seq_len(k)) {
    d <- recvOneResult(cl)
    val[d$tag] <- list(d$value)
    done[d$tag] <- TRUE
  }
  for (i in seq_along(cl))
    if (! done[i]) {
      tools::pskill(pids[[i]], signal = tools::SIGINT)
      recvResult(cl[[i]])
    }
  checkForRemoteErrors(val)
}
environment(clusterAsyncCall) <- getNamespace("parallel")

# export functions from global environment to clusters
clusterExport(cl = cl,
              varlist = c("DEC_Sigma_mat_func", 
                          "Estep", 
                          "DEC_correlation",
                          "delta_func",
                          "rho_func",
                          "ai_func",
                          "bi_func",
                          "ci_func",
                          "diagonal_Sigma_mat_func",
                          "DEC_obj_fun",
                          "Sigma_mat_func",
                          "MV_logdensity",
                          "MV_logpdf",
                          "create_symmetric_matrix",
                          "MV_logpdf_hess",
                          "subset_list",
                          "EM_type",
                          "log_pdf_type",
                          "DEC2_lb",
                          "DEC2_ub"),
              envir = .GlobalEnv)

# export packages to clusters
clusterEvalQ(cl = cl,
             expr = {
               library(MASS)
               library(tidyverse)
               library(numDeriv)
             })


# matrices for evaluating the loglikelihood -------------------------------

Sigma_value <- vector(mode = "list",
                      length = N)
Sigma_inv_value <- vector(mode = "list",
                          length = N)
A_mat_value <- vector(mode = "list",
                      length = N)
M <- vector(mode = "list",
            length = N)


# initial Sigma_value and Sigma_inv_value for diagonal structure ----------
if (Sigma_type == "diagonal") {
  for (i in 1:N) {
    Sigma_value[[i]] <- diag(rep(1,nrow(X[[i]])))
    Sigma_inv_value[[i]] <- Sigma_value[[i]]  
  }
}

# optimization ------------------------------------------------------------

iter <- 1
condition <- TRUE
convergent <- TRUE
while (condition) {

# record the values of parameters before E step ---------------------------
  
  if (Sigma_type == "DEC") {
    parameter_value1 <- c(as.numeric(beta_value),
                          A_value,
                          DEC_value,
                          as.numeric(Psi_value),
                          nu_value)
  } else {
    parameter_value1 <- c(as.numeric(beta_value),
                          A_value,
                          as.numeric(Psi_value),
                          nu_value)
  }
  
  # export parameter values from current environment to clusters
  clusterExport(cl = cl, 
                varlist = c("beta_value",
                            "A_value",
                            "DEC_value",
                            "Psi_value",
                            "nu_value",
                            "Sigma_type",
                            "Psi_inv_value"),
                envir = .GlobalEnv)
  
  # E step ------------------------------------------------------------------
  
  E_step_t1 <- Sys.time()
  
  # the first E step must use all subsets !
  # draw a random variable deciding whether use all subsets in the E step or not
  u_all_subsets <- runif(n = 1, min = 0, max = 1)
  if (iter == 1 | u_all_subsets <= prob_all_subsets) {
    
    # all values in the Estep_output need to be created
    try(Estep_output <- clusterAsyncCall(cl = cl,
                                         k = N_cores,
                                         pids = procIds,
                                         fun = function(pars) {
                                           index <- procIds == Sys.getpid()
                                           subset_index <- pars[index]
                                           subseti <- subset_list[[subset_index]]
                                           Yi <- subseti$Y
                                           Xi <- subseti$X
                                           tii <- subseti$ti
                                           
                                           Ni <- length(Yi)
                                           Sigmai_value <- vector(mode = "list",
                                                                  length = Ni)
                                           Sigmai_inv_value <- vector(mode = "list",
                                                                      length = Ni)
                                           Ai_mat_value <- vector(mode = "list",
                                                                  length = Ni)
                                           
                                           for (i in 1:Ni) {
                                             ni <- nrow(Yi[[i]])
                                             if (Sigma_type == "DEC") {
                                               Sigmai_value[[i]] <- DEC_Sigma_mat_func(ti = tii[[i]],
                                                                                       DEC_vec = DEC_value)
                                               Sigmai_inv_value[[i]] <- solve(Sigmai_value[[i]])
                                             } else {
                                               Sigmai_value[[i]] <- diagonal_Sigma_mat_func(ti = tii[[i]])
                                               Sigmai_inv_value[[i]] <- Sigmai_value[[i]]
                                             }
                                             
                                             Ai_mat_value[[i]] <- matrix(data = rep(A_value,ni),
                                                                         nrow = ni,
                                                                         byrow = TRUE)
                                           }
                                           output <- Estep(Y_list = Yi,
                                                           X_list = Xi,
                                                           ti_list = tii,
                                                           beta_mat = beta_value,
                                                           A_mat_list = Ai_mat_value,
                                                           Sigma_inv_list = Sigmai_inv_value,
                                                           Psi_inv = Psi_inv_value,
                                                           nu = nu_value,
                                                           log_pdf_type = log_pdf_type,
                                                           EM_type = EM_type,
                                                           Sigma_type = Sigma_type,
                                                           DEC2_lb = DEC2_lb,
                                                           DEC2_ub = DEC2_ub,
                                                           DEC_value = DEC_value)
                                           return(output)
                                         },
                                         1:N_cores))
    
    # monitor the usage of subsets
    subsets_usage_value <- data.frame(t(rep(1,N_cores)))
    colnames(subsets_usage_value) <- paste0("subset",1:N_cores)
    subsets_usage_output[[iter]] <- subsets_usage_value
    
  } else {
    # wait until N_wait cores complete the computation
    try(Estep_subset_output <- clusterAsyncCall(cl = cl,
                                                k = N_wait,
                                                pids = procIds,
                                                fun = function(pars) {
                                                  index <- procIds == Sys.getpid()
                                                  subset_index <- pars[index]
                                                  subseti <- subset_list[[subset_index]]
                                                  Yi <- subseti$Y
                                                  Xi <- subseti$X
                                                  tii <- subseti$ti
                                                  
                                                  Ni <- length(Yi)
                                                  Sigmai_value <- vector(mode = "list",
                                                                         length = Ni)
                                                  Sigmai_inv_value <- vector(mode = "list",
                                                                             length = Ni)
                                                  Ai_mat_value <- vector(mode = "list",
                                                                         length = Ni)
                                                  
                                                  for (i in 1:Ni) {
                                                    ni <- nrow(Yi[[i]])
                                                    if (Sigma_type == "DEC") {
                                                      Sigmai_value[[i]] <- DEC_Sigma_mat_func(ti = tii[[i]],
                                                                                              DEC_vec = DEC_value)
                                                      Sigmai_inv_value[[i]] <- solve(Sigmai_value[[i]])
                                                    } else {
                                                      Sigmai_value[[i]] <- diagonal_Sigma_mat_func(ti = tii[[i]])
                                                      Sigmai_inv_value[[i]] <- Sigmai_value[[i]]
                                                    }
                                                    
                                                    Ai_mat_value[[i]] <- matrix(data = rep(A_value,ni),
                                                                                nrow = ni,
                                                                                byrow = TRUE)
                                                  }
                                                  output <- Estep(Y_list = Yi,
                                                                  X_list = Xi,
                                                                  ti_list = tii,
                                                                  beta_mat = beta_value,
                                                                  A_mat_list = Ai_mat_value,
                                                                  Sigma_inv_list = Sigmai_inv_value,
                                                                  Psi_inv = Psi_inv_value,
                                                                  nu = nu_value,
                                                                  log_pdf_type = log_pdf_type,
                                                                  EM_type = EM_type,
                                                                  Sigma_type = Sigma_type,
                                                                  DEC2_lb = DEC2_lb,
                                                                  DEC2_ub = DEC2_ub,
                                                                  DEC_value = DEC_value)
                                                  return(output)
                                                },
                                                1:N_cores))
    
    # monitor the usage of subsets
    subsets_usage_value <- numeric(N_cores)
    for (i in 1:N_cores) {
      valuei <- Estep_subset_output[[i]]
      if (! is.null(valuei)) {
        Estep_output[[i]] <- valuei
        subsets_usage_value[i] <- 1
      }
    }
    subsets_usage_value <- data.frame(t(subsets_usage_value))
    colnames(subsets_usage_value) <- paste0("subset",1:N_cores)
    subsets_usage_output[[iter]] <- subsets_usage_value
  }
  
  E_step_t2 <- Sys.time()
  E_step_time <- E_step_time + E_step_t2 - E_step_t1
  
  # M step ------------------------------------------------------------------
  
  
  # update_beta -------------------------------------------------------------
  
  update_beta_t1 <- Sys.time()
  
  beta_value <- update_beta(Estep_output = Estep_output)
  
  update_beta_t2 <- Sys.time()
  
  update_beta_time <- update_beta_time + update_beta_t2 - update_beta_t1
  
  # update nu -------------------------------------------------------------
  
  update_nu_t1 <- Sys.time()
  
  if (log_pdf_type == "skewT") {
    nu_value <- update_nu(Estep_output = Estep_output)
  } else {
    # only skew T distribution need the degree of freedom
    # for skew N and normal distribution, use nu = 200
    # for numerical stability we can not pick bigger number for nu.
    nu_value <- 200
  }
  
  update_nu_t2 <- Sys.time()
  update_nu_time <- update_nu_time + update_nu_t2 - update_nu_t1
  
  # update A -------------------------------------------------------------
  
  update_A_t1 <- Sys.time()
  
  if (log_pdf_type %in% c("skewT","skewN")) {
    if (EM_type == "ADECM") {
      A_value <- update_A_ADECM(Estep_output = Estep_output)
    } else {
      A_value <- update_A_regular(Estep_output = Estep_output,
                                  beta_mat = beta_value,
                                  subset_list = subset_list)
    }
  } else {
    A_value <- numeric(length = length(A_value))
  }
  
  update_A_t2 <- Sys.time()
  update_A_time <- update_A_time + update_A_t2 - update_A_t1
  
  # update Psi -------------------------------------------------------------
  
  update_Psi_t1 <- Sys.time()
  if (EM_type == "ADECM") {
    Psi_value <- update_Psi_ADECM(Estep_output = Estep_output)
  } else {
    Psi_value <- update_Psi_regular(Estep_output = Estep_output,
                                    beta_mat = beta_value,
                                    A_value = A_value,
                                    subset_list = subset_list)
  }

  Psi_inv_value <- solve(Psi_value)
  
  update_Psi_t2 <- Sys.time()
  update_Psi_time <- update_Psi_time + update_Psi_t2 - update_Psi_t1
  
   # update DEC -------------------------------------------------------------
  
  update_DEC_t1 <- Sys.time()
  
  if (Sigma_type == "DEC") {
    if (EM_type == "ADECM") {
      DEC_value[1] <- update_DEC1_ADECM(Estep_output = Estep_output)
      DEC_value[2] <- update_DEC2_ADECM(Estep_output = Estep_output)
    } else {
      # update DEC1
      DEC_value[1] <- update_DEC1_regular(beta_mat = beta_value,
                                          A_value = A_value,
                                          Psi_inv_value = Psi_inv_value,
                                          nu_value = nu_value,
                                          DEC2 = DEC_value[2],
                                          cl = cl,
                                          subset_list = subset_list,
                                          log_pdf_type = log_pdf_type,
                                          N_cores = N_cores)

      # update DEC2
      DEC_value[2] <- update_DEC2_regular(beta_mat = beta_value,
                                          A_value = A_value,
                                          Psi_inv_value = Psi_inv_value,
                                          nu_value = nu_value,
                                          DEC1 = DEC_value[1],
                                          DEC2_lb = DEC2_lb,
                                          DEC2_ub = DEC2_ub,
                                          cl = cl,
                                          subset_list = subset_list,
                                          log_pdf_type = log_pdf_type,
                                          N_cores = N_cores)
    }
  }
  
  update_DEC_t2 <- Sys.time()
  
  update_DEC_time <- update_DEC_time + update_DEC_t2 - update_DEC_t1


# record the values of parameters after E step ----------------------------
    
  if (Sigma_type == "DEC") {
    parameter_value2 <- c(as.numeric(beta_value),
                          A_value,
                          DEC_value,
                          as.numeric(Psi_value),
                          nu_value)
  } else {
    parameter_value2 <- c(as.numeric(beta_value),
                          A_value,
                          as.numeric(Psi_value),
                          nu_value)
  }
  
  # diff is the absolute difference between point estimations before/after one ECM iteration
  diff <- max(abs(parameter_value2 - parameter_value1))
  
  # print the log
  if (iter %% update_iter == 0) {
    print(paste0("iteration: ", iter, " || ", "max.abs.diff: ", diff))
  }
  
  # end the updating here
  if (diff < tol) {
    condition <- FALSE
  }

  iter <- iter + 1
  if (iter > maxit) {
    condition <- FALSE
    convergent <- FALSE
    message("reach maximum number of iterations.")
  }
  
}

# evaluation of the loglikelihood after ECM -------------------------------

# updating matrices for evaluating the loglikelihood
m_output <- mclapply(1:N,
                     FUN = function(i){
                       ni <- nrow(Y[[i]])
                       if (Sigma_type == "DEC") {
                         Sigma_value <- DEC_Sigma_mat_func(ti = ti[[i]],
                                                           DEC_vec = DEC_value)
                         Sigma_inv_value <- solve(Sigma_value)
                       }

                       A_mat_value <- matrix(data = rep(A_value,ni),
                                             nrow = ni,
                                             byrow = TRUE)
                       M <- X[[i]]%*%beta_value
                       if (Sigma_type == "DEC") {
                         output <- list(Sigma_value = Sigma_value,
                                        Sigma_inv_value = Sigma_inv_value,
                                        A_mat_value = A_mat_value,
                                        M = M)
                       } else {
                         output <- list(A_mat_value = A_mat_value,
                                        M = M)
                       }
                       return(output)
                     },
                     mc.cores = N_cores)

if (Sigma_type == "DEC") {
  for (i in 1:N) {
    Sigma_value[[i]] <- m_output[[i]]$Sigma_value
    Sigma_inv_value[[i]] <- m_output[[i]]$Sigma_inv_value
    A_mat_value[[i]] <-  m_output[[i]]$A_mat_value
    M[[i]] <- m_output[[i]]$M
  }
} else {
  for (i in 1:N) {
    A_mat_value[[i]] <-  m_output[[i]]$A_mat_value
    M[[i]] <- m_output[[i]]$M
  }
}

ll2 <- MST_logpdf_c(Y = Y,
                    M = M,
                    A = A_mat_value,
                    Sigma = Sigma_value,
                    Psi = Psi_value,
                    nu = nu_value,
                    log_pdf_type = log_pdf_type,
                    Sigma_type = Sigma_type)

# Organize output ---------------------------------------------------------

if (log_pdf_type == "skewT" & Sigma_type == "DEC") {
  parameter <- c(as.vector(beta_value),
                 A_value,
                 DEC_value,
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)],
                 nu_value)
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_nu_time = update_nu_time,
                     update_A_time = update_A_time,
                     update_Psi_time = update_Psi_time,
                     update_DEC_time = update_DEC_time)
} else if (log_pdf_type == "skewN" & Sigma_type == "DEC") {
  parameter <- c(as.vector(beta_value),
                 A_value,
                 DEC_value,
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)])
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_A_time = update_A_time,
                     update_Psi_time = update_Psi_time,
                     update_DEC_time = update_DEC_time)
} else if (log_pdf_type == "N" & Sigma_type == "DEC") {
  parameter <- c(as.vector(beta_value),
                 DEC_value,
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)])
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_Psi_time = update_Psi_time,
                     update_DEC_time = update_DEC_time)
} else if (log_pdf_type == "skewT" & Sigma_type == "diagonal") {
  parameter <- c(as.vector(beta_value),
                 A_value,
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)],
                 nu_value)
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_nu_time = update_nu_time,
                     update_A_time = update_A_time,
                     update_Psi_time = update_Psi_time)
} else if (log_pdf_type == "skewN" & Sigma_type == "diagonal") {
  parameter <- c(as.vector(beta_value),
                 A_value,
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)])
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_A_time = update_A_time,
                     update_Psi_time = update_Psi_time)
} else if (log_pdf_type == "N" & Sigma_type == "diagonal") {
  parameter <- c(as.vector(beta_value),
                 Psi_value[upper.tri(x = Psi_value,diag = TRUE)])
  time_value <- list(E_step_time = E_step_time,
                     update_beta_time = update_beta_time,
                     update_Psi_time = update_Psi_time)
}


# observed information matrix ---------------------------------------------  

if (asy_ci) {
  obs_info_mat <- clusterApply(cl,
                               1:N_cores,
                               function(subset_index) {
                                 subseti <- subset_list[[subset_index]]
                                 Yi <- subseti$Y
                                 Xi <- subseti$X
                                 tii <- subseti$ti
                                 output <- numDeriv::hessian(func = MV_logpdf_hess,
                                                             x = parameter,
                                                             method = "Richardson",
                                                             method.args=list(eps=1e-4, d=1e-4, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE),
                                                             Y = Yi,
                                                             X = Xi,
                                                             ti = tii,
                                                             Sigma_type = Sigma_type,
                                                             log_pdf_type = log_pdf_type)
                                 output <- -output
                                 return(output)
                               })
  
  obs_info_mat <- Reduce("+",obs_info_mat)
  obs_info_mat_inv <- tryCatch(expr = {
    solve(obs_info_mat)
  },
  error = function(e) {
    return(solve(obs_info_mat + diag(rep(1e-5,nrow(obs_info_mat)))))
  })
  parameter_ub <- parameter + qnorm(1 - 0.05/2) * sqrt(abs(diag(obs_info_mat_inv)))
  parameter_lb <- parameter - qnorm(1 - 0.05/2) * sqrt(abs(diag(obs_info_mat_inv)))
} else {
  parameter_ub <- NA
  parameter_lb <- NA
}

# Residual calculation ----------------------------------------------------

residual <- pmap(list(Y = Y, M = M), function(Y, M) {return(Y - M)})

# RMSD calculation --------------------------------------------------------

RMSD1 <- sum(sapply(residual,function(x){sum(x^2)}))
RMSD2 <- sum(sapply(residual,function(x){length(x)}))
RMSD <- sqrt(RMSD1/RMSD2)

stopCluster(cl)

if (log_pdf_type == "skewT" & Sigma_type == "DEC") {
  k <- length(beta_value)+
    length(A_value)+
    length(DEC_value)+
    length(Psi_value)+
    length(nu_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 A = A_value,
                 DEC = DEC_value,
                 Psi = Psi_value,
                 nu = nu_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
} else if (log_pdf_type == "skewN" & Sigma_type == "DEC") {
  k <- length(beta_value)+
    length(A_value)+
    length(DEC_value)+
    length(Psi_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 A = A_value,
                 DEC = DEC_value,
                 Psi = Psi_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
} else if (log_pdf_type == "N" & Sigma_type == "DEC") {
  k <- length(beta_value)+
    length(DEC_value)+
    length(Psi_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 DEC = DEC_value,
                 Psi = Psi_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
} else if (log_pdf_type == "skewT" & Sigma_type == "diagonal") {
  k <- length(beta_value)+
    length(A_value)+
    length(Psi_value)+
    length(nu_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 A = A_value,
                 Psi = Psi_value,
                 nu = nu_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
} else if (log_pdf_type == "skewN" & Sigma_type == "diagonal") {
  k <- length(beta_value)+
    length(A_value)+
    length(Psi_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 A = A_value,
                 Psi = Psi_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
} else if (log_pdf_type == "N" & Sigma_type == "diagonal") {
  k <- length(beta_value)+
    length(Psi_value)
  AIC_info <- 2*k-2*ll2
  sample_size <- sum(sapply(Y, length))
  BIC_info <- k*log(sample_size)-2*ll2
  
  output <- list(beta = beta_value, 
                 Psi = Psi_value,
                 point_estimation = parameter,
                 point_estimation_ub = parameter_ub,
                 point_estimation_lb = parameter_lb,
                 residual = residual,
                 RMSD = RMSD,
                 subsets_usage = list_rbind(subsets_usage_output),
                 convergent = convergent,
                 AIC = AIC_info,
                 BIC = BIC_info,
                 time = time_value,
                 num_of_iterations = iter - 1)
}

# remove the CI from the output if the asy_ci option is off.
if (! asy_ci) {
  output$point_estimation_ub <- NULL
  output$point_estimation_lb <- NULL
}

print("beta")
print(output$beta)
print("A")
print(output$A)
print("Psi")
print(output$Psi)
print("nu")
print(output$nu)
print("DEC")
print(output$DEC)
