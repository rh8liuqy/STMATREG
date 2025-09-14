source("./AD_ECM.R")

log_pdf_type <- "skewT"
Sigma_type <- "DEC"
N <- 100000
r1 <- 0.00 ## choose a very small r1

# the number of simulation studies in scheme1

n_simu <- 50

if (! log_pdf_type %in% c("skewT","skewN","N")) {
  stop("log_pdf_type must be one of skewT, skewN or N.")
}

if (! Sigma_type %in% c("DEC","diagonal")) {
  stop("Sigma_type must be DEC or diagonal.")
}

if (! EM_type %in% c("ADECM", "RPECM", "RNECM")) {
  stop("EM_type must be ADECM, RPECM or RNECM")
}

# detect the number of cores

if (EM_type == "RNECM") {
  num_cores <- 1
  r1 <- 1.0
  r2 <- 1.0
} else if (EM_type == "RPECM"){
  if (as.character(Sys.info()[7]) == "kevin_liu") {
    num_cores <- 8
  } else {
    num_cores <- 64
  }
  
  r1 <- 1.0
  r2 <- 1.0
} else {
  if (as.character(Sys.info()[7]) == "kevin_liu") {
    num_cores <- 8
  } else {
    num_cores <- 64
  }
}

# true value of parameters ------------------------------------------------

true_beta <- matrix(c(0.5,1.5,-0.5,0.5,1.5,-0.5),nrow = 3)
if (log_pdf_type == "skewT") {
  true_A <- c(2.0,-2.0)
} else {
  true_A <- c(0.0,0.0)
}
true_DEC <- c(0.9,0.8)
true_Psi <- matrix(c(1.0,-0.5,-0.5,1.0), nrow = 2)
if (log_pdf_type == "skewT") {
  true_nu <- 5 
} else {
  true_nu <- 200
}

print(paste0("num_cores: ", num_cores))

simu_output <- vector("list",n_simu)

j <- 1
while (j <= n_simu) {
  
  set.seed(j)
  
  print(paste0("Simulation Iteration: ", j))
  
  # point estimation part  --------------------------------------------------
  
  simulated_data <- reg_simulation(N = N,
                                   beta_mat = true_beta,
                                   A_vec = true_A,
                                   DEC_vec = true_DEC,
                                   Psi_mat = true_Psi,
                                   nu = true_nu,
                                   Sigma_type = Sigma_type,
                                   log_pdf_type = log_pdf_type)
  
  Y <- simulated_data$Y
  X <- simulated_data$X
  ti <- simulated_data$ti
  
  while (TRUE) {
    
    # initial values for the ECM function
    
    beta_value <- matrix(rnorm(n = 3*2), nrow = 3)
    A_value <- rnorm(n = 2)
    DEC_value <- sample(c(0.5,0.6,0.7,0.8,0.9), size = 2, replace = TRUE)
    Psi_dim <- 2
    Psi_half <- matrix(runif(Psi_dim^2)*2-1, ncol=Psi_dim)
    Psi_value <- t(Psi_half) %*% Psi_half
    nu_value <- 4 + rgamma(n = 1, shape = 2, rate = 1)
    
    N_cores <- num_cores
    N_wait <- floor(num_cores*r2)
    prob_all_subsets <- r1
    DEC2_lb <- 0
    DEC2_ub <- 1
    asy_ci <- FALSE
    tol <- 1e-7
    maxit <- 1000
    update_iter <- 1
    
    source("./AD_ECM_global.R")
    
    # ensure the success of point estimation.
    
    if (all(output$DEC == true_DEC) & output$convergent) {
      simu_output[[j]] <- output
      j <- j + 1
      break
    } 
  }
}
