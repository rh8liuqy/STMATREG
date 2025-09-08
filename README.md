
# STMATREG

<!-- badges: start -->
<!-- badges: end -->

The `STMATREG` package incorporates a linear regression based on the matrix variant skew-t distribution using the asynchronous parallel ECME, regular parallel ECME, and regular non-parallel ECME algorithms.

## Installation

The core dependency, `async` package, does not support the Windows operating system. Therefore, the `STMATREG` package is only available on UNIX operating systems.

You can install the development version of `STMATREG` like so:

``` r
if(require(devtools)) {
  devtools::install_github(repo = "rh8liuqy/STMATREG",build_vignettes = TRUE)
}
```

## Data Simulation

``` r
library(STMATREG)

# Set the sample size as 200
N <- 200
# Set the true values of parameters
true_beta <- matrix(c(0.5,1.5,-0.5,0.5,1.5,-0.5),nrow = 3)
true_A <- c(2.0,-2.0)
true_DEC <- c(0.9,0.8)
true_Psi <- matrix(c(1.0,-0.5,-0.5,1.0), nrow = 2)
true_nu <- 5
# Begin the data generation
simulated_data <- reg_simulation(N = N,
                                 beta_mat = true_beta,
                                 A_vec = true_A,
                                 DEC_vec = true_DEC,
                                 Psi_mat = true_Psi,
                                 nu = true_nu,
                                 Sigma_type = Sigma_type,
                                 log_pdf_type = log_pdf_type)

Y <- simulated_data$Y # A list of response matrices
X <- simulated_data$X # A list of covariate matrices
ti <- simulated_data$ti # A list of records of time lapses
```

## Maximum Likelihood Estimation via the ECME Algorithm

``` r
ECME(EM_type = "ADECME", # Must be one of "ADECME", "RPECME" and "RNECME".
     N_cores = 8, # The number of cores in the computer.
     N_wait = 7, # Must be no larger than N_cores.
     Sigma_type= "DEC", # Must be one of "DEC" and "diagonal"
     log_pdf_type= "skewT", # Must be "skewT" or "N".
     Y = Y,
     X = X,
     ti = ti,
     beta_value = true_beta,
     A_value = true_A,
     DEC_value = true_DEC,
     Psi_value = true_Psi,
     nu_value = true_nu)
print(output)
```

