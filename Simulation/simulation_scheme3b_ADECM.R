EM_type <- "ADECM"

log_pdf_type <- "skewT"

Sigma_type <- "DEC"

N <- 100000

r2 <- 7/8

source("./simulation_scheme3.R")

saveRDS(simu_output,"simulation_scheme3b.RDS")
