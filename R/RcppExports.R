# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

MVskewt_c <- function(Y, M, A, Sigma, Psi, nu, n, p, log_pdf_type, Sigma_type) {
    .Call(`_STMATREG_MVskewt_c`, Y, M, A, Sigma, Psi, nu, n, p, log_pdf_type, Sigma_type)
}

#' Rcpp function to calculate the log-likelihood for all observations
#' @param Y The list of response matrices.
#' @param M The list of location matrices.
#' @param A The list of skewness matrices.
#' @param Sigma The list of scale matrices.
#' @param Psi The column covariance matrix.
#' @param nu The degree of freedom
#' @param log_pdf_type Can be either "skewT" or "N".
#' @param Sigma_type Can be either "DEC" or "diagonal".
#' @return The evaluated log-likelihood.
#' @export
MST_logpdf_c <- function(Y, M, A, Sigma, Psi, nu, log_pdf_type, Sigma_type) {
    .Call(`_STMATREG_MST_logpdf_c`, Y, M, A, Sigma, Psi, nu, log_pdf_type, Sigma_type)
}

DEC1_obj_fun_c <- function(DEC1, Estep_output, beta_mat, A_value, Psi_inv_value, DEC2) {
    .Call(`_STMATREG_DEC1_obj_fun_c`, DEC1, Estep_output, beta_mat, A_value, Psi_inv_value, DEC2)
}

DEC2_obj_fun_c <- function(DEC2, Estep_output, beta_mat, A_value, Psi_inv_value, DEC1) {
    .Call(`_STMATREG_DEC2_obj_fun_c`, DEC2, Estep_output, beta_mat, A_value, Psi_inv_value, DEC1)
}

