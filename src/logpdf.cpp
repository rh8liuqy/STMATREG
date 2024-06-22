#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Rcpp function to calculate the log-likelihood
// [[Rcpp::export]]
double MVskewt_c(const arma::mat& Y, const arma::mat& M, const arma::mat& A,
                 const arma::mat& Sigma, const arma::mat& Psi, double nu, int n, int p,
                 const StringVector& log_pdf_type,
                 const StringVector& Sigma_type) {
  
  // import besselK function from R
  Function besselK("besselK");
  double logDens;
  if (log_pdf_type[0] != "N") {
    double lambda = -(nu + n * p) / 2.0;
    
    // Calculate Delta
    arma::mat YMinusM = Y - M;
    arma::mat SigmaInv;
    if (Sigma_type[0] == "DEC") {
      SigmaInv = solve(Sigma, eye(n, n));
    } else {
      SigmaInv = Sigma;
    }
    arma::mat PsiInv = solve(Psi, eye(p, p));
    double Delta = trace(SigmaInv * YMinusM * PsiInv * YMinusM.t()) + nu;
    
    // Calculate rho
    double rho = trace(SigmaInv * A * PsiInv * A.t());
    
    // Calculate Bessel function
    double sqrtDeltaRho = sqrt(Delta * rho);
    
    NumericVector BessX_vec = besselK(sqrtDeltaRho, lambda, false);
    double BessX = BessX_vec[0];
    
    // For numerical stability
    if (BessX < 1e-100) {
      BessX = 1e-100;
    }
    
    // Calculate the log-likelihood
    logDens = log(2.0) + (nu / 2.0) * log(nu / 2.0) -
      log(2.0 * M_PI) * (n * p / 2.0) - log(det(Sigma)) * (p / 2.0) -
      log(det(Psi)) * (n / 2.0) - lgamma(0.5*nu) + 
      (lambda / 2.0) * (log(Delta) - log(rho)) +
      log(BessX) + trace(SigmaInv * YMinusM * PsiInv * A.t());
  }
  
  else {
    // Calculate Delta
    arma::mat YMinusM = Y - M;
    arma::mat SigmaInv;
    if (Sigma_type[0] == "DEC") {
      SigmaInv = solve(Sigma, eye(n, n));
    } else {
      SigmaInv = Sigma;
    }
    arma::mat PsiInv = solve(Psi, eye(p, p));
    
    double Delta = trace(SigmaInv * YMinusM * PsiInv * YMinusM.t());
    logDens = - log(2.0 * M_PI) * (n * p / 2.0) - log(det(Sigma)) * (p / 2.0)
      - log(det(Psi)) * (n / 2.0) - 0.5 * Delta;
  }
  return logDens;
}

//' Rcpp function to calculate the log-likelihood for all observations
//' @param Y The list of response matrices.
//' @param M The list of location matrices.
//' @param A The list of skewness matrices.
//' @param Sigma The list of scale matrices.
//' @param Psi The column covariance matrix.
//' @param nu The degree of freedom
//' @param log_pdf_type Can be either "skewT" or "N".
//' @param Sigma_type Can be either "DEC" or "diagonal".
//' @return The evaluated log-likelihood.
//' @export
// [[Rcpp::export]]
double MST_logpdf_c(List Y, List M, List A, List Sigma, const arma::mat& Psi, 
                    double nu, 
                    const StringVector& log_pdf_type, 
                    const StringVector& Sigma_type) {
  int N = Y.size();
  double logLikelihood = 0.0;
  
  for (int i = 0; i < N; ++i) {
    arma::mat Yi = as<arma::mat>(Y[i]);
    arma::mat Mi = as<arma::mat>(M[i]);
    arma::mat Ai = as<arma::mat>(A[i]);
    arma::mat Sigmai = as<arma::mat>(Sigma[i]);
    int ni = Yi.n_rows;
    int pi = Yi.n_cols;
    
    // Call the MVskewt Rcpp function
    double mvskewtLogLikelihood = MVskewt_c(Yi, Mi, Ai, Sigmai, Psi, 
                                            nu, ni, pi, 
                                            log_pdf_type, Sigma_type);
    logLikelihood += mvskewtLogLikelihood;
  }
  
  return logLikelihood;
}

arma::mat DEC_correlation_c(NumericVector x, double rho, double theta) {
  if (rho < 0 || rho >= 1) {
    stop("0 <= rho < 1 must hold!");
  }
  if (theta < 0) {
    stop("theta must be positive");
  }
  
  int n = x.size();
  arma::mat output(n, n, fill::none);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      output(i, j) = pow(rho,pow(std::abs(x[i] - x[j]),theta));
    }
  }
  
  return output;
}

// Function to calculate the Sigma matrix
arma::mat Sigma_mat_func_c(NumericVector ti, NumericVector DEC_vec) {
  arma::mat output = DEC_correlation_c(ti, DEC_vec[0], DEC_vec[1]);
  return output;
}

// Function to calculate DEC1_obj_fun
// [[Rcpp::export]]
double DEC1_obj_fun_c(double DEC1, List Estep_output, const NumericMatrix& beta_mat,
                      NumericVector A_value, const arma::mat& Psi_inv_value, double DEC2) {
  double p1 = 0.0;
  double p2 = 0.0;
  double p3 = 0.0;
  double p4 = 0.0;
  double p5 = 0.0;
  int N_subset = Estep_output.size();
  int p = beta_mat.ncol();
  arma::mat beta_mat_value = as<arma::mat>(beta_mat);
  
  List Yi;
  for (int i = 0; i < N_subset; ++i) {
    List Estep_output_subseti = Estep_output[i];
    Yi = Estep_output_subseti["Y"];
    int Ni = Yi.size();
    for (int j = 0; j < Ni; ++j) {
      List Estep_output_subseti_Y = Estep_output_subseti["Y"];
      List Estep_output_subseti_X = Estep_output_subseti["X"];
      List Estep_output_subseti_ti = Estep_output_subseti["ti"];
      List Estep_output_subseti_ni = Estep_output_subseti["ni"];
      List Estep_output_subseti_a = Estep_output_subseti["a"];
      List Estep_output_subseti_b = Estep_output_subseti["b"];
      
      arma::mat Yj = as<arma::mat>(Estep_output_subseti_Y[j]);
      arma::mat Xj = as<arma::mat>(Estep_output_subseti_X[j]);
      NumericVector tj = Estep_output_subseti_ti[j];
      arma::mat Mj = Xj * beta_mat_value;
      arma::mat yjstar = Yj - Mj;
      arma::mat tyjstar = trans(yjstar);
      int nj = as<int>(Estep_output_subseti_ni[j]);
      arma::mat Amatj(nj, p);
      for (int row = 0; row < nj; ++row) {
        Amatj.row(row) = as<rowvec>(A_value);
      }
      arma::mat tAmatj = trans(Amatj);
      double aj = as<double>(Estep_output_subseti_a[j]);
      double bj = as<double>(Estep_output_subseti_b[j]);
      
      // Calculate Sigmaj using the provided Sigma_mat_func
      arma::mat Sigmaj = Sigma_mat_func_c(tj, NumericVector::create(DEC1, DEC2));
      arma::mat Sigmaj_inv = solve(Sigmaj,eye(nj, nj));
      arma::mat Sigmaj_inv_yjstar_Psi_inv = Sigmaj_inv * yjstar * Psi_inv_value;
      arma::mat Sigmaj_inv_Amatj_Psi_inv = Sigmaj_inv * Amatj * Psi_inv_value;
      
      p1 -= 0.5 * p * log(det(Sigmaj));
      p2 += 0.5 * trace(Sigmaj_inv_yjstar_Psi_inv * tAmatj);
      p3 += 0.5 * trace(Sigmaj_inv_Amatj_Psi_inv * tyjstar);
      p4 -= 0.5 * bj * trace(Sigmaj_inv_yjstar_Psi_inv * tyjstar);
      p5 -= 0.5 * aj * trace(Sigmaj_inv_Amatj_Psi_inv * tAmatj);
    }
  }
  
  double output = p1 + p2 + p3 + p4 + p5;
  return output;
}

// Function to calculate DEC2_obj_fun
// [[Rcpp::export]]
double DEC2_obj_fun_c(double DEC2, List Estep_output, const NumericMatrix& beta_mat,
                      NumericVector A_value, const arma::mat& Psi_inv_value, double DEC1) {
  return DEC1_obj_fun_c(DEC1, Estep_output, beta_mat, A_value, Psi_inv_value, DEC2);
}
