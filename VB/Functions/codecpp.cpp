#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;
const double log2pi = std::log(2.0 * M_PI);

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// [[Rcpp::export]]
arma::vec mvrnormArma(arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return mu + arma::trans(Y * arma::chol(Sigma));
}

// [[Rcpp::export]]
double log_dmvnorm_cpp(arma::vec data, arma::vec m, arma::mat Sigma){
  int xdim = data.size();
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  double constants = -(xdim/2) * log2pi;
  arma::vec z = rooti * ( data - m) ;  
  
  return (constants - 0.5 * arma::sum(z%z) + rootisum);   
}

///////// MATRIX PRODUCT FUNCTIONS


double XiSigmapsiXi(arma::vec &X_i, int X_y_index_i, 
                    int X_s_index_i, arma::mat &Sigma, 
                    int Y, int centers, int ncov){
  
  
  int p = X_i.size();
  int numOfTerms = 1;
  if(centers > 0){
    numOfTerms += 1;
  }
  
  arma::vec out = arma::zeros(numOfTerms + ncov);
  
  // time
  out[0] = Sigma(X_y_index_i - 1, X_y_index_i - 1);
  if(centers > 0){
    out[0] += Sigma(Y + X_s_index_i - 1, X_y_index_i - 1);
  }
  for(int k = 0; k < ncov; k++){
    out[0] += X_i[Y + centers + k] * Sigma(Y + centers + k, X_y_index_i - 1);
  }
  
  // space
  if(centers > 0){
    out[1] = Sigma(X_y_index_i - 1, Y + X_s_index_i - 1) + 
      Sigma(Y + X_s_index_i - 1, Y + X_s_index_i - 1);
    for(int k = 0; k < ncov; k++){
      out[1] += X_i[Y + centers + k] * Sigma(Y + centers + k, Y + X_s_index_i - 1);
    }
  }
  
  // covariates
  for(int k = 0; k < ncov; k++){
    out[k + numOfTerms] = Sigma(X_y_index_i - 1, Y + centers + k);
    if(centers > 0){
      out[k + numOfTerms] += Sigma(Y + X_s_index_i - 1, Y + centers + k);
    }
    for(int k2 = 0; k2 < ncov; k2++){
      out[k + numOfTerms] += X_i[Y + centers + k2] * Sigma(Y + centers + k2, Y + centers + k);
    }
  }
  
  double result2 = 0;
  result2 += out[0];
  if(centers > 0){
    result2 += out[1];
  }
  
  for(int k = 0; k < ncov; k++){
    result2 += out[k + numOfTerms] * X_i[Y + centers + k];
  }
  
  return(result2);
}

double XiSigmapsiXi2(arma::rowvec &X_i, int X_y_index_i, 
                     int X_s_index_i, arma::mat &Sigma, 
                     int Y, int centers, int ncov){
  
  
  int p = X_i.size();
  int numOfTerms = 1;
  if(centers > 0){
    numOfTerms += 1;
  }
  
  arma::vec out = arma::zeros(numOfTerms + ncov);
  
  // time
  out[0] = Sigma(X_y_index_i - 1, X_y_index_i - 1);
  if(centers > 0){
    out[0] += Sigma(Y + X_s_index_i - 1, X_y_index_i - 1);
  }
  for(int k = 0; k < ncov; k++){
    out[0] += X_i[Y + k] * Sigma(Y + centers + k, X_y_index_i - 1);
  }
  
  // space
  if(centers > 0){
    out[1] = Sigma(X_y_index_i - 1, Y + X_s_index_i - 1) + 
      Sigma(Y + X_s_index_i - 1, Y + X_s_index_i - 1);
    for(int k = 0; k < ncov; k++){
      out[1] += X_i[Y + k] * Sigma(Y + centers + k, Y + X_s_index_i - 1);
    }
  }
  
  // covariates
  for(int k = 0; k < ncov; k++){
    out[k + numOfTerms] = Sigma(X_y_index_i - 1, Y + centers + k);
    if(centers > 0){
      out[k + numOfTerms] += Sigma(Y + X_s_index_i - 1, Y + centers + k);
    }
    for(int k2 = 0; k2 < ncov; k2++){
      out[k + numOfTerms] += X_i[Y + k2] * Sigma(Y + centers + k2, Y + centers + k);
    }
  }
  
  double result2 = 0;
  result2 += out[0];
  if(centers > 0){
    result2 += out[1];
  }
  
  for(int k = 0; k < ncov; k++){
    result2 += out[k + numOfTerms] * X_i[Y + k];
  }
  
  return(result2);
}

// [[Rcpp::export]]
arma::vec XSigmapsiX(arma::mat &X, IntegerVector X_y_index, 
                     IntegerVector X_s_index, arma::mat &Sigma, 
                     int Y, int centers, int ncov){
  
  
  arma::vec out = arma::zeros(X_y_index.size());
  
  for(int i = 0; i < X_y_index.size(); i++){
    
    arma::rowvec X_l = X.row(i);
    out[i] =  XiSigmapsiXi2(X_l, X_y_index[i],  X_s_index[i],
                            Sigma, Y, centers, ncov);
    
  }
  
  return(out);
}

double XiSigmapXi(arma::rowvec X_i, int X_y_index_p, 
                  arma::mat &Sigma, 
                  int p_intercepts, int ncov){
  
  
  int p = X_i.size();
  int numOfTerms = 1;
  
  arma::vec out = arma::zeros(numOfTerms + ncov);
  
  // time
  out[0] = Sigma(X_y_index_p - 1, X_y_index_p - 1);
  for(int k = 0; k < ncov; k++){
    out[0] += X_i[p_intercepts + k] * Sigma(p_intercepts + k, X_y_index_p - 1);
  }
  
  // covariates
  for(int k = 0; k < ncov; k++){
    out[k + numOfTerms] = Sigma(X_y_index_p - 1, p_intercepts + k);
    for(int k2 = 0; k2 < ncov; k2++){
      out[k + numOfTerms] += X_i[p_intercepts + k2] * Sigma(p_intercepts + k2, p_intercepts + k);
    }
  }
  
  double result2 = 0;
  result2 += out[0];
  
  for(int k = 0; k < ncov; k++){
    result2 += out[k + numOfTerms] * X_i[p_intercepts + k];
  }
  
  return(result2);
}

// [[Rcpp::export]]
arma::vec XSigmapX(arma::mat &X, IntegerVector X_y_index_p, 
                   arma::mat &Sigma, 
                   int p_intercepts, int ncov){
  
  
  arma::vec out = arma::zeros(X_y_index_p.size());
  
  for(int i = 0; i < X_y_index_p.size(); i++){
    
    arma::rowvec X_l = X.row(i);
    out[i] =  XiSigmapXi(X_l, X_y_index_p[i],
                         Sigma, p_intercepts, ncov);
    
  }
  
  return(out);
}

///////// SAMPLE LATENT OCCUPANCIES

arma::vec update_gamma_z_cpp_old(arma::vec gamma_z, arma::mat k_s, arma::vec mu_eps, 
                             arma::vec mu_beta_psi, arma::mat X_psi,
                             arma::vec mu_beta_p, arma::mat X_p,
                             arma::vec tX_psi_mu_beta_psi,
                             arma::vec tX_p_mu_beta_p){
  
  
  int l = 0;
  
  for(int j = 0; (unsigned)j < gamma_z.size(); j++){
    
    if(k_s(j, 2) == 0){
      
      double meanlogphi = (tX_psi_mu_beta_psi[j] + mu_eps[j]) -
        log(1 + exp(tX_psi_mu_beta_psi[j] + mu_eps[j]));
      
      double meanlog1mphi = - log(1 + exp(tX_psi_mu_beta_psi[j] + mu_eps[j]));
      
      double sum_c_i = 0;
      for(int k = 0; k < k_s(j,3); k++){
        sum_c_i += (- log(1 + exp(tX_p_mu_beta_p[l + k])));
      }
      
      double a_star = sum_c_i + meanlogphi;
      double b_star = meanlog1mphi;
      double c_star = a_star - b_star;
      
      gamma_z[j] = exp(c_star) / (1 + exp(c_star));
      
    } else {
      
      gamma_z[j] = 1;
      
    }
    
    l += k_s(j, 3);
    
  }
  
  return(gamma_z);
}  

double hessianF(double xmu){
  
  return - exp(xmu) / ((1 + exp(xmu) * (1 + exp(xmu))));
  
}

double meanlogoneminuspsi(double xmueps, 
                          double XiiSigmaXi,
                          double mu_eps_i,
                          double sd_eps_i){ 
  
  // arma::rowvec X_psi_l,
  // int X_y_index_l, //X_y_index[l]
  // int X_s_index_l, // X_s_index[l]
  // arma::mat &Sigma_beta_psi,
  // int Y,
  // int X_centers,
  // int ncov_psi_all, // numTimeSpaceCov + ncov_psi
  
  
  // arma::vec X_l = arma::conv_to<arma::vec>::from(X_psi_l);
  
  // double XiiSigmaXi = XiSigmapsiXi(X_l, X_y_index_l, X_s_index_l,
  // Sigma_beta_psi,
  // Y, X_centers, ncov_psi_all);
  
  double out = - log(1 + exp(xmueps)) + 
    .5 * hessianF(xmueps) * (XiiSigmaXi + mu_eps_i * sd_eps_i * sd_eps_i * mu_eps_i);
  
  return(out);
  
}

double meanlogoneminusp(double xmu,
                        double X_pSigmaX_p){
  // arma::rowvec X_l,
  // arma::mat &Sigma_beta_p_current,
  // int X_y_index_p_l,
  // int p_intercepts,
  // int ncov){
  
  // double X_pSigmaX_p = arma::as_scalar(X_l * Sigma_beta_p_current * arma::trans(X_l));
  
  // double X_pSigmaX_p = XiSigmapXi(X_l, X_y_index_p_l, 
  // Sigma_beta_p_current, 
  // p_intercepts, ncov);
  
  return(- log(1 + exp(xmu)) +
         .5 * hessianF(xmu) * X_pSigmaX_p);
  
}

// [[Rcpp::export]]
arma::vec update_gamma_z_cpp(arma::vec gamma_z, arma::mat k_s, 
                                    arma::vec mu_eps, arma::vec sd_eps, 
                                    arma::vec XSigmapsiX_all, arma::vec XSigmapX_all,
                                    IntegerVector X_y_index, IntegerVector X_s_index,
                                    int Y, int X_centers, int ncov_psi_all,
                                    IntegerVector &X_y_index_p,
                                    int p_intercepts, int ncov_p,
                                    arma::mat Sigma_beta_psi, arma::mat &Sigma_beta_p,
                                    arma::vec tX_psi_mu_beta_psi,
                                    arma::vec tX_p_mu_beta_p){
  
  
  int l = 0;
  
  for(int j = 0; (unsigned)j < gamma_z.size(); j++){
    
    if(k_s(j, 2) == 0){
      
      double meanlogoneminuspsi_i = meanlogoneminuspsi(tX_psi_mu_beta_psi[j] + mu_eps[j],
                                                       XSigmapsiX_all[j],
                                                                     // X_psi.row(j),
                                                                     // X_y_index[j],
                                                                     // X_s_index[j],
                                                                     // Sigma_beta_psi,
                                                                     // Y, X_centers, ncov_psi_all,
                                                                     mu_eps[j], sd_eps[j]);
      
      double meanlogphi = (tX_psi_mu_beta_psi[j] + mu_eps[j]) +
        // ( - log(1 + exp(tX_psi_mu_beta_psi[j] + mu_eps[j])));
        meanlogoneminuspsi_i;
      
      // double meanlog1mphi = - log(1 + exp(tX_psi_mu_beta_psi[j] + mu_eps[j]));
      double meanlog1mphi = meanlogoneminuspsi_i;
      
      double sum_c_i = 0;
      for(int k = 0; k < k_s(j,3); k++){
        // sum_c_i += (- log(1 + exp(tX_p_mu_beta_p[l + k])));
        double meanlogoneminusp_i = meanlogoneminusp(tX_p_mu_beta_p[l + k],
                                                     XSigmapX_all[l + k]);
        // X_p.row(l + k),
        // Sigma_beta_p,
        // X_y_index_p[l + k],
        // p_intercepts,
        // ncov_p);
        
        sum_c_i += meanlogoneminusp_i;
      }
      
      double a_star = sum_c_i + meanlogphi;
      double b_star = meanlog1mphi;
      double c_star = a_star - b_star;
      
      gamma_z[j] = exp(c_star) / (1 + exp(c_star));
      
    } else {
      
      gamma_z[j] = 1;
      
    }
    
    l += k_s(j, 3);
    
  }
  
  return(gamma_z);
}  

////////////// VB P

arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){
  
  arma::mat result(X.n_rows, D.size());
  for(int i = 0; i < result.n_rows; i++){
    for(int j = 0; j < result.n_cols; j++){
      result(i, j) = X(i,j) * D(j);
    }  
  }
  
  return(result);
}

arma::mat XtOmegaX_betap(arma::mat &X, int p_intercepts, int ncov_p, arma::vec Omega,
                         IntegerVector X_y_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
  
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  // year covariates times standard covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    for(int j = 0; j < ncov_p; j++){
      
      XtOmegaX2(p_intercepts + j, X_y_index[i] - 1) +=  X(i, p_intercepts + j) * Omega[i];
      
    }
  }
  
  for (int i = 0; i < p_intercepts; i++) {
    for (int j = 0; j < ncov_p; j++){
      XtOmegaX2(i, p_intercepts + j) = XtOmegaX2(p_intercepts + j, i);
    }
  }
  
  // standard covariates times standard covariates 
  
  for(int i = 0; i < ncov_p; i++){
    for (int j = 0; j <= i; j++) {
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX2(p_intercepts + i, p_intercepts + j) += Omega[l] * 
          X(l, p_intercepts + i) * X(l, p_intercepts + j);
      }
    }
  }
  
  for (int i = 0; i < (ncov_p- 1); i++) {
    for (int j = i; j < ncov_p; j++) {
      XtOmegaX2(p_intercepts + i, p_intercepts + j) = 
        XtOmegaX2(p_intercepts + j, p_intercepts + i);
    }
  }
  
  
  
  return(XtOmegaX2);
}

arma::vec XtransposeK_betap(arma::mat &X, IntegerVector X_y_index, 
                            arma::vec &k, int p_intercepts, int ncov){
  
  
  arma::vec Xk = arma::zeros(X.n_cols);
  
  for(int i = 0; i < X_y_index.size(); i++){
    Xk(X_y_index[i] - 1) += k[i];
  }
  for(int i = 0; i < ncov; i++){
    Xk(p_intercepts + i) = as_scalar(k.t() * X.col(p_intercepts + i));
  }
  
  return(Xk);
}

// [[Rcpp::export]]
List update_p_cpp(arma::vec mu_beta_p_current, arma::mat Sigma_beta_p_current,
                  arma::vec y, arma::vec gamma_z_all,
                  arma::vec tX_p_mu_beta_p, arma::vec XSigmapX_all,
                  arma::mat &X_p, arma::vec b_p, arma::mat inv_B_p,
                  int ncov_p,
                  int p_intercepts, IntegerVector X_y_index_p,
                  bool usingYearDetProb){
  
  
  // update pg variational params
  // arma::vec csi_p = arma::zeros(y.size());
  arma::vec diagZeta = arma::zeros(y.size());
  
  for(int l = 0; l < diagZeta.size(); l++){
    // double term1 = arma::as_scalar(X_p.row(l) * Sigma_beta_p_current * arma::trans(X_p.row(l)));
    
    // arma::rowvec X_l = X_p.row(l);
    double term1 = XSigmapX_all[l];
    // double term1 = XiSigmapXi(X_l, X_y_index_p[l], 
    // Sigma_beta_p_current, 
    // p_intercepts, ncov_p);
    double term2 = tX_p_mu_beta_p[l];// arma::as_scalar(X_p.row(l) * mu_beta_p_current);
    double phi_p = - .5 * (term1 + term2 * term2);
    double csi_p = sqrt(-2 * phi_p);
    diagZeta[l] = .5 * gamma_z_all[l] * tanh(.5 * csi_p) / csi_p;
  }
  
  // arma::vec diagZeta =  .5 * gamma_z_all % tanh(.5 * csi_p) / csi_p;
  // arma::mat tX = arma::trans(X_p);
  // arma::mat X_pZeta = diagMatrixProd(tX, diagZeta);
  // arma::mat X_pZetaX_p = X_pZeta * X_p;
  
  arma::mat X_pZetaX_p = XtOmegaX_betap(X_p, p_intercepts, 
                                        ncov_p, diagZeta,
                                        X_y_index_p);
  
  arma::mat Lambda_2 = - .5 * (inv_B_p + X_pZetaX_p);
  
  arma::mat Sigma_beta_p = arma::inv(- 2 * Lambda_2);
  
  // arma::mat Lambda_1 = arma::trans(X_p) * (gamma_z_all % (y - .5)) + inv_B_p * b_p;
  arma::vec y_tilde = gamma_z_all % (y - .5);
  arma::mat Lambda_1 = XtransposeK_betap(X_p, X_y_index_p,
                                         y_tilde, 
                                         p_intercepts, ncov_p);//arma::trans(X_p) * ) + inv_B_p * b_p;
  arma::vec mu_beta_p = Sigma_beta_p * Lambda_1;
  
  return(List::create(_["mu_beta_p"] = mu_beta_p,
                      _["Sigma_beta_p"] = Sigma_beta_p));
  // list("mu_beta_p" = mu_beta_p, "Sigma_beta_p" = Sigma_beta_p)
  // return X_pZetaX_p;
}

double mean_pg(double c){
  return tanh(c / 2) / (2 * c);
}

arma::mat tXtOmegaX_betapsi(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
                            IntegerVector X_y_index, IntegerVector X_s_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X.n_cols + X_centers, X.n_cols + X_centers);
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  if(X_centers > 0){
    // spatial  covariates times spatial covariates
    for (int i = 0; (unsigned)i < X_s_index.size(); i++){
      
      XtOmegaX2(Y + X_s_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
      
    }  
  }
  
  // year covariates times standard covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    for(int j = 0; j < ncov_psi; j++){
      
      XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + j) * Omega[i];
      
    }
  }
  
  for (int i = 0; i < Y; i++) {
    for (int j = 0; j < ncov_psi; j++){
      XtOmegaX2(i, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i);
    }
  }
  // for (int i = 1; i <= Y; i++) {
  //   for (int j = 0; j < ncov_psi; j++){
  //     XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
  //   }
  // }
  
  // spatial  covariates times year covariates
  if(Y > 0){
    if(X_centers > 0){
      for (int i = 0; (unsigned)i < X_s_index.size(); i++){
        
        XtOmegaX2(X_y_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
        
      }  
    }
  }
  
  for (int i = 0; i < Y; i++) {
    for (int j = 0; j < X_centers; j++){
      XtOmegaX2(Y + j, i) = XtOmegaX2(i, Y + j);
    }
  }
  // for (int i = 1; i <= Y; i++) {
  //   for (int j = 0; j < X_centers; j++){
  //     XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
  //   }
  // }
  
  // spatial covariates times standard covariates
  if(X_centers > 0){
    for (int i = 0; (unsigned)i < X_s_index.size(); i++){
      
      for(int j = 0; j < ncov_psi; j++){
        
        XtOmegaX2(Y + X_centers + j, Y + X_s_index[i] - 1) +=  X(i, Y + j) * Omega[i];
        
      }
    }  
    
    for (int i = 1; i <= X_centers; i++) {
      
      for (int j = 0; j < ncov_psi; j++){
        XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
      }
      
    }
  }
  
  
  
  // standard covariates times standard covariates 
  
  for(int i = 0; i < ncov_psi; i++){
    for (int j = 0; j <= i; j++) {
      // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
      // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
      // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + i) * X(l, Y + j);
      }
    }
  }
  
  for (int i = 0; i < (ncov_psi- 1); i++) {
    for (int j = i; j < ncov_psi; j++) {
      XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
        XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
    }
  }
  
  return(XtOmegaX2);
}


arma::vec XtransposeK_betapsi(arma::mat &X, IntegerVector X_y_index, 
                              IntegerVector X_s_index, arma::vec &k, 
                              int Y, int centers, int ncov){
  
  
  arma::vec Xk = arma::zeros(X.n_cols + centers);
  
  for(int i = 0; i < X_y_index.size(); i++){
    Xk(X_y_index[i] - 1) += k[i];
  }
  if(centers > 0){
    for(int i = 0; i < X_s_index.size(); i++){
      Xk(Y + X_s_index[i] - 1) += k[i];
    }  
  }
  
  for(int i = 0; i < ncov; i++){
    Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + i));
  }
  
  return(Xk);
}

// [[Rcpp::export]]
arma::vec tXmubetapsi(arma::mat &X, IntegerVector X_y_index, 
                      IntegerVector X_s_index, arma::vec mu, 
                      int Y, int centers, int ncov){
  
  int Xrows = X.n_rows;
  arma::vec Xmu = arma::zeros(Xrows);
  
  for(int i = 0; i < Xrows; i++){
    
    Xmu[i] += mu[X_y_index[i] - 1];
    
    for(int j = 0; j < ncov; j++){
      Xmu[i] += X(i, Y + j) * mu[Y + j];
    }
    
  }
  
  if(centers > 0){
    for(int i = 0; i < Xrows; i++){
      
      Xmu[i] += mu[Y + X_s_index[i] - 1];
      
    }
    
  }
  
  return(Xmu);
}

// [[Rcpp::export]]
arma::vec XiSigmapsi(arma::vec X_i, int X_y_index_i, 
                     int X_s_index_i, arma::mat Sigma, 
                     int Y, int centers, int ncov){
  
  
  int p = X_i.size();
  arma::vec out = arma::zeros(p);
  
  for(int i = 0; i < p; i++){
    double result = 0;
    result += Sigma(X_y_index_i - 1, i);
    if(centers > 0){
      result += Sigma(Y + X_s_index_i - 1, i);
    }
    for(int k = 0; k < ncov; k++){
      result += X_i[Y + k] * Sigma(Y + centers + k, i);
    }
    
    // // product for second matrix
    // for(int j = 0; j < n; j++){
    //   
    //   double result2 = 0;
    //   s
    //   
    //   Sigma[i, ]
    //   
    // }
    out[i] = result;
  }
  
  
  
  // for(int i = 0; i < X_y_index.size(); i++){
  //   Xk(X_y_index[i] - 1) += k[i];
  // }
  // for(int i = 0; i < X_s_index.size(); i++){
  //   Xk(Y + X_s_index[i] - 1) += k[i];
  // }
  // for(int i = 0; i < ncov; i++){
  //   Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
  // }
  
  return(out);
}


// [[Rcpp::export]]
List update_psi_cpp_old(arma::vec mu_beta_psi_current, 
                        arma::mat &Sigma_beta_psi_current,
                        arma::vec gamma_z, 
                        arma::mat X_psi, 
                        arma::vec b_psi, 
                        arma::mat inv_B_psi,
                        int Y,
                        int X_centers,
                        int ncov_psi,
                        int numTimeSpaceCov,
                        IntegerVector X_y_index,
                        IntegerVector X_s_index,
                        arma::vec mu_eps_current, 
                        arma::vec sd_eps_current,
                        double sigma_eps){
  
  
  // update pg variational params
  // arma::vec csi_p = arma::zeros(y.size());
  arma::vec diagZeta = arma::zeros(gamma_z.size());
  arma::vec meanpg_csi_psi = arma::zeros(gamma_z.size());
  
  for(int l = 0; l < diagZeta.size(); l++){
    // double term1 = arma::as_scalar(X_psi.row(l) * Sigma_beta_psi_current * arma::trans(X_psi.row(l)));
    
    arma::vec X_l = arma::conv_to<arma::vec>::from(X_psi.row(l));
    
    // arma::vec XiiSigma = XiSigmapsi(X_l, X_y_index[l], X_s_index[l], 
    //                                 Sigma_beta_psi_current, 
    //                                 Y, X_centers, numTimeSpaceCov + ncov_psi);
    // double XiiSigmaXi = arma::as_scalar(arma::trans(XiiSigma) * (X_l));
    double XiiSigmaXi = XiSigmapsiXi(X_l, X_y_index[l], X_s_index[l],
                                     Sigma_beta_psi_current,
                                     Y, X_centers, numTimeSpaceCov + ncov_psi);
    
    
    double Xit_mu = arma::as_scalar(X_psi.row(l) * mu_beta_psi_current);
    double term3 = 2 * Xit_mu * mu_eps_current[l];
    double term4 = sd_eps_current[l] * sd_eps_current[l] + mu_eps_current[l] * mu_eps_current[l];
    double phi_psi = - .5 * (XiiSigmaXi + Xit_mu * Xit_mu + term3 + term4);
    double csi_psi = sqrt(-2 * phi_psi);
    meanpg_csi_psi[l] = mean_pg(csi_psi);
    diagZeta[l] = .5 * tanh(.5 * csi_psi) / csi_psi;
  }
  
  // update beta psi
  
  // arma::mat tX = arma::trans(X_psi);
  // arma::mat X_psiZeta = diagMatrixProd(tX, diagZeta);
  // 
  // arma::mat X_psiZetaX_psi = X_psiZeta * X_psi;
  
  arma::mat X_psiZetaX_psi = tXtOmegaX_betapsi(X_psi, Y,
                                               X_centers,
                                               numTimeSpaceCov + ncov_psi,
                                               diagZeta,
                                               X_y_index, X_s_index);
  
  arma::mat Lambda_2 = - .5 * (inv_B_psi + X_psiZetaX_psi);
  
  // X_psiZetaX_psi <- diagMatrixProd(t(X_psi), diagZeta) %*% X_psi
  // 
  // sims_invB <- 5000
  // inv_B_psi_list <- lapply(1:sims_invB, function(i){
  //   l_T <- rgamma(1, a_l_T, b_l_T)
  //   sigma_T <- rgamma(1, a_sigma_T, b_sigma_T)
  //   l_S <- rgamma(1, a_l_s, b_l_s)
  //   sigma_S <- rgamma(1, a_sigma_s, b_sigma_s)
  //   
  //   inv_B_psi_new <- inv_B_psi
  //   
  //   K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-8), nrow = Y)
  //   
  //   inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
  //   inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]
  //   
  //   inv_B_psi_new[1:Y, 1:Y] <- inverse_Kl
  //   
  //   if(usingSpatial){
  //     K_S <- K(X_centers, X_centers, sigma_S^2, l_S) + diag(exp(-8), nrow = Y)
  //     
  //     inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_S)
  //     inverse_Kl[lower.tri(K_S)] <- t(inverse_Kl)[lower.tri(K_S)]
  //     
  //     inv_B_psi_new[Y + 1:X_centers, Y + 1:X_centers] <- inverse_Kl
  //   }
  //   
  //   inv_B_psi_new
  // })  
  // 
  // invB_psi_mean <- Reduce("+", inv_B_psi_list) / sims_invB
  
  // Lambda_2 = - .5 * (inv_B_psi + X_psiZetaX_psi);
  // Lambda_2 <- - .5 * (invB_psi_mean + X_psiZetaX_psi)
  arma::mat Sigma_beta_psi = arma::inv(- 2 * Lambda_2);
  
  arma::vec y_tilde = (gamma_z - .5) - meanpg_csi_psi % mu_eps_current;
  arma::vec Xpsiytilde = XtransposeK_betapsi(X_psi,
                                             X_y_index,
                                             X_s_index,
                                             y_tilde,
                                             Y, X_centers, numTimeSpaceCov + ncov_psi);
  // arma::vec Lambda_1 = arma::trans(X_psi) * y_tilde + inv_B_psi * b_psi;
  arma::vec Lambda_1 = Xpsiytilde + inv_B_psi * b_psi;
  // Lambda_1 <- t(X_psi) %*% y_tilde + invB_psi_mean %*% b_psi
  arma::vec mu_beta_psi = Sigma_beta_psi * Lambda_1;
  
  // update mu eps
  
  arma::vec X_psi_mu_beta_psi = X_psi * mu_beta_psi;
  arma::vec mu_eps = arma::zeros(X_psi.n_rows);
  arma::vec sd_eps = arma::zeros(X_psi.n_rows);
  for(int j = 0; j < mu_eps.size(); j++) {
    
    double mean_xib = X_psi_mu_beta_psi[j];
    double b = (gamma_z[j] - .5) - meanpg_csi_psi[j] * mean_xib;
    double a = meanpg_csi_psi[j] + 1 / (sigma_eps * sigma_eps);
    
    mu_eps[j] = b / a;
    sd_eps[j] = sqrt(1 / a);
    
  }
  
  return(List::create(_["mu_beta_psi"] = mu_beta_psi,
                      _["Sigma_beta_psi"] = Sigma_beta_psi,
                      _["mu_eps"] = mu_eps,
                      _["sd_eps"] = sd_eps));
  // list("mu_beta_p" = mu_beta_p, "Sigma_beta_p" = Sigma_beta_p)
  // return X_pZetaX_p;
}

// [[Rcpp::export]]
List update_psi_cpp(arma::vec mu_beta_psi_current, 
                    arma::mat &Sigma_beta_psi_current,
                    arma::vec XiSigmapsiXi_all,
                    arma::vec X_psi_mu_beta_psi_current,
                    arma::vec gamma_z, 
                    arma::mat &X_psi, 
                    arma::vec b_psi, 
                    arma::mat inv_B_psi,
                    arma::vec k_s, 
                    arma::vec sites,
                    int Y,
                    int X_centers,
                    int ncov_psi,
                    int numTimeSpaceCov,
                    IntegerVector X_y_index,
                    IntegerVector X_s_index,
                    arma::vec mu_eps_current, 
                    arma::vec sd_eps_current,
                    double sigma_eps){
  
  
  // update pg variational params
  arma::vec diagZeta = arma::zeros(gamma_z.size());
  arma::vec meanpg_csi_psi = arma::zeros(gamma_z.size());
  
  for(int l = 0; l < diagZeta.size(); l++){
    // double term1 = arma::as_scalar(X_psi.row(l) * Sigma_beta_psi_current * arma::trans(X_psi.row(l)));
    
    // arma::vec X_l = arma::conv_to<arma::vec>::from(X_psi.row(l));
    
    // arma::vec XiiSigma = XiSigmapsi(X_l, X_y_index[l], X_s_index[l],
    //                                 Sigma_beta_psi_current,
    //                                 Y, X_centers, numTimeSpaceCov + ncov_psi);
    // double XiiSigmaXi = arma::as_scalar(arma::trans(XiiSigma) * (X_l));
    // double XiiSigmaXi = XiSigmapsiXi(X_l, X_y_index[l], X_s_index[l],
    // Sigma_beta_psi_current,
    // Y, X_centers, numTimeSpaceCov + ncov_psi);
    double XiiSigmaXi = XiSigmapsiXi_all[l];
    
    // double Xit_mu = arma::as_scalar(X_psi.row(l) * mu_beta_psi_current); 
    double Xit_mu = X_psi_mu_beta_psi_current[l]; 
    double term3 = 2 * Xit_mu * mu_eps_current[l];
    double term4 = sd_eps_current[l] * sd_eps_current[l] + mu_eps_current[l] * mu_eps_current[l];
    double phi_psi = - .5 * (XiiSigmaXi + Xit_mu * Xit_mu + term3 + term4);
    double csi_psi = sqrt(-2 * phi_psi);
    meanpg_csi_psi[l] = mean_pg(csi_psi);
    diagZeta[l] = .5 * tanh(.5 * csi_psi) / csi_psi;
  }
  
  // update beta psi
  
  // arma::mat tX = arma::trans(X_psi);
  // arma::mat X_psiZeta = diagMatrixProd(tX, diagZeta);
  // arma::mat X_psiZetaX_psi = X_psiZeta * X_psi;
  
  arma::mat X_psiZetaX_psi = tXtOmegaX_betapsi(X_psi, Y,
                                               X_centers,
                                               numTimeSpaceCov + ncov_psi,
                                               diagZeta,
                                               X_y_index, X_s_index);
  
  arma::mat Lambda_2 = - .5 * (inv_B_psi + X_psiZetaX_psi);
  
  arma::mat Sigma_beta_psi = arma::inv(- 2 * Lambda_2);
  
  arma::vec y_tilde = (gamma_z - .5) - meanpg_csi_psi % mu_eps_current;
  arma::vec Xpsiytilde = XtransposeK_betapsi(X_psi,
                                             X_y_index,
                                             X_s_index,
                                             y_tilde,
                                             Y, X_centers, numTimeSpaceCov + ncov_psi);
  arma::vec Lambda_1 = Xpsiytilde + inv_B_psi * b_psi;
  
  arma::vec mu_beta_psi = Sigma_beta_psi * Lambda_1;
  // arma::vec mu_beta_psi = mu_beta_psi_current;
  // arma::mat Sigma_beta_psi = Sigma_beta_psi_current;
  
  // eps
  
  // arma::vec X_psi_mu_beta_psi = X_psi * mu_beta_psi;
  arma::vec X_psi_mu_beta_psi = tXmubetapsi(X_psi, X_y_index, X_s_index,
                                            mu_beta_psi, Y, X_centers,
                                            ncov_psi + numTimeSpaceCov);
  
  arma::vec mu_eps = arma::zeros(k_s.size());
  arma::vec sd_eps = arma::zeros(k_s.size());
  
  int index_site = 0;
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    int l = 1;
    
    // find rows associated with current site
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    double sum_a = 0;
    double sum_b = 0;
    for(int j = 0; j < l; j++){
      sum_a += meanpg_csi_psi[indexes_site[j]];
      sum_b += (gamma_z[indexes_site[j]] - .5) -
        meanpg_csi_psi[indexes_site[j]] * X_psi_mu_beta_psi[indexes_site[j]];
    }
    
    double b = sum_b;
    double a = sum_a + 1 / (sigma_eps * sigma_eps);
    
    for(int j = 0; j < l; j++){
      mu_eps[indexes_site[j]] = b / a;
      sd_eps[indexes_site[j]] = sqrt(1 / a);
    }
    
  }
  
  
  // return(List::create(_["mu_beta_psi"] = mu_beta_psi_current,
  // _["Sigma_beta_psi"] = Sigma_beta_psi_current,
  return(List::create(_["mu_beta_psi"] = mu_beta_psi,
                      _["Sigma_beta_psi"] = Sigma_beta_psi,
                      _["mu_eps"] = mu_eps,
                      _["sd_eps"] = sd_eps));
}

///////// GOODNESS OF FIT

// [[Rcpp::export]]
arma::vec simulateDetections(arma::vec p, arma::vec z_all){
  
  arma::vec y = arma::zeros(p.size());
  
  // this loops over p
  for(int i = 0; (unsigned)i < y.size(); i++){
    
    if(z_all(i) == 1){
      
      y[i] = R::rbinom(1, p[i]);
      
    } 
    
  }
  
  return(y);
}  

arma::vec logit(arma::vec x){
  return(1 / (1 + exp(-x)));  
}

// [[Rcpp::export]]
arma::vec computeYearEffect(int Y, arma::vec a_s_unique, arma::vec beta_psi){
  
  arma::vec yearEffect = arma::zeros(Y);
  
  for(int y = 0; (unsigned)y < Y; y++){
    
    arma::vec occProbs = logit(beta_psi[y] + a_s_unique);
    
    yearEffect[y] = mean(occProbs);
  }
  
  return(yearEffect);
}  

// GAUSSIAN PROCESS FUNCTIONS

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());
  
  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }  
  }
  
  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);
  
  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }  
  }
  
  return res;
}

// CREATE GRID

// [[Rcpp::export]]
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) < x_grid[i + 1]){
        return(true);
      }
    } 
    
  }
  
  return(false);
}

// [[Rcpp::export]]
bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j) {
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) > x_grid[i - 1]){
        return(true);
      }
    } 
    
  }
  
  return(false);
}

// [[Rcpp::export]]
bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) > y_grid[j-1]){
        return(true);
      }
    }
    
  }
  
  return(false);
  
}

// [[Rcpp::export]]
bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) < y_grid[j+1]){
        return(true);
      }
    }
    
  }
  
  return(false);
  
}

// [[Rcpp::export]]
IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde){
  
  IntegerVector closestPoint(XY_sp.n_rows);
  
  for(int k = 0; k < XY_sp.n_rows; k++){
    
    double newDistance = 0;
    double minDistance = exp(50);
    int bestIndex = 0;
    
    for(int i = 0; i < X_tilde.n_rows; i++){
      newDistance = pow(X_tilde(i, 0) - XY_sp(k, 0), 2) + pow(X_tilde(i, 1) - XY_sp(k, 1), 2);
      
      if(newDistance < minDistance){
        minDistance = newDistance;
        bestIndex = i + 1;
      }
    }
    
    closestPoint[k] = bestIndex;
    
  }
  
  return(closestPoint);
}

/// HYPERPARAMETER

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// // [[Rcpp::export]]
// double dmvnrm_arma_fast(arma::vec const &x,  
//                            arma::rowvec const &mean,  
//                            arma::mat const &sigma, 
//                            bool const logd = false) { 
//   using arma::uword;
//   uword const xdim = x.size();
//   double out;
//   arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
//   double const rootisum = arma::sum(log(rooti.diag())), 
//     constants = -(double)xdim/2.0 * log2pi, 
//     other_terms = rootisum + constants;
//   
//   arma::rowvec x_row = arma::conv_to<arma::rowvec>::from(x);
//   arma::rowvec z = (x_row - mean);
//   inplace_tri_mat_mult(z, rooti);
//   out = other_terms - 0.5 * arma::dot(z, z);     
//   
//   if (logd)
//     return out;
//   return exp(out);
// }
// [[Rcpp::export]]
double dmvnrm_arma_fast(arma::vec const &x,   
                        arma::mat const &sigma, 
                        bool const logd = false) { 
  using arma::uword;
  uword const xdim = x.size();
  double out;
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z = arma::conv_to<arma::rowvec>::from(x);
  inplace_tri_mat_mult(z, rooti);
  out = other_terms - 0.5 * arma::dot(z, z);     
  
  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
double logp_T(arma::vec y, double l, double sig,
              arma::vec x_gp, double sig2){
  
  int ydim = y.size();
  
  arma::mat K_l = K(x_gp, x_gp, sig, l);
  
  if(sig2 > 0){
    arma::mat sig2_mat = sig2 * sig2 * arma::ones(ydim, ydim);
    K_l += sig2_mat;
  } 
  
  arma::mat diagMat = arma::eye(ydim, ydim);
  K_l += .001 * diagMat;
  
  return(dmvnrm_arma_fast(y, K_l, true));
}

// [[Rcpp::export]]
double logp_s(arma::vec y, double l, double sig,
              arma::mat x_gp, double sig2){
  
  int ydim = y.size();
  
  arma::mat K_l = K2(x_gp, x_gp, sig, l);
  
  if(sig2 > 0){
    arma::mat sig2_mat = sig2 * sig2 * arma::ones(ydim, ydim);
    K_l += sig2_mat;
  } 
  
  arma::mat diagMat = arma::eye(ydim, ydim);
  K_l += .001 * diagMat;
  
  return(dmvnrm_arma_fast(y, K_l, true));
}

double logp0(double x, double a_x, double b_x){
  return R::dgamma(x, a_x, 1 / b_x, 1);
}

double logq(double x, double a, double b){
  return R::dgamma(x, a, 1 / b, 1);
}

arma::vec Grad_logq_cpp(double x, double a, double b){
  
  arma::vec out = arma::zeros(2);
  out[0] = log(b) - R::digamma(a) + log(x);
  out[1] = a / b - x;
  
  return(out);
}

double cov_rcpp(arma::vec x, 
                arma::vec y){
  
  double mean_x = mean(x);  
  double mean_y = mean(y);  
  
  double out = 0;
  
  int n = x.size();
  for(int i = 0; i < n; i++){
    out += (x[i] - mean_x) * (y[i] - mean_y);
  }
  
  return out / (n - 1);
}

arma::mat G_lambda(double a, double b){
  
  arma::mat out = arma::zeros(2,2);
  out(0,0) = R::trigamma(a);
  out(1, 0) = - 1 / b;
  out(0, 1) = - 1 / b;
  out(1, 1) = a / (b * b);
  
  return out;
}

// [[Rcpp::export]]
List update_sigmasq_T_cpp(arma::vec y,
                          double l_T,
                          double sig2,
                          arma::vec x_gp,
                          double a_sig_0,
                          double b_sig_0){
  
  double a_sig = a_sig_0;
  double b_sig = b_sig_0;
  
  int Y = y.size();
  
  a_sig += Y / 2.0;
  
  arma::mat K_l = K(x_gp, x_gp, 1, l_T);// + diag(exp(-10), nrow = nrow(X_tilde))
  
  if(sig2 > 0){
    arma::mat sig2_mat = sig2 * sig2 * arma::ones(Y, Y);
    K_l += sig2_mat;
  } 
  
  arma::mat diagMat = arma::eye(Y, Y);
  K_l += .001 * diagMat;
  
  // arma::mat inv_chol_Kl = arma::inv(arma::chol(K_l));
  // arma::vec Ltbt = inv_chol_Kl * y;
  double b_ig = arma::as_scalar(arma::trans(y) * arma::inv(K_l) * y) / 2;
  
  b_sig += b_ig;
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

// [[Rcpp::export]]
List update_sigmasq_S_cpp(arma::vec y,
                          double l_S,
                          arma::mat x_gp,
                          double sig2,
                          double a_sig_0,
                          double b_sig_0){
  
  double a_sig = a_sig_0;
  double b_sig = b_sig_0;
  
  int centers = y.size();
  
  a_sig += centers / 2.0;
  
  arma::mat K_l = K2(x_gp, x_gp, 1, l_S);// + diag(exp(-10), nrow = nrow(X_tilde))
  
  if(sig2 > 0){
    arma::mat sig2_mat = sig2 * sig2 * arma::ones(centers, centers);
    K_l += sig2_mat;
  } 
  
  arma::mat diagMat = arma::eye(centers, centers);
  K_l += .001 * diagMat;
  
  // arma::mat inv_chol_Kl = arma::inv(arma::chol(K_l));
  // arma::vec Ltbt = inv_chol_Kl * y;
  double b_ig = arma::as_scalar(arma::trans(y) * arma::inv(K_l) * y) / 2;
  
  b_sig += b_ig;
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

// [[Rcpp::export]]
List update_sigmasq_S_cpp_fast(arma::vec y,
                               double l_S,
                               arma::mat inv_Kl,
                               double a_sig_0,
                               double b_sig_0){
  
  double a_sig = a_sig_0;
  double b_sig = b_sig_0;
  
  int centers = y.size();
  
  a_sig += centers / 2.0;
  
  double b_ig = arma::as_scalar(arma::trans(y) * inv_Kl * y) / 2;
  
  b_sig += b_ig;
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

List update_sig_T_cpp(arma::vec y,
                      double a_l,
                      double b_l,
                      double a_sig,
                      double b_sig,
                      double sig2,
                      arma::vec x_gp,
                      double a_sig_0,
                      double b_sig_0,
                      int sims_grad,
                      double rho){
  
  arma::mat grad_sims = arma::zeros(sims_grad, 5);
  for(int s = 0; s < sims_grad; s++){
    double l_s = R::rgamma(a_l, 1 / b_l);
    double sig_s = R::rgamma(a_sig, 1/ b_sig);
    arma::vec h_i = Grad_logq_cpp(sig_s, a_sig, b_sig);
    double d_i = logp_T(y, l_s, sig_s * sig_s, x_gp, sig2) +
      logp0(sig_s, a_sig_0, b_sig_0) -
      logq(sig_s, a_sig, b_sig);
    
    grad_sims(s, 0) = h_i[0];
    grad_sims(s, 1) = h_i[1];
    grad_sims(s, 2) = d_i;
    grad_sims(s, 3) = d_i * h_i[0];
    grad_sims(s, 4) = d_i * h_i[1];
    
  }
  
  double sum_covhifi = 0;
  double sum_varhi = 0;
  for(int d = 0; d < 2; d++){
    arma::vec h_i = arma::conv_to<arma::vec>::from(grad_sims.col(d));
    arma::vec f_i = arma::conv_to<arma::vec>::from(grad_sims.col(d + 3));
    sum_covhifi += cov_rcpp(h_i, f_i);
    sum_varhi += cov_rcpp(h_i, h_i);
  }
  double a_i = sum_covhifi / sum_varhi;
  
  arma::vec grad_est = arma::zeros(2);
  for(int d = 0; d < 2; d++){
    grad_est[d] = mean(grad_sims.col(d + 3) - a_i * grad_sims.col(d));
  }
  
  arma::mat Glambdamat = G_lambda(a_sig, b_sig);
  
  arma::vec nat_grad = arma::inv(Glambdamat) * grad_est;
  
  a_sig += rho * nat_grad[0];
  b_sig += rho * nat_grad[1];
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

List update_sig_T_delta_cpp(arma::vec y,
                            double l_s,
                            double a_sig,
                            double b_sig,
                            double sig2,
                            arma::vec x_gp,
                            double a_sig_0,
                            double b_sig_0,
                            int sims_grad,
                            double rho){
  
  arma::mat grad_sims = arma::zeros(sims_grad, 5);
  for(int s = 0; s < sims_grad; s++){
    double sig_s = R::rgamma(a_sig, 1/ b_sig);
    arma::vec h_i = Grad_logq_cpp(sig_s, a_sig, b_sig);
    double d_i = logp_T(y, l_s, sig_s * sig_s, x_gp, sig2) +
      logp0(sig_s, a_sig_0, b_sig_0) -
      logq(sig_s, a_sig, b_sig);
    
    grad_sims(s, 0) = h_i[0];
    grad_sims(s, 1) = h_i[1];
    grad_sims(s, 2) = d_i;
    grad_sims(s, 3) = d_i * h_i[0];
    grad_sims(s, 4) = d_i * h_i[1];
    
  }
  
  double sum_covhifi = 0;
  double sum_varhi = 0;
  for(int d = 0; d < 2; d++){
    arma::vec h_i = arma::conv_to<arma::vec>::from(grad_sims.col(d));
    arma::vec f_i = arma::conv_to<arma::vec>::from(grad_sims.col(d + 3));
    sum_covhifi += cov_rcpp(h_i, f_i);
    sum_varhi += cov_rcpp(h_i, h_i);
  }
  double a_i = sum_covhifi / sum_varhi;
  
  arma::vec grad_est = arma::zeros(2);
  for(int d = 0; d < 2; d++){
    grad_est[d] = mean(grad_sims.col(d + 3) - a_i * grad_sims.col(d));
  }
  
  arma::mat Glambdamat = G_lambda(a_sig, b_sig);
  
  arma::vec nat_grad = arma::inv(Glambdamat) * grad_est;
  
  a_sig += rho * nat_grad[0];
  b_sig += rho * nat_grad[1];
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

List update_sig_S_cpp(arma::vec y,
                      double a_l,
                      double b_l,
                      double a_sig,
                      double b_sig,
                      double sig2,
                      arma::mat x_gp,
                      double a_sig_0,
                      double b_sig_0,
                      int sims_grad,
                      double rho){
  
  arma::mat grad_sims = arma::zeros(sims_grad, 5);
  for(int s = 0; s < sims_grad; s++){
    double l_s = R::rgamma(a_l, 1 / b_l);
    double sig_s = R::rgamma(a_sig, 1/ b_sig);
    arma::vec h_i = Grad_logq_cpp(sig_s, a_sig, b_sig);
    double d_i = logp_s(y, l_s, sig_s * sig_s, x_gp, sig2) +
      logp0(sig_s, a_sig_0, b_sig_0) -
      logq(sig_s, a_sig, b_sig);
    
    grad_sims(s, 0) = h_i[0];
    grad_sims(s, 1) = h_i[1];
    grad_sims(s, 2) = d_i;
    grad_sims(s, 3) = d_i * h_i[0];
    grad_sims(s, 4) = d_i * h_i[1];
    
  }
  
  double sum_covhifi = 0;
  double sum_varhi = 0;
  for(int d = 0; d < 2; d++){
    arma::vec h_i = arma::conv_to<arma::vec>::from(grad_sims.col(d));
    arma::vec f_i = arma::conv_to<arma::vec>::from(grad_sims.col(d + 3));
    sum_covhifi += cov_rcpp(h_i, f_i);
    sum_varhi += cov_rcpp(h_i, h_i);
  }
  double a_i = sum_covhifi / sum_varhi;
  
  arma::vec grad_est = arma::zeros(2);
  for(int d = 0; d < 2; d++){
    grad_est[d] = mean(grad_sims.col(d + 3) - a_i * grad_sims.col(d));
  }
  
  arma::mat Glambdamat = G_lambda(a_sig, b_sig);
  
  arma::vec nat_grad = arma::inv(Glambdamat) * grad_est;
  
  a_sig += rho * nat_grad[0];
  b_sig += rho * nat_grad[1];
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

List update_sig_S_delta_cpp(arma::vec y,
                            double l_s,
                            double a_sig,
                            double b_sig,
                            double sig2,
                            arma::mat x_gp,
                            double a_sig_0,
                            double b_sig_0,
                            int sims_grad,
                            double rho){
  
  arma::mat grad_sims = arma::zeros(sims_grad, 5);
  for(int s = 0; s < sims_grad; s++){
    
    double sig_s = R::rgamma(a_sig, 1/ b_sig);
    arma::vec h_i = Grad_logq_cpp(sig_s, a_sig, b_sig);
    double d_i = logp_s(y, l_s, sig_s * sig_s, x_gp, sig2) +
      logp0(sig_s, a_sig_0, b_sig_0) -
      logq(sig_s, a_sig, b_sig);
    
    grad_sims(s, 0) = h_i[0];
    grad_sims(s, 1) = h_i[1];
    grad_sims(s, 2) = d_i;
    grad_sims(s, 3) = d_i * h_i[0];
    grad_sims(s, 4) = d_i * h_i[1];
    
  }
  
  double sum_covhifi = 0;
  double sum_varhi = 0;
  for(int d = 0; d < 2; d++){
    arma::vec h_i = arma::conv_to<arma::vec>::from(grad_sims.col(d));
    arma::vec f_i = arma::conv_to<arma::vec>::from(grad_sims.col(d + 3));
    sum_covhifi += cov_rcpp(h_i, f_i);
    sum_varhi += cov_rcpp(h_i, h_i);
  }
  double a_i = sum_covhifi / sum_varhi;
  
  arma::vec grad_est = arma::zeros(2);
  for(int d = 0; d < 2; d++){
    grad_est[d] = mean(grad_sims.col(d + 3) - a_i * grad_sims.col(d));
  }
  
  arma::mat Glambdamat = G_lambda(a_sig, b_sig);
  
  arma::vec nat_grad = arma::inv(Glambdamat) * grad_est;
  
  a_sig += rho * nat_grad[0];
  b_sig += rho * nat_grad[1];
  
  return(List::create(_["a_sig"] = a_sig,
                      _["b_sig"] = b_sig));
}

List update_l_s_cpp(arma::vec y,
                    double a_l,
                    double b_l,
                    double a_sig,
                    double b_sig,
                    double sig2,
                    arma::mat x_gp,
                    double a_l_0,
                    double b_l_0,
                    int sims_grad,
                    double rho){
  
  arma::mat grad_sims = arma::zeros(sims_grad, 5);
  for(int s = 0; s < sims_grad; s++){
    double l_s = R::rgamma(a_l, 1 / b_l);
    double sig_s = R::rgamma(a_sig, 1 / b_sig);
    arma::vec h_i = Grad_logq_cpp(l_s, a_l, b_l);
    double d_i = logp_s(y, l_s, sig_s * sig_s, x_gp, sig2) +
      logp0(l_s, a_l_0, b_l_0) -
      logq(l_s, a_l, b_l);
    
    grad_sims(s, 0) = h_i[0];
    grad_sims(s, 1) = h_i[1];
    grad_sims(s, 2) = d_i;
    grad_sims(s, 3) = d_i * h_i[0];
    grad_sims(s, 4) = d_i * h_i[1];
    
  }
  
  double sum_covhifi = 0;
  double sum_varhi = 0;
  for(int d = 0; d < 2; d++){
    arma::vec h_i = arma::conv_to<arma::vec>::from(grad_sims.col(d));
    arma::vec f_i = arma::conv_to<arma::vec>::from(grad_sims.col(d + 3));
    sum_covhifi += cov_rcpp(h_i, f_i);
    sum_varhi += cov_rcpp(h_i, h_i);
  }
  double a_i = sum_covhifi / sum_varhi;
  
  arma::vec grad_est = arma::zeros(2);
  for(int d = 0; d < 2; d++){
    grad_est[d] = mean(grad_sims.col(d + 3) - a_i * grad_sims.col(d));
  }
  
  arma::mat Glambdamat = G_lambda(a_l, b_l);
  
  arma::vec nat_grad = arma::inv(Glambdamat) * grad_est;
  
  a_l += rho * nat_grad[0];
  b_l += rho * nat_grad[1];
  
  return(List::create(_["a_l"] = a_l,
                      _["b_l"] = b_l));
}


// [[Rcpp::export]]
arma::vec computeElboVals(arma::vec y,
                          arma::mat pointsToSearch,
                          double sig2,
                          arma::mat x_gp,
                          double a_l_S_0,
                          double b_l_S_0,
                          double a_sigma_S_0,
                          double b_sigma_S_0){
  
  int pointsGrid = pointsToSearch.n_rows;
  arma::vec elbo_vals = arma::zeros(pointsGrid);
  for(int i = 0; i < pointsGrid; i++){
    
    double sigma_S_current = pointsToSearch(i, 0);
    double l_S_current = pointsToSearch(i, 1);
    
    double logp_term = logp_s(y,
                              l_S_current, 
                              sigma_S_current * sigma_S_current, 
                              x_gp, 
                              sig2);
    
    elbo_vals[i] = logp_term + 
      R::dgamma(l_S_current, a_l_S_0, 1 / b_l_S_0, 1) +
      R::dgamma(sigma_S_current, a_sigma_S_0, 1 / b_sigma_S_0, 1);
    
  }
  
  return(elbo_vals);
}

//// HYPERPARAMETER SAMPLING FUNCTIONS

// [[Rcpp::export]]
arma::mat matrixProductXtOmegaX_year(int Y, arma::vec Omega,
                                     arma::vec X_y_index){
  
  arma::mat XtOmegaX2 = arma::zeros(Y, Y);
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  return(XtOmegaX2);
}

// [[Rcpp::export]]
arma::mat matrixProductXtOmegaX_spatial(int X_centers, arma::vec Omega,
                                        arma::vec X_s_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X_centers, X_centers);
  
  // spatial  covariates times spatial covariates
  for (int i = 0; (unsigned)i < X_s_index.size(); i++){
    
    XtOmegaX2(X_s_index[i] - 1, X_s_index[i] - 1) += Omega[i];
    
  }
  
  return(XtOmegaX2);
}

// [[Rcpp::export]]
arma::vec XpsiYz(arma::vec X_y_index, arma::vec z, int Y){
  
  arma::vec out = arma::zeros(Y);
  
  for(int i = 0; i < X_y_index.size(); i++){
    out(X_y_index[i] - 1) += z[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
arma::vec XpsinoYbetaz(arma::vec X_s_index, int X_centers, int ncov,
                       arma::mat &X_cov, arma::vec& beta){
  
  // X_psi[,-(1:Y),drop=F] %*% beta_psi[-(1:Y)]
  
  return(X_cov * beta);
  // arma::vec out = arma::zeros(X.n_cols);
  
  // for(int i = 0; i < X_s_index.size(); i++){
  //   Xk(Y + X_s_index[i] - 1) += k[i];
  // }
  // for(int i = 0; i < ncov; i++){
  //   Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
  // }
  
}


// [[Rcpp::export]]
double rcpp_log_dmvnorm_fast_cpp(arma::mat &inv_S, arma::vec &diag_S, 
                                 double sigma_s, arma::vec &x) {
  int n = x.size();
  
  double term1 = (-n / 2.0) * log(2 * M_PI) - sum(log(diag_S) + log(sigma_s));
  
  arma::vec term2 = inv_S * x;
  
  arma::vec term3 = arma::trans(x) * term2;
  
  double term4 = (.5) * (1 / (sigma_s*sigma_s)) * term3[0];
  
  // return(0);
  return(term1 + term4);
}

// [[Rcpp::export]]
arma::vec sample_l_grid_cpp(arma::vec l_s_grid, double sigma_s, 
                            arma::cube &inv_K_s_grid, arma::mat &diag_K_s_grid,
                            double a_l_S, double b_l_S, arma::vec a_s){
  
  arma::vec posterior_val = arma::zeros(l_s_grid.size());
  
  for(int j = 0; j < l_s_grid.size(); j++){
    // for (j in 1:length(l_s_grid)) {
    
    double l_s = l_s_grid[j];
    
    arma::mat inv_K_s_grid_j = inv_K_s_grid.subcube(arma::span(), arma::span(), arma::span(j));
    arma::vec diag_K_s_grid_j = arma::conv_to<arma::vec>::from(diag_K_s_grid.col(j));
    // # K_s_grid_j <- K_s_grid[,,j] * sigma_s^2
    // # inv_K_s_grid_j <- inv_K_s_grid[,,j] / sigma_s^2
    // # diag_K_s_grid_j <- diag_K_s_grid[,j] * sigma_s
    
    double loglikelihood = 0;//rcpp_log_dmvnorm_fast_cpp(inv_K_s_grid_j, diag_K_s_grid_j, 
    //                      sigma_s, a_s);
    
    // loglikelihood <- rcpp_log_dmvnorm_fast(inv_K_s_grid[,,j], 
    // diag_K_s_grid[,j], sigma_s, a_s)
    
    // # loglikelihood <- rcpp_log_dmvnorm_fast(1, inv_K_s_grid_j, 
    // #                                        diag_K_s_grid_j,  a_s)
    //     
    // # Sigma_l <- K2(X_tilde, X_tilde, sigma_s^2, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
    //     
    // # (loglikelihood2 <- rcpp_log_dmvnorm( Sigma_l, rep(0, X_centers), a_s, F))
    
    double logPrior = R::dgamma(l_s, a_l_S, 1 / b_l_S, 1);
    
    posterior_val[j] = logPrior + loglikelihood;
    
  }
  
  return(posterior_val);
}


///////// ELBO CALCULATIONS

// [[Rcpp::export]]
double log_q_z(arma::vec gamma_z, 
               arma::vec z){
  
  double out = 0;
  for(int i = 0; i < gamma_z.size(); i++){
    
    
    out += R::dbinom(z[i], 1, gamma_z[i], 1);
    
  }
  
  return out;  
} 

// [[Rcpp::export]]
double log_q_eps(arma::vec eps_s, 
                 arma::vec mu_eps,
                 arma::vec sd_eps){
  
  double out = 0;
  for(int i = 0; i < mu_eps.size(); i++){
    
    out += R::dnorm(eps_s[i], mu_eps[i], sd_eps[i], 1);
  }
  
  return out;  
} 

// [[Rcpp::export]]
arma::vec simulate_q_z(arma::vec gamma_z){
  
  arma::vec z = arma::zeros(gamma_z.size());
  for(int i = 0; i < gamma_z.size(); i++){
    
    z[i] = R::rbinom(1, gamma_z[i]);
    
  }
  
  return z;  
} 

// [[Rcpp::export]]
arma::vec simulate_eps(arma::vec mu_eps,
                       arma::vec sd_eps){
  
  arma::vec eps_s = arma::zeros(mu_eps.size());
  for(int i = 0; i < mu_eps.size(); i++){
    
    eps_s[i] = R::rnorm(mu_eps[i], sd_eps[i]);
    
  }
  
  return eps_s;  
} 

// [[Rcpp::export]]
double y_lik(arma::vec y,
             arma::vec p,
             arma::vec z_all){
  
  double out = 0;
  for(int i = 0; i < y.size(); i++){
    if(z_all[i] == 1){
      out += R::dbinom(y[i], 1, p[i], 1);
    }
  }
  
  return out;
}

// // [[Rcpp::export]]
// arma::vec XpsinoYbetaz(arma::mat X_tilde_star, arma::vec indexesSite,
//                        arma::mat beta_psi_output, int year){
//   
//   arma::vec sitesProb = arma::zeros(X_tilde_star.n_rows); 
//   
//   int niter = beta_psi_output.n_rows;
//   
//   for(int i = 0; i < sitesProb.n_rows; i++){
//     
//     arma::vec siteProbs = arma::zeros(niter);
//     
//     for(int j = 0; j < niter; j++){
//       
//       siteProbs[j] = beta_psi_output(j, year) + beta_psi_output(j,Y + indexesSite) + 
//         X_tilde_star[indexesSite,1] * X_psi_yearcov_values[year] * 
//         beta_psi_output[,Y + X_centers + 1] +
//         X_tilde_star[indexesSite,2] * X_psi_yearcov_values[year] *
//         beta_psi_output[,Y + X_centers + 2])
//       
//     }
//     
//     sitesProb[i] = mean(logit(siteProbs));
//     
//   }
//   
//   <- sapply(1:nrow(X_tilde_star), function(i){
//   mean(logit(beta_psi_output[,year] + beta_psi_output[,Y + indexesSite] + 
//     X_tilde_star[indexesSite,1] * X_psi_yearcov_values[year] * 
//     beta_psi_output[,Y + X_centers + 1] +
//     X_tilde_star[indexesSite,2] * X_psi_yearcov_values[year] *
//     beta_psi_output[,Y + X_centers + 2]))
// })
//   
// }

//// OLD STUFF

// // [[Rcpp::export]]
// arma::vec computelikelihood_cpp(arma::vec Occs, arma::vec p, arma::vec z_all){
//   
//   arma::vec likelihoods = arma::zeros(Occs.size());
//   
//   for(int i = 0; i < Occs.size(); i++){
//     if(z_all[i] == 0){
//       likelihoods[i] = 1;
//     } else {
//       likelihoods[i] = R::dbinom(Occs[i], 1, p[i], 0);
//     }
//   }
//   
//   return(likelihoods);
// }
// 
// List sampler_beta_sp_old(arma::vec beta,
//                          arma::vec a_s,
//                          arma::mat& X, 
//                          arma::vec b, 
//                          arma::mat invB, 
//                          arma::vec n, 
//                          arma::vec k, 
//                          int Y, 
//                          int X_centers,
//                          int ncov_psi,
//                          IntegerVector X_y_index,
//                          IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_sp(X, Y, X_centers, 5 + ncov_psi, Omega,
//                                                 X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast(X, invB, b, knew, XtOmegaX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// List sampler_beta_sp_old2(arma::vec beta,
//                           arma::vec a_s,
//                           arma::mat& X, 
//                           arma::vec b, 
//                           arma::mat invB, 
//                           arma::vec n, 
//                           arma::vec k,
//                           int Y, 
//                           int X_centers,
//                           int ncov_psi,
//                           IntegerVector X_y_index,
//                           IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_sp(X, Y, X_centers, 5 + ncov_psi, Omega,
//                                                 X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
//                                      X_y_index, X_s_index, Y, X_centers, 5 + ncov_psi);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// List sampler_beta_integrated(arma::vec beta,
//                              arma::vec a_s,
//                              arma::vec mu_i,
//                              arma::vec sigmasq_i,
//                              arma::mat X, 
//                              arma::vec b, 
//                              arma::mat B, 
//                              arma::vec n, 
//                              arma::vec k, 
//                              int Y, 
//                              int ncov_psi,
//                              arma::vec X_y_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, ncov_psi, Omega,
//                                              X_y_index);
//   
//   // create new variables for integrated posterior
//   arma::vec ktilde = arma::zeros(k.size());
//   for(int i = 0; i < ktilde.size(); i++){
//     
//     ktilde[i] = k[i] - (Omega[i] / (sigmasq_i[i] * Omega[i] + 1)) * (sigmasq_i[i] * k[i] + mu_i[i]);
//     
//   }
//   
//   arma::vec Omegatilde = arma::zeros(Omega.size());
//   for(int i = 0; i < Omegatilde.size(); i++){
//     
//     Omegatilde[i] = Omega[i] / (sigmasq_i[i] * Omega[i] + 1);
//     
//   } 
//   
//   arma::mat XtOmegaTildeX = matrixProductXtOmegaX(X, Y, ncov_psi, Omegatilde,
//                                                   X_y_index);
//   
//   beta = sample_beta_cpp_fast(X, B, b, ktilde, XtOmegaTildeX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega,
//                       _["XtOmegaX"] = XtOmegaX));
// }
// 
// Eigen::VectorXd sample_beta_cpp_fast_eigen(Eigen::MatrixXd& X, Eigen::MatrixXd& invB, Eigen::VectorXd& b, 
//                                            Eigen::VectorXd& k, Eigen::MatrixXd& XtOmegaX){
//   
//   
//   Eigen::MatrixXd Lambda_B = XtOmegaX + invB;
//   Eigen::VectorXd mu_B =  X.transpose() * k + invB * b;
//   
//   Eigen::MatrixXd L = Lambda_B.llt().matrixL();
//   Eigen::VectorXd tmp = L.triangularView<Eigen::Lower>().solve(mu_B);
//   Eigen::VectorXd mu = L.transpose().triangularView<Eigen::Upper>().solve(tmp);
//   
//   Eigen::VectorXd z = as<Eigen::VectorXd>(rnorm(L.cols()));
//   
//   Eigen::VectorXd v = L.transpose().triangularView<Eigen::Upper>().solve(z);
//   Eigen::VectorXd out = mu + v;
//   
//   // Eigen::MatrixXd id_matrix = Eigen::MatrixXd::Identity(L.cols(), L.cols());
//   // Eigen::MatrixXd cholSigma = L.triangularView<Eigen::Lower>().solve(id_matrix);
//   
//   
//   
//   // Eigen::VectorXd out = mu + cholSigma * z;
//   
//   
//   return out;
//   
//   //
//   
//   // Eigen::MatrixXd tX = X.transpose();
//   // 
//   // Eigen::MatrixXd Lambda_B = XtOmegaX + invB;
//   // Eigen::MatrixXd L = Lambda_B.llt().matrixL().transpose();
//   // Eigen::VectorXd tmp = L.triangularView<Eigen::Lower>().solve(tX * k + invB * b);
//   // Eigen::VectorXd alpha = L.transpose().triangularView<Eigen::Upper>().solve(tmp);
//   // 
//   // // Eigen::MatrixXd L = arma::trans(arma::chol());
//   // // arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + invB * b);
//   // // arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   // Eigen::MatrixXd id_matrix = Eigen::MatrixXd::Identity(invB_arma.n_cols,invB_arma.n_cols);
//   // Eigen::MatrixXd tinvL = L.triangularView<Eigen::Lower>().solve(id_matrix).transpose();
//   // 
//   // // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   // 
//   // arma::mat tinvL_arma = arma::mat(tinvL.data(), tinvL.rows(), tinvL.cols(),
//   //                              false, false);
//   // 
//   // 
//   // arma::vec mu = arma::vec(alpha.data(), alpha.size(), false, false);
//   // 
//   // int ncols = invB_arma.n_cols;
//   // arma::vec Y = arma::randn(ncols);
//   // return (mu + tinvL_arma * Y);
//   
//   //
//   
//   // Eigen::MatrixXd Lambda_B = Eigen::Map<Eigen::MatrixXd>(Lambda_B_arma.memptr(),
//   //                                                        Lambda_B_arma.n_rows,
//   //                                                        Lambda_B_arma.n_cols);
//   // 
//   // Eigen::VectorXd mu_B = Eigen::Map<Eigen::VectorXd>(mu_B_arma.memptr(),
//   //                                                    mu_B_arma.size());
//   // 
//   // Eigen::MatrixXd L = Lambda_B.llt().matrixL();//rcppeigen_get_chol(Lambda_B);
//   // 
//   // Eigen::VectorXd w = L.triangularView<Eigen::Lower>().solve(mu_B);
//   // 
//   // Eigen::VectorXd mu = L.transpose().triangularView<Eigen::Upper>().solve(w);
//   // 
//   // arma::vec z_arma = arma::randn(L.cols());
//   // 
//   // Eigen::VectorXd z = Eigen::Map<Eigen::VectorXd>(z_arma.memptr(),
//   //                                                 z_arma.size());
//   // 
//   // Eigen::VectorXd v = L.transpose().triangularView<Eigen::Upper>().solve(z);
//   // 
//   // Eigen::VectorXd out = mu + v;
//   // 
//   // arma::vec out_arma = arma::vec(out.data(), out.size(), false, false);
//   // 
//   // return(out_arma);
// }
// 
// double integratedTermz(double i, double c_s, double mu_s, double sigma_s, double omega_s){
//   
//   double k_i = i - .5;
//   
//   double mu_tilde = k_i - omega_s * c_s + (mu_s / (sigma_s * sigma_s));
//   double sigmasq_tilde = omega_s + (1 / (sigma_s * sigma_s));
//   
//   return(exp(k_i * c_s)) * (1 / sigmasq_tilde) * exp(mu_tilde * mu_tilde / (2 * sigmasq_tilde));
// }
// 
// arma::vec sample_z_cpp_integrated(arma::vec c_s, arma::vec p, arma::mat k_s,
//                                   arma::vec mu_clusters, arma::vec sigma_clusters, 
//                                   arma::vec Omega){
//   
//   arma::vec z = arma::zeros(k_s.n_rows);
//   
//   // this loops over p
//   int l = 0;
//   
//   for(int i = 0; (unsigned)i < z.size(); i++){
//     
//     if(k_s(i, 0) == 1){
//       
//       z[i] = 1;
//       l += k_s(i,1);
//       
//     } else {
//       
//       double prod1mp = 1; //prod(1 - p[i + seq_len(k_s[i,4])])
//       for(int k = 0; k < k_s(i,1); k++){
//         prod1mp *= (1 - p(l + k));
//       }
//       
//       double integratedTerm0 = integratedTermz(0, c_s[i], mu_clusters[i],
//                                                sigma_clusters[i], Omega[i]);
//       double integratedTerm1 = integratedTermz(1, c_s[i], mu_clusters[i],
//                                                sigma_clusters[i], Omega[i]);
//       
//       double p_zsequal1 = (prod1mp * integratedTerm1) / (prod1mp * integratedTerm1 + integratedTerm0);
//       z[i] = R::rbinom(1, p_zsequal1);
//       l += k_s(i,1);
//       Rcout << p_zsequal1 << std::endl;
//     }
//     
//   }
//   
//   return(z);
// }  
// 
// arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b, 
//                           arma::vec& Omega, arma::vec& k){
//   
//   arma::mat tX = arma::trans(X);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X;
//   
//   arma::mat L = arma::trans(arma::chol(XtOmegaX + arma::inv(B)));
//   arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   return(result);
// }

// arma::vec sampleNormFast(arma::mat Lambda_B, arma::vec mu_B){
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   
//   arma::vec w = arma::solve(arma::trimatl(L), mu_B);
//   
//   arma::vec mu = arma::solve(arma::trimatu(arma::trans(L)), w);
//   
//   arma::vec z = arma::randn(L.n_cols);
//   
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec out = mu + v;
//   
//   return(out);
// }
// 
// arma::vec sampleNormFast_eigen(arma::mat Lambda_B_arma, arma::vec mu_B_arma){
//   
//   Eigen::MatrixXd Lambda_B = Eigen::Map<Eigen::MatrixXd>(Lambda_B_arma.memptr(),
//                                                          Lambda_B_arma.n_rows,
//                                                          Lambda_B_arma.n_cols);
//   
//   Eigen::VectorXd mu_B = Eigen::Map<Eigen::VectorXd>(mu_B_arma.memptr(),
//                                                      mu_B_arma.size());
//   
//   Eigen::MatrixXd L = Lambda_B.llt().matrixL();//rcppeigen_get_chol(Lambda_B);
//   
//   Eigen::VectorXd w = L.triangularView<Eigen::Lower>().solve(mu_B);
//   
//   Eigen::VectorXd mu = L.transpose().triangularView<Eigen::Upper>().solve(w);
//   
//   arma::vec z_arma = arma::randn(L.cols());
//   
//   Eigen::VectorXd z = Eigen::Map<Eigen::VectorXd>(z_arma.memptr(),
//                                                   z_arma.size());
//   
//   Eigen::VectorXd v = L.transpose().triangularView<Eigen::Upper>().solve(z);
//   
//   Eigen::VectorXd out = mu + v;
//   
//   arma::vec out_arma = arma::vec(out.data(), out.size(), false, false);
//   
//   return(out_arma);
// }
// 
// arma::vec sample_beta_cpp_fast_new(arma::mat& X, arma::mat& B, arma::vec& b, 
//                                    arma::vec& k, arma::mat XtOmegaX){
//   
//   arma::mat tX = arma::trans(X);
//   arma::mat invB = arma::inv(B);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::mat mu_B = tX * k + invB * b;
//   
//   arma::vec result = sampleNormFast_eigen(Lambda_B, mu_B);
//   
//   // arma::mat L = arma::trans(arma::chol(XtOmegaX + arma::inv(B)));
//   // arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
//   // arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   return(result);
// }
// 
// arma::mat matrixProductXtOmegaX_old(arma::mat& X, int Y, int ncov_psi, arma::vec Omega,
//                                     arma::vec X_y_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
//   
//   // first element
//   XtOmegaX2(0, 0) = sum(Omega);
//   
//   // year covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i], X_y_index[i]) += Omega[i];
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + 1 + j, X_y_index[i]) +=  X(i, Y + 1 + j) * Omega[i];
//       
//     }
//   }
//   
//   // transpose
//   for (int i = 1; i <= Y; i++) {
//     
//     XtOmegaX2(0, i) = XtOmegaX2(i, i);
//     XtOmegaX2(i, 0) = XtOmegaX2(i, i);
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(i, Y + 1 + j) = XtOmegaX2(Y + 1 + j, i);
//       
//     }
//     
//   }
//   
//   for(int i = 0; i < ncov_psi; i++){
//     
//     double result = arma::as_scalar(Omega.t() * X.col(Y + 1 + i));
//     
//     XtOmegaX2(Y + 1 + i, 0) = result;
//     XtOmegaX2(0, Y + 1 + i) = XtOmegaX2(Y + 1 + i, 0);
//   }
//   
//   for(int i = 0; i < ncov_psi; i++){
//     for (int j = 0; j <= i; j++) {
//       arma::vec firstProduct = Omega % X.col(Y + 1 + i);
//       arma::vec secondProduct = firstProduct % X.col(Y + 1 + j);
//       XtOmegaX2(Y + 1 + i, Y + 1 + j) = sum(secondProduct);
//     }
//   }
//   
//   for (int i = 0; i < (ncov_psi- 1); i++) {
//     for (int j = i; j < ncov_psi; j++) {
//       XtOmegaX2(Y + 1 + i, Y + 1 + j) = XtOmegaX2(Y + 1 + j, Y + 1 + i);
//     }
//   }
//   
//   return(XtOmegaX2);
// }
// 
// 
// arma::vec XtransposeKarma_old(arma::mat &X, IntegerVector X_y_index, 
//                               IntegerVector X_s_index, arma::vec &k, 
//                               int Y, int centers, int ncov){
//   
//   
//   arma::vec Xk = arma::zeros(X.n_cols);
//   
//   Xk[0] = as_scalar(k.t() * X.col(0));
//   for(int i = 0; i < X_y_index.size(); i++){
//     Xk(X_y_index[i]) += k[i];
//   }
//   for(int i = 0; i < X_s_index.size(); i++){
//     Xk(Y + X_s_index[i]) += k[i];
//   }
//   for(int i = 0; i < ncov; i++){
//     Xk(1 + Y + centers + i) = as_scalar(k.t() * X.col(1 + Y + centers + i));
//   }
//   
//   return(Xk);
// }
// 
// Eigen::VectorXd XtransposeK(Eigen::MatrixXd &X, IntegerVector X_y_index, 
//                             IntegerVector X_s_index, Eigen::VectorXd &k, 
//                             int Y, int centers, int ncov){
//   
//   
//   VectorXd Xk = VectorXd(X.cols());
//   Xk.setZero();
//   
//   Xk(0) = X.col(0).adjoint() * k;
//   for(int i = 0; i < X_y_index.size(); i++){
//     Xk(X_y_index[i]) += k[i];
//   }
//   for(int i = 0; i < X_s_index.size(); i++){
//     Xk(Y + X_s_index[i]) += k[i];
//   }
//   for(int i = 0; i < ncov; i++){
//     Xk(1 + Y + centers + i) = X.col(1 + Y + centers + i).adjoint() * k;
//   }
//   
//   return(Xk);
// }
// 
// List sampler_beta_nonspatial(arma::vec beta,
//                              arma::vec a_s,
//                              arma::mat X, 
//                              arma::vec b, 
//                              arma::mat invB, 
//                              arma::vec n, 
//                              arma::vec k, 
//                              int Y, 
//                              int ncov_psi,
//                              arma::vec X_y_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, ncov_psi + 5, Omega,
//                                              X_y_index);
//   
//   // arma::mat tX = arma::trans(X);
//   // arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   // arma::mat XtOmegaX = tXOmega * X;
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast(X, invB, b, knew, XtOmegaX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega,
//                       _["XtOmegaX"] = XtOmegaX));
// }
// 
// List sampler_beta_sp(arma::vec beta,
//                      arma::vec a_s,
//                      arma::mat& X, 
//                      arma::vec b, 
//                      arma::mat invB, 
//                      arma::vec n, 
//                      arma::vec k,
//                      int Y, 
//                      int X_centers,
//                      int ncov_psi,
//                      int numTimeSpaceCov,
//                      IntegerVector X_y_index,
//                      IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_sp(X, Y, X_centers, ncov_psi + numTimeSpaceCov, Omega,
//                                                 X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
//                                      X_y_index, X_s_index, Y, X_centers, ncov_psi + numTimeSpaceCov);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// 
// double prior_as(double a_s, double sum_as, double beta, double a){
//   return(pow(a_s*a_s + sum_as + beta,-a));
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_as_mh_t_cpp(arma::vec a_s, arma::vec k_s, arma::vec sites, arma::vec Xbeta,
//                              arma::vec z, double a_sigma, double b_sigma){
//   
//   
//   int index_site = 0;
//   
//   double sum_as = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int l = 1;
//     
//     int site = sites[i];
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;  
//     
//     sum_as += a_s[indexes_site[0]] * a_s[indexes_site[0]];
//     
//   }
//   
//   index_site = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int site = sites[i];
//     
//     int l = 1;
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;  
//     
//     arma::vec c_site = arma::zeros(l); 
//     arma::vec z_site = arma::zeros(l);
//     for(int j = 0; j < l; j++){
//       c_site[j] = Xbeta[indexes_site[j]];
//       z_site[j] = z[indexes_site[j]];
//     }  
//     
//     // sample x
//     arma::vec k_site = (z_site - .5);
//     arma::vec n_site = arma::ones(l);
//     
//     double a_s_old = a_s[indexes_site[0]];
//     double a_s_new = R::rnorm(a_s_old, .05);
//     
//     double loglikelihood_current = 0;
//     for(int j = 0; j < l; j++){
//       loglikelihood_current += z[indexes_site[j]] * (Xbeta[indexes_site[j]] + a_s_old) - 
//         log(1 + exp(Xbeta[indexes_site[j]] + a_s_old));
//     }
//     
//     double loglikelihood_new = 0;
//     for(int j = 0; j < l; j++){
//       loglikelihood_new += z[indexes_site[j]] * (Xbeta[indexes_site[j]] + a_s_new) - 
//         log(1 + exp(Xbeta[indexes_site[j]] + a_s_new));
//     }
//     
//     sum_as -= (a_s_old * a_s_old);
//     
//     double prior_current = prior_as(a_s_old, sum_as, b_sigma, (a_sigma + sites.size())/ 2);
//     double prior_new = prior_as(a_s_new, sum_as, b_sigma, (a_sigma + sites.size())/ 2);
//     
//     double logposterior_current = dt_cpp(a_s_old, 2 * a_sigma, 0, b_sigma / a_sigma, 1) +
//       loglikelihood_current;
//     double logposterior_new = dt_cpp(a_s_new, 2 * a_sigma, 0, b_sigma / a_sigma, 1) +
//       loglikelihood_new;
//     
//     if(R::runif(0,1) < exp(logposterior_new - logposterior_current)){
//       a_s_old = a_s_new;
//     }
//     
//     sum_as += (a_s_old * a_s_old);
//     
//     for(int j = 0; j < l; j++){
//       a_s[indexes_site[j]] = a_s_old;
//     }
//     
//     // mu_clusters[site - 1] = b * a;
//     // sigma_clusters[site - 1] = sqrt(a);
//     
//   }
//   
//   return(a_s);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_as_cpp_old(arma::vec k_s, arma::vec sites,
//                             arma::vec beta_psi, arma::mat &X_psi,
//                             arma::vec k, 
//                             arma::vec z, arma::vec Omega,
//                             double sigma_a){
//   
//   arma::vec c_psi = X_psi * beta_psi;
//   arma::vec b_psi = k - Omega % c_psi;
//   
//   arma::vec a_s = arma::zeros(k_s.size());
//   
//   int index_site = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int site = sites[i];
//     
//     int l = 1;
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;
//     
//     // arma::mat data_psi_site(l, X_psi.n_cols);
//     // arma::vec z_site = arma::zeros(l);
//     arma::vec w_site = arma::zeros(l);
//     // arma::vec c_site = arma::zeros(l);
//     arma::vec b_site = arma::zeros(l);
//     for(int j = 0; j < l; j++){
//       // data_psi_site.row(j) = X_psi.row(indexes_site[j]);
//       // z_site[j] = z[indexes_site[j]];
//       w_site[j] = Omega[indexes_site[j]];
//       // c_site[j] = c_psi[indexes_site[j]];
//       b_site[j] = b_psi[indexes_site[j]];
//     }
//     
//     // sample x
//     // arma::vec k_site = (z_site - .5);
//     // arma::vec n_site = arma::ones(l);
//     // arma::vec c_site = data_psi_site * beta_psi;
//     
//     double a = 1 / (sum(w_site) + 1 / (sigma_a*sigma_a));
//     // double b = sum(k_site - w_site % c_site);
//     double b = sum(b_site);
//     
//     double x = R::rnorm(b * a, sqrt(a));
//     
//     for(int j = 0; j < l; j++){
//       a_s[indexes_site[j]] = x;
//     }
//     
//   }
//   
//   return(a_s);
// }
// 
// // [[Rcpp::export]]
// double dt_cpp(double x, double nu, double mu, double sigma, bool returnLog){
//   
//   double ratio = exp(R::lgammafn(( nu + 1 ) / 2) - R::lgammafn( nu / 2 ));
//   double product = (x - mu) * (x - mu) / (sigma * sigma);
//   double num = pow(1 + (1 / nu) * product, - ( nu + 1 ) / 2);
//   double den = sigma * sqrt(M_PI * nu);
//   if(returnLog){
//     return log(ratio * num / den);
//   } else {
//     return ratio * num / den;
//   }
// }
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector arma_setdiff(arma::uvec& x, arma::uvec& y){
//   
//   x = arma::unique(x);
//   y = arma::unique(y);
//   for (size_t j = 0; j < y.n_elem; j++) {
//     arma::uvec q1 = arma::find(x == y[j]);
//     x.shed_row(q1(0));
//   }
//   
//   Rcpp::NumericVector x2 = Rcpp::wrap(x);
//   x2.attr("dim") = R_NilValue;
//   return x2;
// }
// 
// // [[Rcpp::export]]
// IntegerVector stl_sort(IntegerVector x) {
//   IntegerVector y = clone(x);
//   std::sort(y.begin(), y.end());
//   return y;
// }
// 
// // [[Rcpp::export]]
// IntegerVector setdiffna(IntegerVector S_c){
//   IntegerVector S_cnew = stl_sort(unique(S_c));
//   IntegerVector NAvec = IntegerVector::create(NA_INTEGER);
//   arma::uvec temp1 = as<arma::uvec>(S_cnew);
//   arma::uvec temp2 = as<arma::uvec>(NAvec);
//   IntegerVector temp = as<IntegerVector>(arma_setdiff(temp1, temp2));
//   return temp;
// }
// 
// 
// 
// arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
//   int ncols = cholsigma.n_cols;
//   arma::vec Y = arma::randn(ncols);
//   return mu + cholsigma * Y;
// }
// 
// 
// // [[Rcpp::export]]
// double loglikelihood_l_gp_cpp(double l, 
//                               double a, 
//                               int Y, 
//                               arma::mat XtOmegaX, 
//                               arma::mat Xtz){
//   
//   arma::vec years(Y);
//   std::iota(years.begin(), years.end(), 1);
//   arma::mat K_l = K(years, years, a, l);
//   
//   arma::mat invKl = arma::inv(K_l);
//   
//   arma::mat XtOmegaXpsolveK_l = XtOmegaX + invKl;
//   arma::mat identityMatrix = arma::eye(invKl.n_rows, invKl.n_rows);
//   arma::mat invKlXtOmegaX = invKl * XtOmegaX;
//   arma::mat IpXtOmegaXpsolveK_l = identityMatrix + invKlXtOmegaX;
//   
//   double logdet1det2;
//   double sign;
//   
//   log_det(logdet1det2, sign, IpXtOmegaXpsolveK_l);
//   
//   arma::vec term = arma::trans(Xtz) * arma::inv(XtOmegaXpsolveK_l) * Xtz;
//   
//   return (- .5 * logdet1det2 + .5 * term[0]);
// }
// 
// 
// // [[Rcpp::export]]
// arma::mat vec_subset_mat(const arma::mat& x, const arma::uvec& idx) {
//   return x.cols(idx);
// }
// 
// // [[Rcpp::export]]
// arma::vec subset_vec(const arma::vec& x, const arma::uvec& idx) {
//   return x.elem(idx);
// }
// 
// // [[Rcpp::export]]
// double sampler_l(double l, double a, arma::vec beta, arma::mat X,
//                  double a_l, double b_l, double sd_l,
//                  arma::vec a_s, int Y,
//                  arma::mat XtOmegaX, 
//                  arma::vec k, arma::vec Omega,
//                  arma::vec X_y_index){
//   
//   arma::vec idxes = arma::zeros(X.n_cols - Y);
//   idxes[0] = 0;
//   for(int i = 0; (unsigned)i < (idxes.size() - 1); i++) idxes[i + 1] = 1 + Y + i;
//   arma::uvec uidxes = arma::conv_to<arma::uvec>::from(idxes);
//   arma::mat X_no_y = vec_subset_mat(X, uidxes);
//   arma::vec beta_no_y = subset_vec(beta, uidxes);
//   
//   arma::vec c_i = X_no_y * beta_no_y + a_s;
//   
//   arma::vec y = k - Omega % c_i;
//   
//   arma::vec Xtz = arma::zeros(Y);
//   for(int i = 0; (unsigned)i < X_y_index.size(); i++){
//     Xtz[X_y_index[i] - 1] += y[i];
//   }
//   
//   double l_star = R::rnorm(l, sd_l);
//   
//   arma::mat XtOmegaX_subset = arma::zeros(Y, Y);
//   for(int i = 0; i < Y; i++){
//     XtOmegaX_subset(i,i) = XtOmegaX(i + 1, i + 1);
//   }
//   
//   if(l_star > 0){
//     
//     double loglikelihood_star = loglikelihood_l_gp_cpp(l_star, a, Y, XtOmegaX_subset, Xtz);
//     double loglikelihood_current = loglikelihood_l_gp_cpp(l, a, Y, XtOmegaX_subset, Xtz);
//     
//     double logprior_star = R::dgamma(l_star, a_l, 1 / b_l, 1);
//     double logprior_current = R::dgamma(l, a_l, 1 / b_l, 1);
//     
//     double logposterior_star = loglikelihood_star + logprior_star;
//     double logposterior_current = loglikelihood_current + logprior_current;
//     
//     if(R::runif(0, 1) < exp(logposterior_star - logposterior_current)){
//       
//       l = l_star;
//       
//     }
//     
//   }
//   
//   return l;
// }
//
// // [[Rcpp::depends(RcppParallel)]]
// #include <RcppParallel.h>
// using namespace RcppParallel;
//
//
// #include <RcppEigen.h>
// 
// // [[Rcpp::depends(RcppEigen)]]
//
//
// struct Sampleomega : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomega(NumericMatrix X, NumericVector beta, NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// }; 
// 
// NumericVector sampleOmegaParallel(NumericMatrix &X, NumericVector beta, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomega sampleomega(X, beta, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
// 
// struct Sampleomegaas : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> as;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomegaas(NumericMatrix X, NumericVector beta, NumericVector as,
//                 NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), as(as), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       Xibeta += as[i];
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// };
// 
// NumericVector sampleOmegaParallelas(NumericMatrix &X, NumericVector beta, 
//                                     NumericVector as, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomegaas sampleomega(X, beta, as, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
//
// 
// arma::vec XtransposeKarma_SoR(arma::mat &X, IntegerVector X_y_index, 
//                               arma::mat &X_s_index, arma::vec &k, 
//                               int Y, int centers, int ncov){
//   
//   
//   arma::vec Xk = arma::zeros(X.n_cols);
//   
//   for(int i = 0; i < X_y_index.size(); i++){
//     Xk(X_y_index[i] - 1) += k[i];
//   }
//   
//   for(int i = 0; i < k.size(); i++){
//     for(int l = 0; l < X_s_index.n_cols; l++){
//       Xk(Y + X_s_index(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * k[i];
//     }
//   }
//   
//   for(int i = 0; i < ncov; i++){
//     Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
//   }
//   
//   return(Xk);
// }
// 
// arma::vec sample_beta_cpp_fast_sparse_SoR(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                           arma::vec &k, arma::mat XtOmegaX,
//                                           IntegerVector X_y_index, arma::mat &X_s_index,
//                                           int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma_SoR(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   return(result);
// }
//
//
// struct Sampleomega : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomega(NumericMatrix X, NumericVector beta, NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// }; 
// 
// NumericVector sampleOmegaParallel(NumericMatrix &X, NumericVector beta, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomega sampleomega(X, beta, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
// 
// struct Sampleomegaas : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> as;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomegaas(NumericMatrix X, NumericVector beta, NumericVector as,
//                 NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), as(as), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       Xibeta += as[i];
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// };
// 
// NumericVector sampleOmegaParallelas(NumericMatrix &X, NumericVector beta, 
//                                     NumericVector as, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomegaas sampleomega(X, beta, as, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
//
// 
// arma::mat matrixProductXtOmegaX_SoR(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
//                                     IntegerVector X_y_index, arma::mat X_s_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
//   
//   int maxPoints = X_s_index.n_cols;
//   
//   // year covariates times year covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
//     
//   }
//   
//   // spatial  covariates times spatial covariates
//   for(int l = 0; l < maxPoints; l++){
//     
//     for(int l2 = 0; l2 < maxPoints; l2++){
//       
//       for (int i = 0; (unsigned)i < Omega.size(); i++){
//         
//         XtOmegaX2(Y + X_s_index(i,l) - 1, Y + X_s_index(i,l2) - 1) += X(i, Y + X_s_index(i, l) - 1) * 
//           Omega[i] * X(i, Y + X_s_index(i, l2) - 1);    
//         
//       }
//       
//       // XtOmegaX2(Y + X_s_index[i,l2] - 1, Y + X_s_index[i,l] - 1) = 
//       // XtOmegaX2(Y + X_s_index[i,l] - 1, Y + X_s_index[i,l2] - 1);
//       
//     }
//     
//   }
//   
//   // year covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
//       
//     }
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
//     }
//     
//   }
//   
//   // spatial  covariates times year covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int l = 0; l < maxPoints; l++){
//       
//       XtOmegaX2(X_y_index[i] - 1, Y + X_s_index(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//       
//     }
//     
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < X_centers; j++){
//       XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
//     }
//     
//   }
//   
//   // spatial covariates times standard covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       for(int l = 0; l < maxPoints; l++){
//         
//         XtOmegaX2(Y + X_centers + j, Y + X_s_index(i,l) - 1) +=  
//           X(i, Y + X_centers + j) * X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//         
//       }
//     }
//   }
//   
//   for (int i = 1; i <= X_centers; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
//     }
//     
//   }
//   
//   // standard covariates times standard covariates 
//   
//   for(int i = 0; i < ncov_psi; i++){
//     for (int j = 0; j <= i; j++) {
//       // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
//       // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
//       // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
//       for(int l = 0; l < Omega.size(); l++){
//         XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + X_centers + i) * X(l, Y + X_centers + j);
//       }
//     }
//   }
//   
//   for (int i = 0; i < (ncov_psi- 1); i++) {
//     for (int j = i; j < ncov_psi; j++) {
//       XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
//         XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
//     }
//   }
//   
//   return(XtOmegaX2);
// }
//
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp(arma::vec beta,
//                            arma::mat &X, 
//                            arma::vec b, 
//                            arma::mat B, 
//                            arma::vec n, 
//                            arma::vec k){
//   
//   // sample Omega
//   arma::vec Omega = sample_Omega_cpp(X, beta, n);
//   
//   arma::mat tX = arma::trans(X);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X;
//   
//   beta = sample_beta_cpp_fast(X, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));
// }
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp_parallel(SEXP beta,
//                                     SEXP X, 
//                                     arma::vec b, 
//                                     arma::mat B, 
//                                     SEXP n, 
//                                     arma::vec k){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallel(X_num, beta_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   
//   arma::mat tX = arma::trans(X_arma);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast(X_arma, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));
// }
// 
// arma::vec sample_Omega_cpp(arma::mat &X, arma::vec &beta, arma::vec &n){
//   
//   int nsize = n.size();
//   arma::vec Omega_vec(nsize);
//   
//   arma::vec Xbeta = X * beta;
//   
//   for(int i = 0; i < nsize; i++){
//     
//     Omega_vec[i] = rpg(n[i], Xbeta[i]);
//     
//   }
//   
//   return(Omega_vec);
// }
//
// 
// arma::mat diagMatrixProd(arma::mat &X, arma::vec &D){
//   // this is slow
//   arma::mat result(X.n_rows, D.size());
//   for(int j = 0; (unsigned)j < result.n_cols; j++){
//     result.col(j) = X.col(j) * D(j);
//   }
//   
//   // RMatrix<double>::Row rowi = A.row(i);
//   
//   // out(i,1) = std::inner_product(rowi.begin(), rowi.end(), x.begin(), 0.0);
//   
//   return(result);
// }
//
// arma::vec sample_beta_cpp_fast_sparse_constrained(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                                   arma::vec &k, arma::mat XtOmegaX,
//                                                   IntegerVector X_y_index, IntegerVector X_s_index,
//                                                   int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // constraint
//   
//   // arma::mat invQ = Lambda_B;
//   arma::mat invQ = arma::inv(Lambda_B);
//   
//   // sQ <- sum(Lambda_B[Y + 1:X_centers, Y + 1:X_centers])
//   double sQ = 0;
//   for(int l = 0; l < centers; l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       sQ += invQ(Y + l, Y + l2);
//     }
//   }
//   
//   // columns <- apply(Lambda_B[,Y + 1:X_centers], 1, sum) / sQ
//   arma::vec columns = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       columns[l] += invQ(l, Y + l2);
//     }
//     columns[l] = columns[l] / sQ;
//   }
//   
//   // allTerm <- sapply(1:length(x), function(i){
//   //sum(columns[i] * x[Y + 1:X_centers])
//   //})
//   arma::vec secondTerm = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       secondTerm[l] += columns[l] * result[Y + l2];
//     }
//   }
//   
//   // Rcout << columns << " - " << sQ << std::endl;
//   
//   result = result - secondTerm;
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// List sampler_beta(arma::vec beta,
//                   arma::vec a_s,
//                   arma::mat &X, 
//                   arma::vec b, 
//                   arma::mat invB, 
//                   arma::vec n, 
//                   arma::vec k,
//                   int Y, 
//                   int X_centers,
//                   int ncov_psi,
//                   int numTimeSpaceCov,
//                   IntegerVector X_y_index,
//                   IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, X_centers, ncov_psi + numTimeSpaceCov, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
//                                      X_y_index, X_s_index, Y, X_centers, ncov_psi + numTimeSpaceCov);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));
// }
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel(SEXP beta,
//                            SEXP a_s,
//                            SEXP &X, 
//                            arma::vec b, 
//                            arma::mat invB, 
//                            SEXP n, 
//                            arma::vec k,
//                            int Y, 
//                            int X_centers,
//                            int ncov_psi,
//                            int numTimeSpaceCov,
//                            IntegerVector X_y_index,
//                            IntegerVector X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   // arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//   // X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//                                                   X_y_index, X_s_index, Y, X_centers, 
//                                                   numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
//
// 
// arma::vec sample_beta_cpp_fast(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                arma::vec &k, arma::mat XtOmegaX){
//   
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = arma::trans(X) * k + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   return(result);
// }
//
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel_SoR(SEXP beta,
//                                SEXP a_s,
//                                SEXP &X, 
//                                arma::vec b, 
//                                arma::mat invB, 
//                                SEXP n, 
//                                arma::vec k,
//                                int Y, 
//                                int X_centers,
//                                int ncov_psi,
//                                int numTimeSpaceCov,
//                                IntegerVector X_y_index,
//                                arma::mat &X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_SoR(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                                  X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast_sparse_SoR(X_arma, invB, b, knew, XtOmegaX,
//                                                       X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * k[i];
//     }
//   }
//   
//   for(int i = 0; i < ncov; i++){
//     Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
//   }
//   
//   return(Xk);
// }
// 
// arma::vec sample_beta_cpp_fast_sparse_SoR(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                           arma::vec &k, arma::mat XtOmegaX,
//                                           IntegerVector X_y_index, arma::mat &X_s_index,
//                                           int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma_SoR(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   return(result);
// }
//
//
// struct Sampleomega : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomega(NumericMatrix X, NumericVector beta, NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// }; 
// 
// NumericVector sampleOmegaParallel(NumericMatrix &X, NumericVector beta, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomega sampleomega(X, beta, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
// 
// struct Sampleomegaas : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> as;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomegaas(NumericMatrix X, NumericVector beta, NumericVector as,
//                 NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), as(as), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       Xibeta += as[i];
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// };
// 
// NumericVector sampleOmegaParallelas(NumericMatrix &X, NumericVector beta, 
//                                     NumericVector as, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomegaas sampleomega(X, beta, as, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
//
// 
// arma::mat matrixProductXtOmegaX_SoR(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
//                                     IntegerVector X_y_index, arma::mat X_s_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
//   
//   int maxPoints = X_s_index.n_cols;
//   
//   // year covariates times year covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
//     
//   }
//   
//   // spatial  covariates times spatial covariates
//   for(int l = 0; l < maxPoints; l++){
//     
//     for(int l2 = 0; l2 < maxPoints; l2++){
//       
//       for (int i = 0; (unsigned)i < Omega.size(); i++){
//         
//         XtOmegaX2(Y + X_s_index(i,l) - 1, Y + X_s_index(i,l2) - 1) += X(i, Y + X_s_index(i, l) - 1) * 
//           Omega[i] * X(i, Y + X_s_index(i, l2) - 1);    
//         
//       }
//       
//       // XtOmegaX2(Y + X_s_index[i,l2] - 1, Y + X_s_index[i,l] - 1) = 
//       // XtOmegaX2(Y + X_s_index[i,l] - 1, Y + X_s_index[i,l2] - 1);
//       
//     }
//     
//   }
//   
//   // year covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
//       
//     }
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
//     }
//     
//   }
//   
//   // spatial  covariates times year covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int l = 0; l < maxPoints; l++){
//       
//       XtOmegaX2(X_y_index[i] - 1, Y + X_s_index(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//       
//     }
//     
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < X_centers; j++){
//       XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
//     }
//     
//   }
//   
//   // spatial covariates times standard covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       for(int l = 0; l < maxPoints; l++){
//         
//         XtOmegaX2(Y + X_centers + j, Y + X_s_index(i,l) - 1) +=  
//           X(i, Y + X_centers + j) * X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//         
//       }
//     }
//   }
//   
//   for (int i = 1; i <= X_centers; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
//     }
//     
//   }
//   
//   // standard covariates times standard covariates 
//   
//   for(int i = 0; i < ncov_psi; i++){
//     for (int j = 0; j <= i; j++) {
//       // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
//       // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
//       // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
//       for(int l = 0; l < Omega.size(); l++){
//         XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + X_centers + i) * X(l, Y + X_centers + j);
//       }
//     }
//   }
//   
//   for (int i = 0; i < (ncov_psi- 1); i++) {
//     for (int j = i; j < ncov_psi; j++) {
//       XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
//         XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
//     }
//   }
//   
//   return(XtOmegaX2);
// }
//
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp(arma::vec beta,
//                            arma::mat &X, 
//                            arma::vec b, 
//                            arma::mat B, 
//                            arma::vec n, 
//                            arma::vec k){
//   
//   // sample Omega
//   arma::vec Omega = sample_Omega_cpp(X, beta, n);
//   
//   arma::mat tX = arma::trans(X);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X;
//   
//   beta = sample_beta_cpp_fast(X, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));
// }
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp_parallel(SEXP beta,
//                                     SEXP X, 
//                                     arma::vec b, 
//                                     arma::mat B, 
//                                     SEXP n, 
//                                     arma::vec k){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallel(X_num, beta_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   
//   arma::mat tX = arma::trans(X_arma);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast(X_arma, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));
// }
// 
// arma::vec sample_Omega_cpp(arma::mat &X, arma::vec &beta, arma::vec &n){
//   
//   int nsize = n.size();
//   arma::vec Omega_vec(nsize);
//   
//   arma::vec Xbeta = X * beta;
//   
//   for(int i = 0; i < nsize; i++){
//     
//     Omega_vec[i] = rpg(n[i], Xbeta[i]);
//     
//   }
//   
//   return(Omega_vec);
// }
//
// 
// arma::mat diagMatrixProd(arma::mat &X, arma::vec &D){
//   // this is slow
//   arma::mat result(X.n_rows, D.size());
//   for(int j = 0; (unsigned)j < result.n_cols; j++){
//     result.col(j) = X.col(j) * D(j);
//   }
//   
//   // RMatrix<double>::Row rowi = A.row(i);
//   
//   // out(i,1) = std::inner_product(rowi.begin(), rowi.end(), x.begin(), 0.0);
//   
//   return(result);
// }
//
// arma::vec sample_beta_cpp_fast_sparse_constrained(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                                   arma::vec &k, arma::mat XtOmegaX,
//                                                   IntegerVector X_y_index, IntegerVector X_s_index,
//                                                   int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // constraint
//   
//   // arma::mat invQ = Lambda_B;
//   arma::mat invQ = arma::inv(Lambda_B);
//   
//   // sQ <- sum(Lambda_B[Y + 1:X_centers, Y + 1:X_centers])
//   double sQ = 0;
//   for(int l = 0; l < centers; l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       sQ += invQ(Y + l, Y + l2);
//     }
//   }
//   
//   // columns <- apply(Lambda_B[,Y + 1:X_centers], 1, sum) / sQ
//   arma::vec columns = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       columns[l] += invQ(l, Y + l2);
//     }
//     columns[l] = columns[l] / sQ;
//   }
//   
//   // allTerm <- sapply(1:length(x), function(i){
//   //sum(columns[i] * x[Y + 1:X_centers])
//   //})
//   arma::vec secondTerm = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       secondTerm[l] += columns[l] * result[Y + l2];
//     }
//   }
//   
//   // Rcout << columns << " - " << sQ << std::endl;
//   
//   result = result - secondTerm;
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// List sampler_beta(arma::vec beta,
//                   arma::vec a_s,
//                   arma::mat &X, 
//                   arma::vec b, 
//                   arma::mat invB, 
//                   arma::vec n, 
//                   arma::vec k,
//                   int Y, 
//                   int X_centers,
//                   int ncov_psi,
//                   int numTimeSpaceCov,
//                   IntegerVector X_y_index,
//                   IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, X_centers, ncov_psi + numTimeSpaceCov, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
//                                      X_y_index, X_s_index, Y, X_centers, ncov_psi + numTimeSpaceCov);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));
// }
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel(SEXP beta,
//                            SEXP a_s,
//                            SEXP &X, 
//                            arma::vec b, 
//                            arma::mat invB, 
//                            SEXP n, 
//                            arma::vec k,
//                            int Y, 
//                            int X_centers,
//                            int ncov_psi,
//                            int numTimeSpaceCov,
//                            IntegerVector X_y_index,
//                            IntegerVector X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   // arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//   // X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//                                                   X_y_index, X_s_index, Y, X_centers, 
//                                                   numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
//
// 
// arma::vec sample_beta_cpp_fast(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                arma::vec &k, arma::mat XtOmegaX){
//   
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = arma::trans(X) * k + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   return(result);
// }
//
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel_SoR(SEXP beta,
//                                SEXP a_s,
//                                SEXP &X, 
//                                arma::vec b, 
//                                arma::mat invB, 
//                                SEXP n, 
//                                arma::vec k,
//                                int Y, 
//                                int X_centers,
//                                int ncov_psi,
//                                int numTimeSpaceCov,
//                                IntegerVector X_y_index,
//                                arma::mat &X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_SoR(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                                  X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast_sparse_SoR(X_arma, invB, b, knew, XtOmegaX,
//                                                       X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }t << columns << " - " << sQ << std::endl;
                            //   
                              //   result = result - secondTerm;
                              //   
                                //   return(result);
                              // }
  // 
    // // [[Rcpp::export]]
  // List sampler_beta(arma::vec beta,
                       //                   arma::vec a_s,
                       //                   arma::mat &X, 
                       //                   arma::vec b, 
                       //                   arma::mat invB, 
                       //                   arma::vec n, 
                       //                   arma::vec k,
                       //                   int Y, 
                       //                   int X_centers,
                       //                   int ncov_psi,
                       //                   int numTimeSpaceCov,
                       //                   IntegerVector X_y_index,
                       //                   IntegerVector X_s_index){
    //   
      //   arma::vec Xbeta = X * beta + a_s;
      //   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
      //   
        //   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, X_centers, ncov_psi + numTimeSpaceCov, Omega,
                                                        //                                              X_y_index, X_s_index);
        //   
          //   arma::vec knew = k - Omega % a_s;
          //   
            //   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
                                                    //                                      X_y_index, X_s_index, Y, X_centers, ncov_psi + numTimeSpaceCov);
            //   
              //   return(List::create(_["beta"] = beta,
                                        //                       _["Omega"] = Omega));
            // }
  // 
    // // [[Rcpp::export]]
  // List sampler_beta_parallel(SEXP beta,
                                //                            SEXP a_s,
                                //                            SEXP &X, 
                                //                            arma::vec b, 
                                //                            arma::mat invB, 
                                //                            SEXP n, 
                                //                            arma::vec k,
                                //                            int Y, 
                                //                            int X_centers,
                                //                            int ncov_psi,
                                //                            int numTimeSpaceCov,
                                //                            IntegerVector X_y_index,
                                //                            IntegerVector X_s_index){
    //   
      //   NumericMatrix X_num = as<NumericMatrix>(X);
      //   NumericVector beta_num = as<NumericVector>(beta);
      //   NumericVector n_num = as<NumericVector>(n);
      //   NumericVector as_num = as<NumericVector>(a_s);
      //   
        //   // sample Omega
      //   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
      //   
        //   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
      //   arma::vec as_arma(as_num.begin(), as_num.size(), false);
      //   
        //   // arma::vec Xbeta = X * beta + a_s;
        //   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
        //   
          //   arma::mat XtOmegaX = matrixProductXtOmegaX(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
                                                          //                                              X_y_index, X_s_index);
          //   
            //   arma::vec knew = k - Omega % as_arma;
            //   
              //   // arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
                                                                      //   // X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
              //   arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
                                                                   //                                                   X_y_index, X_s_index, Y, X_centers, 
                                                                   //                                                   numTimeSpaceCov + ncov_psi);
              //   
                //   return(List::create(_["beta"] = betaout,
                                          //                       _["Omega"] = Omega));//,
              //   // _["XtOmegaX"] = XtOmegaX));
// }
  //
    // 
    // arma::vec sample_beta_cpp_fast(arma::mat &X, arma::mat &invB, arma::vec &b, 
                                      //                                arma::vec &k, arma::mat XtOmegaX){
      //   
        //   
        //   arma::mat Lambda_B = XtOmegaX + invB;
        //   arma::vec mu_B = arma::trans(X) * k + invB * b;
        //   
          //   arma::mat L = arma::trans(arma::chol(Lambda_B));
          //   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
          //   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
          //   
            //   arma::vec z = arma::randn(invB.n_cols);
            //   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
            //   
              //   arma::vec result = v + alpha;
              //   
                //   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
                //   
                  //   return(result);
                // }
  //
    // 
    // // [[Rcpp::export]]
  // List sampler_beta_parallel_SoR(SEXP beta,
                                    //                                SEXP a_s,
                                    //                                SEXP &X, 
                                    //                                arma::vec b, 
                                    //                                arma::mat invB, 
                                    //                                SEXP n, 
                                    //                                arma::vec k,
                                    //                                int Y, 
                                    //                                int X_centers,
                                    //                                int ncov_psi,
                                    //                                int numTimeSpaceCov,
                                    //                                IntegerVector X_y_index,
                                    //                                arma::mat &X_s_index){
    //   
      //   NumericMatrix X_num = as<NumericMatrix>(X);
      //   NumericVector beta_num = as<NumericVector>(beta);
      //   NumericVector n_num = as<NumericVector>(n);
      //   NumericVector as_num = as<NumericVector>(a_s);
      //   
        //   // sample Omega
      //   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
      //   
        //   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
      //   arma::vec as_arma(as_num.begin(), as_num.size(), false);
      //   
        //   // arma::vec Xbeta = X * beta + a_s;
        //   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
        //   
          //   arma::mat XtOmegaX = matrixProductXtOmegaX_SoR(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
                                                              //                                                  X_y_index, X_s_index);
          //   
            //   arma::vec knew = k - Omega % as_arma;
            //   
              //   arma::vec betaout = sample_beta_cpp_fast_sparse_SoR(X_arma, invB, b, knew, XtOmegaX,
                                                                       //                                                       X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
              //   
                //   return(List::create(_["beta"] = betaout,
                                          //                       _["Omega"] = Omega));//,
              //   // _["XtOmegaX"] = XtOmegaX));
// }