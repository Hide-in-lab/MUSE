#include <iostream>
#include <RcppArmadillo.h>
#include <armadillo>
#include <Rcpp.h>
#include <boost/random.hpp>
#include <random>
#include <vector>
#include <thread>
#include <mutex>
#include <ctime>
#include <cmath>
#include <RcppEigen.h>
#include <Eigen/Dense>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

NumericVector rnorm_cpp(int n, double mean, double sd);
  
double dnorm_cpp(double x, double mean, double sd); 
  
NumericMatrix mvrnorm_cpp(int n, NumericVector mean, NumericMatrix cov_matrix);
  
double runi_cpp(double min_val, double max_val); 
  
NumericVector rinvgamma_cpp(int n, double shape, double scale);
  
CharacterVector rsample(CharacterVector characters, NumericVector probs, int sample_size);
  
NumericMatrix rDirichlet( NumericVector alpha, int times );
  
List lm_cpp(NumericVector y, NumericMatrix X);
   
NumericMatrix cor_cpp(NumericMatrix input_matrix);
  
IntegerVector whichIn(CharacterVector x, CharacterVector values);
  
double vtv(NumericVector x, NumericVector y = NumericVector());
  
NumericVector qt_cpp(NumericVector x, NumericVector probs);
  
NumericVector res_output(NumericVector x);
   
NumericVector index_cpp(NumericVector x , IntegerVector index = IntegerVector());
  
NumericVector minor_genotype(NumericVector x);
  
List genotype_generation(List size_sample);

IntegerVector order_abs(NumericVector x, bool decreasing = true);

NumericVector rmindex(NumericVector x, IntegerVector indices = IntegerVector::create(0));

NumericVector arm_to_num(arma::vec x);

// [[Rcpp::export]]
List MUSE( arma::vec &data_mean_gamma, arma::vec &data_s2_gamma, arma::vec &data_mean_Gamma, arma::vec &data_s2_Gamma, 
           int n_snps, int iter_times = 1000, Rcpp::Nullable<arma::mat> R_SNP = R_NilValue ){

  NumericMatrix R_SNP_matrix;
  if(R_SNP.isNotNull()){ R_SNP_matrix = as<NumericMatrix>(R_SNP); };
  

  CharacterVector scenario = CharacterVector::create("A", "B", "C");
  IntegerVector snp_A, snp_B, snp_C, snp_AC, snp_BC;

  double p_A, p_B, p_C; p_A = 1.0; p_B = 1.0; p_C = 1.0;
  NumericVector alpha = NumericVector::create( p_A, p_B, p_C );
  NumericMatrix pai = rDirichlet( alpha, n_snps );
  NumericVector pai1, pai2, pai3;
  pai1 = pai(_, 0); pai2 = pai(_, 1); pai3 = pai(_, 2);

  CharacterVector yita(n_snps); 
  for(int i = 0; i < n_snps; ++i){
    NumericVector prob = NumericVector::create(pai1[i], pai2[i], pai3[i]);
    CharacterVector a = rsample( scenario, prob, 1 );
    yita[i] = a[0];
  };

  double beta = 0.3; NumericVector beta_new(iter_times); 
  NumericVector delta, kappax, kappay;
  delta = NumericVector(n_snps, 0.1); kappax = NumericVector(n_snps, 0.1); kappay = NumericVector(n_snps, 0.1);  
  double a_delta  = 0.0, b_delta  = 0.0, sigma_delta_sq  = 0.1; NumericVector sigma_delta_sq_new(iter_times); 
  double a_kappax = 0.0, b_kappax = 0.0, sigma_kappax_sq = 0.1; NumericVector sigma_kappax_sq_new(iter_times); 
  double a_kappay = 0.0, b_kappay = 0.0, sigma_kappay_sq = 0.1; NumericVector sigma_kappay_sq_new(iter_times); 
  double a_theta  = 1.0, b_theta  = 1.0, sigma_theta_sq  = 0.1; NumericVector sigma_theta_sq_new(iter_times);
  double miu_beta, sigma_beta_sq, sigma_beta;


  NumericVector Gamma = NumericVector(n_snps, 0.1); 
  NumericVector miu1_Gamma, miu2_Gamma, miu3_Gamma, sigma1_sq_Gamma, sigma2_sq_Gamma, sigma3_sq_Gamma;
                miu1_Gamma = NumericVector(n_snps); miu2_Gamma = NumericVector(n_snps); miu3_Gamma = NumericVector(n_snps);
                sigma1_sq_Gamma = NumericVector(n_snps); sigma2_sq_Gamma = NumericVector(n_snps); sigma3_sq_Gamma = NumericVector(n_snps);
  NumericVector gamma = NumericVector(n_snps, 0.1);
  NumericVector miu1_gamma, miu2_gamma, miu3_gamma, sigma1_sq_gamma, sigma2_sq_gamma, sigma3_sq_gamma;
                miu1_gamma = NumericVector(n_snps); miu2_gamma = NumericVector(n_snps); miu3_gamma = NumericVector(n_snps);
                sigma1_sq_gamma = NumericVector(n_snps); sigma2_sq_gamma = NumericVector(n_snps); sigma3_sq_gamma = NumericVector(n_snps);
  NumericVector miu1_kappay, miu2_kappay, miu3_kappay, sigma1_sq_kappay, sigma2_sq_kappay, sigma3_sq_kappay;
                miu1_kappay = NumericVector(n_snps); miu2_kappay = NumericVector(n_snps); miu3_kappay = NumericVector(n_snps);
                sigma1_sq_kappay = NumericVector(n_snps); sigma2_sq_kappay = NumericVector(n_snps); sigma3_sq_kappay = NumericVector(n_snps);
  
  NumericMatrix probs(n_snps, 3); 
  colnames(probs) = CharacterVector::create("p1", "p2", "p3");
  
 
  for( int iter = 0; iter < iter_times; ++iter ){
    for( int snp = 0; snp < n_snps; ++snp ){
      
      if(R_SNP.isNotNull() == false){
        

        sigma1_sq_Gamma[snp] = 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu1_Gamma[snp] = ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + beta*delta[snp]/sigma_theta_sq )*sigma1_sq_Gamma[snp];
        sigma2_sq_Gamma[snp] = 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu2_Gamma[snp] = ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + (beta*kappax[snp] + kappay[snp])/sigma_theta_sq )*sigma2_sq_Gamma[snp];
        sigma3_sq_Gamma[snp] = 1/( 1/data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu3_Gamma[snp] = ( data_mean_Gamma[snp]/data_s2_Gamma[snp] + (beta*(delta[snp] + kappax[snp]) + kappay[snp])/sigma_theta_sq )*sigma3_sq_Gamma[snp];
        NumericVector miu_Gamma = NumericVector::create( miu1_Gamma[snp], miu2_Gamma[snp], miu3_Gamma[snp] );
        NumericVector sigma_Gamma = NumericVector::create( sqrt(sigma1_sq_Gamma[snp]), sqrt(sigma2_sq_Gamma[snp]), sqrt(sigma3_sq_Gamma[snp]) );
        if( yita[snp] == "A" ){
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[0], sigma_Gamma[0]));
          Gamma[snp] = random_Gamma;
        } else if( yita[snp] == "B" ){
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[1], sigma_Gamma[1]));
          Gamma[snp] = random_Gamma;
        } else { 
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[2], sigma_Gamma[2]));
          Gamma[snp] = random_Gamma;
        };
        
        

        sigma1_sq_gamma[snp] = 1/( 1/data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/sigma_delta_sq );
        miu1_gamma[snp] = ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*Gamma[snp]/sigma_theta_sq )*sigma1_sq_gamma[snp];
        sigma2_sq_gamma[snp] = 1/( 1/data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/sigma_kappax_sq );
        miu2_gamma[snp] = ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma2_sq_gamma[snp];
        sigma3_sq_gamma[snp] = 1/( 1/data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/(sigma_delta_sq + sigma_kappax_sq) );
        miu3_gamma[snp] = ( data_mean_gamma[snp]/data_s2_gamma[snp] + beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma3_sq_gamma[snp];
        NumericVector miu_gamma = NumericVector::create(miu1_gamma[snp], miu2_gamma[snp], miu3_gamma[snp]);
        NumericVector sigma_gamma = NumericVector::create( sqrt(sigma1_sq_gamma[snp]), sqrt(sigma2_sq_gamma[snp]), sqrt(sigma3_sq_gamma[snp]) );
        if( yita[snp] == "A" ){
          double random_value = as<double>(rnorm_cpp(1, miu_gamma[0], sigma_gamma[0]));
          gamma[snp]  = random_value;
          delta[snp]  = random_value;
          kappax[snp] = 0.0;
        } else if( yita[snp] == "B" ){
          double random_value = as<double>(rnorm_cpp(1, miu_gamma[1], sigma_gamma[1]));
          gamma[snp]  = random_value;
          kappax[snp] = random_value;
          delta[snp]  = 0.0;
        } else{ 
          double random_delta  = as<double>(rnorm_cpp(1, miu_gamma[2] - kappax[snp], sigma_gamma[2]));
          delta[snp]  = random_delta;
          double random_kappax = as<double>(rnorm_cpp(1, miu_gamma[2] - delta[snp],  sigma_gamma[2]));
          kappax[snp] = random_kappax; 
          gamma[snp]  = delta[snp] + kappax[snp];
        };
        

        sigma1_sq_kappay[snp] = 0;
        miu1_kappay[snp] = 0;
        sigma2_sq_kappay[snp] = 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq );
        miu2_kappay[snp] = ( ( Gamma[snp] - beta*kappax[snp] )/sigma_theta_sq )*sigma2_sq_kappay[snp]; 
        sigma3_sq_kappay[snp] = 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq );
        miu3_kappay[snp] = ( ( Gamma[snp] - beta*(delta[snp] + kappax[snp]) )/sigma_theta_sq )*sigma3_sq_kappay[snp];
        NumericVector miu_kappay = NumericVector::create( miu1_kappay[snp], miu2_kappay[snp], miu3_kappay[snp] );
        NumericVector sigma_kappay = NumericVector::create( sqrt(sigma1_sq_kappay[snp]), sqrt(sigma2_sq_kappay[snp]), sqrt(sigma3_sq_kappay[snp]) );
        if( yita[snp] == "A" ){
          kappay[snp] = 0;
        } else if( yita[snp] == "B" ){
          double random_kappay = as<double>(rnorm_cpp(1, miu_kappay[1], sigma_kappay[1]) );
          kappay[snp] = random_kappay;
        } else{
          double random_kappay = as<double>(rnorm_cpp(1, miu_kappay[2], sigma_kappay[2]) );
        };
        

        NumericVector pai_elements(3);
        pai_elements = NumericVector::create( pai1[snp] * dnorm_cpp(gamma[snp], miu_gamma[0], sigma_gamma[0]),
                                              pai2[snp] * dnorm_cpp(gamma[snp], miu_gamma[1], sigma_gamma[1]),
                                              pai3[snp] * dnorm_cpp(gamma[snp], miu_gamma[2], sigma_gamma[2]));
        probs(snp, _) = pai_elements/(Rcpp::sum(pai_elements));
        yita[snp] = rsample( scenario, probs(snp, _), 1 )[0]; 
        
      } else{
        
        NumericVector Data_mean_gamma, Data_s2_gamma, Data_mean_Gamma, Data_s2_Gamma;
        Data_mean_gamma = arm_to_num(data_mean_gamma); Data_s2_gamma = arm_to_num(data_s2_gamma);
        Data_mean_Gamma = arm_to_num(data_mean_Gamma); Data_s2_Gamma = arm_to_num(data_s2_Gamma);
        NumericVector R = R_SNP_matrix(snp, _);
        IntegerVector cor_snps_ID = order_abs( R );

        
        sigma1_sq_Gamma[snp] = 1/( 1/Data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu1_Gamma[snp] = ( (Data_mean_Gamma/Data_s2_Gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_Gamma/sqrt(Data_s2_Gamma), cor_snps_ID ) ) )/sqrt(Data_s2_Gamma[snp]) +
          beta*delta[snp]/sigma_theta_sq )*sigma1_sq_Gamma[snp];
        sigma2_sq_Gamma[snp] = 1/( 1/Data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu2_Gamma[snp] = ( (Data_mean_Gamma/Data_s2_Gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_Gamma/sqrt(Data_s2_Gamma), cor_snps_ID ) ) )/sqrt(Data_s2_Gamma[snp]) +
          beta*(kappax[snp] + kappay[snp])/sigma_theta_sq )*sigma2_sq_Gamma[snp];
         
        sigma3_sq_Gamma[snp] = 1/( 1/Data_s2_Gamma[snp] + 1/sigma_theta_sq );
        miu3_Gamma[snp] = ( (Data_mean_Gamma/Data_s2_Gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_Gamma/sqrt(Data_s2_Gamma), cor_snps_ID ) ) )/sqrt(Data_s2_Gamma[snp]) +
                            (beta*(delta[snp] + kappax[snp]) + kappay[snp])/sigma_theta_sq )*sigma3_sq_Gamma[snp]; 
        
        NumericVector miu_Gamma = NumericVector::create( miu1_Gamma[snp], miu2_Gamma[snp], miu3_Gamma[snp] );
        NumericVector sigma_Gamma = NumericVector::create( sqrt(sigma1_sq_Gamma[snp]), sqrt(sigma2_sq_Gamma[snp]), sqrt(sigma3_sq_Gamma[snp]) );
        if( yita[snp] == "A" ){
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[0], sigma_Gamma[0]));
          Gamma[snp] = random_Gamma;
        } else if( yita[snp] == "B" ){  
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[1], sigma_Gamma[1]));
          Gamma[snp] = random_Gamma;
        } else {   
          double random_Gamma = as<double>(rnorm_cpp(1, miu_Gamma[2], sigma_Gamma[2]));
          Gamma[snp] = random_Gamma;
        };  
        
        
 
        sigma1_sq_gamma[snp] = 1/( 1/Data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/sigma_delta_sq );
        miu1_gamma[snp] = ( (Data_mean_gamma/Data_s2_gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_gamma/sqrt(Data_s2_gamma), cor_snps_ID ) ) )/sqrt(Data_s2_gamma[snp]) +
          beta*Gamma[snp]/sigma_theta_sq )*sigma1_sq_gamma[snp];
        sigma2_sq_gamma[snp] = 1/( 1/Data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/sigma_kappax_sq );
        miu2_gamma[snp] = ( (Data_mean_gamma/Data_s2_gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_gamma/sqrt(Data_s2_gamma), cor_snps_ID ) ) )/sqrt(Data_s2_gamma[snp]) +
          beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma2_sq_gamma[snp];
        sigma3_sq_gamma[snp] = 1/( 1/Data_s2_gamma[snp] + beta*beta/sigma_theta_sq + 1/(sigma_delta_sq + sigma_kappax_sq) );
        miu3_gamma[snp] = ( (Data_mean_gamma/Data_s2_gamma)[snp] - sum( rmindex( index_cpp( R*Data_mean_gamma/sqrt(Data_s2_gamma), cor_snps_ID ) ) )/sqrt(Data_s2_gamma[snp]) +
          beta*(Gamma[snp] - kappay[snp])/sigma_theta_sq )*sigma3_sq_gamma[snp];
        NumericVector miu_gamma = NumericVector::create(miu1_gamma[snp], miu2_gamma[snp], miu3_gamma[snp]);
        NumericVector sigma_gamma = NumericVector::create( sqrt(sigma1_sq_gamma[snp]), sqrt(sigma2_sq_gamma[snp]), sqrt(sigma3_sq_gamma[snp]) );
        if( yita[snp] == "A" ){
          double random_value = as<double>(rnorm_cpp(1, miu_gamma[0], sigma_gamma[0]));
          gamma[snp]  = random_value;
          delta[snp]  = random_value;
          kappax[snp] = 0.0;
        } else if( yita[snp] == "B" ){
          double random_value = as<double>(rnorm_cpp(1, miu_gamma[1], sigma_gamma[1]));
          gamma[snp]  = random_value;
          kappax[snp] = random_value;
          delta[snp]  = 0.0;
        } else{ 
          double random_delta  = as<double>(rnorm_cpp(1, miu_gamma[2] - kappax[snp], sigma_gamma[2]));
          delta[snp]  = random_delta;
          double random_kappax = as<double>(rnorm_cpp(1, miu_gamma[2] - delta[snp],  sigma_gamma[2]));
          kappax[snp] = random_kappax; 
          gamma[snp]  = delta[snp] + kappax[snp];
        };
        

        sigma1_sq_kappay[snp] = 0;
        miu1_kappay[snp] = 0;
        sigma2_sq_kappay[snp] = 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq );
        miu2_kappay[snp] = ( ( Gamma[snp] - beta*kappax[snp] )/sigma_theta_sq )*sigma2_sq_kappay[snp]; 
        sigma3_sq_kappay[snp] = 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq );
        miu3_kappay[snp] = ( ( Gamma[snp] - beta*(delta[snp] + kappax[snp]) )/sigma_theta_sq )*sigma3_sq_kappay[snp];
        NumericVector miu_kappay = NumericVector::create( miu1_kappay[snp], miu2_kappay[snp], miu3_kappay[snp] );
        NumericVector sigma_kappay = NumericVector::create( sqrt(sigma1_sq_kappay[snp]), sqrt(sigma2_sq_kappay[snp]), sqrt(sigma3_sq_kappay[snp]) );
        if( yita[snp] == "A" ){
          kappay[snp] = 0;
        } else if( yita[snp] == "B" ){
          double random_kappay = as<double>(rnorm_cpp(1, miu_kappay[1], sigma_kappay[1]) );
          kappay[snp] = random_kappay;
        } else{
          double random_kappay = as<double>(rnorm_cpp(1, miu_kappay[2], sigma_kappay[2]) );
        };
        
 
        NumericVector pai_elements(3);
        pai_elements = NumericVector::create( pai1[snp] * dnorm_cpp(gamma[snp], miu_gamma[0], sigma_gamma[0]),
                                              pai2[snp] * dnorm_cpp(gamma[snp], miu_gamma[1], sigma_gamma[1]),
                                              pai3[snp] * dnorm_cpp(gamma[snp], miu_gamma[2], sigma_gamma[2]));
        probs(snp, _) = pai_elements/(Rcpp::sum(pai_elements));
        yita[snp] = rsample( scenario, probs(snp, _), 1 )[0]; 
        
        
      }; 
      
    };
    
    
    CharacterVector A, B, C, AC, BC;
    A  = CharacterVector::create("A"); B = CharacterVector::create("B"); C = CharacterVector::create("C");
    AC = CharacterVector::create("A", "C"); BC = CharacterVector::create("B", "C");
    snp_A  = whichIn( yita, A );  snp_B  = whichIn( yita, B); snp_C = whichIn( yita, C );
    snp_AC = whichIn( yita, AC ); snp_BC = whichIn( yita, BC );
    
    double p_A_new, p_B_new, p_C_new;
    p_A_new = p_A + snp_A.size(); p_B_new = p_B + snp_B.size(); p_C_new = p_C + snp_C.size();
    pai = rDirichlet( NumericVector::create(p_A_new, p_B_new, p_C_new), n_snps );
    pai1 = pai(_, 0); pai2 = pai(_, 1); pai3 = pai(_, 2);
    
    

    sigma_beta_sq = 1/( ( vtv( index_cpp(delta, snp_A) ) + vtv( index_cpp(kappax, snp_B) ) + vtv( index_cpp(delta + kappax, snp_C) ) )/sigma_theta_sq );
    sigma_beta = sqrt( sigma_beta_sq );
    miu_beta = ( ( vtv( index_cpp(delta, snp_A), index_cpp(Gamma, snp_A) ) +
      vtv( index_cpp(kappax, snp_B), index_cpp( Gamma - kappay, snp_B) ) +
      vtv( index_cpp(delta + kappax, snp_C), index_cpp(Gamma - kappay, snp_C) ) ) / sigma_theta_sq ) * sigma_beta_sq;
    double random_beta = as<double>(rnorm_cpp(1, miu_beta, sigma_beta));
    beta_new[iter] = random_beta;
    beta = random_beta;
    

    double a_delta_new, b_delta_new;
    a_delta_new = a_delta + (snp_AC.size())/2;
    b_delta_new = b_delta + ( vtv( index_cpp(delta, snp_AC ) ) )/2;
    double random_sigma_delta_sq = as<double>( rinvgamma_cpp(1, a_delta_new, b_delta_new) );
    sigma_delta_sq_new[iter] = random_sigma_delta_sq; sigma_delta_sq = random_sigma_delta_sq;
    

    double a_kappax_new, b_kappax_new;
    a_kappax_new = a_kappax + (snp_BC.size())/2;
    b_kappax_new = b_kappax + vtv( index_cpp(kappax, snp_BC) )/2;
    double random_sigma_kappax_sq = as<double>( rinvgamma_cpp(1, a_kappax_new, b_kappax_new) );
    sigma_kappax_sq_new[iter] = random_sigma_kappax_sq; sigma_kappax_sq = random_sigma_kappax_sq;
    

    double a_kappay_new, b_kappay_new;
    a_kappay_new = a_kappay + (snp_BC.size())/2;
    b_kappay_new = b_kappay + vtv( index_cpp( kappay, snp_BC ) );
    double random_sigma_kappay_sq = as<double>( rinvgamma_cpp( 1, a_kappay_new, b_kappay_new ) );
    sigma_kappay_sq_new[iter] = random_sigma_kappay_sq; sigma_kappay_sq = random_sigma_kappay_sq;
    

    double a_theta_new, b_theta_new;
    a_theta_new = a_theta + n_snps/2;
    b_theta_new = b_theta + ( vtv( index_cpp(Gamma - beta*delta, snp_A) ) +
      vtv( index_cpp(Gamma - beta*kappax - kappay, snp_B) ) +
      vtv( index_cpp(Gamma - beta*(delta + kappax) - kappay, snp_C) ) )/2;
    double random_sigma_theta_sq = as<double>( rinvgamma_cpp( 1, a_theta_new, b_theta_new ) );
    sigma_theta_sq_new[iter] = random_sigma_theta_sq; sigma_theta_sq = random_sigma_theta_sq;
    
  }; 
  
  
  NumericVector beta_new_small = beta_new[Range(iter_times/2, iter_times - 1)];
  NumericVector output = res_output(beta_new_small);
  
  
  List result;
  result["delta"] = delta; result["kappax"] = kappax; result["kappay"] = kappay;
  result["gamma"] = gamma; result["Gamma"] = Gamma;
  result["probs"] = probs;
  result["senario"] = yita;
  result["beta_new"] = beta_new;
  result["output"] = output;

  return result;
}


















// [[Rcpp::export]]

NumericVector rnorm_cpp(int n, double mean, double sd) {
  NumericVector result = rnorm(n, mean, sd);
  return result;
}  

// [[Rcpp::export]]

double dnorm_cpp(double x, double mean, double sd) {
  double exponent = -0.5 * pow((x - mean) / sd, 2);
  double coefficient = 1.0 / (sd * sqrt(2 * M_PI));
  double pdf = coefficient * exp(exponent);
  return pdf;
} // 核验正确 

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

NumericMatrix mvrnorm_cpp(int n, NumericVector mean, NumericMatrix cov_matrix) {
  int p = mean.size(); 
  
  Eigen::VectorXd eigen_mean = as<Eigen::VectorXd>(mean);
  Eigen::MatrixXd eigen_cov_matrix = as<Eigen::MatrixXd>(cov_matrix);
  Eigen::MatrixXd random_matrix = Eigen::MatrixXd::Zero(n, p);
  
  Eigen::LLT<Eigen::MatrixXd> lltOfCovMatrix(eigen_cov_matrix);
   
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd random_vector(p);
    for (int j = 0; j < p; ++j) {
      random_vector(j) = R::rnorm(0.0, 1.0); 
    }   
    Eigen::VectorXd correlated_vector = lltOfCovMatrix.matrixL() * random_vector + eigen_mean;
    random_matrix.row(i) = correlated_vector;
  }   
  return Rcpp::wrap(random_matrix);
}   

// [[Rcpp::export]]
double runi_cpp(double min_val, double max_val) {
  std::uniform_real_distribution<double> distribution(min_val, max_val);
  std::mt19937 generator(std::random_device{}());
  double result = distribution(generator);
  return result;
}  

// [[Rcpp::export]]
NumericVector rinvgamma_cpp(int n, double shape, double scale) {
  NumericVector result = Rcpp::rgamma(n, shape, 1.0 / scale);
  return 1.0 / result;
} 

// [[Rcpp::export]]
CharacterVector rsample(CharacterVector characters, NumericVector probs, int sample_size) {
  int n = characters.size();
  CharacterVector output(sample_size);
  
  NumericVector cum_probs(n);
  double sum_probs = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_probs += probs[i];
    cum_probs[i] = sum_probs;
  }  
  
  for (int i = 0; i < sample_size; ++i) {
    double rand_val = R::runif(0, 1);
    for (int j = 0; j < n; ++j) {
      if (rand_val < cum_probs[j]) {
        output[i] = characters[j];
        break;
      }
    } 
  }  
  
  return output;
} 

// [[Rcpp::export]] 
NumericMatrix rDirichlet( NumericVector alpha, int times ){
  int num_dimensions = alpha.size();
  NumericMatrix rmatrix(times, num_dimensions); 
  
  for( int i = 0; i < times; ++i ){
    for( int j = 0; j < num_dimensions; ++j ){
      rmatrix(i, j) = R::rgamma( alpha(j), 1.0 );
    }  
    rmatrix.row(i) = rmatrix.row(i)/sum(rmatrix.row(i));
  }  
  return rmatrix; 
} 

// [[Rcpp::export]]
List lm_cpp(NumericVector y, NumericMatrix X) {
  int n = y.size(); int times = X.ncol(); 
  NumericVector coef(times); NumericVector std(times);
  for(int i = 0; i < times; ++i){
    NumericVector X_i = X(_, i);
    Eigen::MatrixXd matrix(n, 2);
    for(int j = 0; j < n; ++j){
      matrix(j, 0) = 1.0;
      matrix(j, 1) = X_i[j];
    } 
    
   
    Eigen::VectorXd Y(n);
    for(int j = 0; j < n; ++j){
      Y(j) = y[j];
    }  
    
    
    Eigen::VectorXd coef_vec = matrix.colPivHouseholderQr().solve(Y);
    Eigen::VectorXd res = Y - matrix*coef_vec;
   
    double rss = res.squaredNorm();
    
    double std_error = std::sqrt(rss/(n - 2));
   
    Eigen::VectorXd se_coef = std_error * (matrix.transpose() * matrix).inverse().diagonal().array().sqrt();
    
    coef[i] = coef_vec(1); std[i] = se_coef(1);
  }  
  
  List result;
  result["coef"] = coef; result["std"] = std;
  return result; 
} 

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix cor_cpp(NumericMatrix input_matrix) {
  const Eigen::Map<Eigen::MatrixXd> eigen_matrix(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(input_matrix));
  Eigen::MatrixXd centered_matrix = eigen_matrix.rowwise() - eigen_matrix.colwise().mean();
  Eigen::MatrixXd cov_matrix = (centered_matrix.adjoint() * centered_matrix) / (eigen_matrix.rows() - 1);
  
  Eigen::MatrixXd correlation_matrix = cov_matrix;
  for (int i = 0; i < cov_matrix.rows(); ++i) {
    double sd_i = std::sqrt(cov_matrix(i, i));
    for (int j = 0; j < cov_matrix.cols(); ++j) {
      double sd_j = std::sqrt(cov_matrix(j, j));
      correlation_matrix(i, j) = cov_matrix(i, j) / (sd_i * sd_j);
    }
  }  
  
  NumericMatrix result = Rcpp::wrap(correlation_matrix);
  return result;
} 

// [[Rcpp::export]]
IntegerVector whichIn(CharacterVector x, CharacterVector values) {
  IntegerVector indices;
  
  for( int i = 0; i < values.size(); ++i ){
    for( int j = 0; j < x.size(); ++j ){
      
      if( x[j] == values[i] ) indices.push_back(j + 1);  
      
    }
  } 
  
  return indices;
} 



// [[Rcpp::export]]
double vtv(NumericVector x, NumericVector y){
  int len_x = x.size();
  int len_y = y.size();
  double result_sum = 0.0;
  if(len_y == 0){
    for(int i = 0; i < len_x; ++i){
      result_sum += x[i]*x[i];
    }
  } else{
    for(int i = 0; i < len_x; ++i){
      result_sum += x[i]*y[i];
    }
  }
  return result_sum;
} 


// [[Rcpp::export]]
NumericVector qt_cpp(NumericVector x, NumericVector probs){
  int n_probs = probs.size();
  NumericVector result(n_probs);
  
  NumericVector clean_x = na_omit(x);

  Function myqt("quantile"); 
  
  for(int i = 0; i < n_probs; ++i){
    result[i] = as<double>(myqt(clean_x, probs[i]));
  }
  
  return result;
} 

// [[Rcpp::export]]
NumericVector res_output(NumericVector x) {
  double q0025 = 0.025, median = 0.5, q0975 = 0.975;
  Function mysd("sd");
  
  double std_dev = as<double>( mysd(x) );
  
  NumericVector result = qt_cpp(x, NumericVector::create(median, q0025, q0975) );
  result.push_back(std_dev);
  
  return result;
}

// [[Rcpp::export]]
NumericVector index_cpp(NumericVector x , IntegerVector index){
  int len = index.size();
  NumericVector result(len);
  
  if( len == 0 ){
    return result = NumericVector();
  } else{
    for( int i = 0; i < len; ++ i){
      result[i] = x[index[i] - 1];
    }
    return result;
  }
  
} 

// [[Rcpp::export]]
NumericVector minor_genotype(NumericVector x){
  int len = x.size();
  double minor = runi_cpp(0.05, 0.5);
  double Q1 = minor*minor;
  double Q2 = minor*(2 - minor);
  NumericVector qt = qt_cpp( x, NumericVector::create(Q1, Q2) );
  NumericVector result(len);
  for(int i = 0; i < len; ++i){
    if(x[i] <= qt[0]){
      result[i] = 0;
    } else if(x[i] <= qt[1]){
      result[i] = 1;
    } else{
      result[i] = 2;
    }
  }
  return result;
} 

// [[Rcpp::export]]
List genotype_generation(List size_sample){
  int nrow_x = size_sample["nrow_x"]; int nrow_y = size_sample["nrow_y"]; Nullable<int> nrow_r = size_sample["nrow_r"]; 
  int n_snps = size_sample["n_snps"]; double auto_r = size_sample["auto_r"];
  
  NumericVector miu_x = NumericVector(n_snps);
  NumericVector miu_y = NumericVector(n_snps);
  NumericMatrix sigma2(n_snps, n_snps);
  for(int i = 0; i < n_snps; ++i){
    for(int j = 0; j < n_snps; ++j){
      if(i == j){
        sigma2(i, j) = 1.0;
      } else{
        sigma2(i, j) = std::pow(auto_r, std::abs(i - j));
      }
    }
  }
  
  NumericMatrix data_xx; NumericMatrix G_x(nrow_x, n_snps);
  data_xx = mvrnorm_cpp(nrow_x, miu_x, sigma2);
  for(int i = 0; i < n_snps; ++i){
    NumericVector x = data_xx(_, i);
    NumericVector y = minor_genotype(x);
  
    G_x(_, i) = y;
  }
  
  NumericMatrix data_yy; NumericMatrix G_y(nrow_y, n_snps);
  data_yy = mvrnorm_cpp(nrow_y, miu_y, sigma2);
  for(int i = 0; i < n_snps; ++i){
    NumericVector x = data_yy(_, i);
    NumericVector y = minor_genotype(x);
   
    G_y(_, i) = y;
  }
  
  List output;
  
  if( nrow_r.isNotNull() ){
    int nrow = as<int>(nrow_r); NumericVector miu_r = NumericVector(n_snps);
    NumericMatrix data_rr; NumericMatrix G_r(nrow, n_snps);
    data_rr = mvrnorm_cpp(nrow, miu_r, sigma2);
    for(int i = 0; i < n_snps; ++i){
      NumericVector x = data_rr(_, i);
      NumericVector y = minor_genotype(x);
      
      G_r(_, i) = y;
    };
    NumericMatrix R_r = cor_cpp(G_r);
    output["G_x"] = G_x; 
    output["G_y"] = G_y; 
    output["R_r"] = R_r;
    
  } else{
    output["G_x"] = G_x; 
    output["G_y"] = G_y;
  }
  
  return output;
  
}  

// [[Rcpp::export]]
IntegerVector order_abs(NumericVector x, bool decreasing) {
  int n = x.size();
  IntegerVector indices(n);
  
  
  for (int i = 0; i < n; ++i) {
    indices[i] = i + 1;
  }
  
  
  if( decreasing == true ){
    std::sort(indices.begin(), indices.end(),
              [&x](int i, int j) { return std::abs(x[i - 1]) >= std::abs(x[j - 1]); } );
  } else{ 
    std::sort(indices.begin(), indices.end(),
              [&x](int i, int j) { return std::abs(x[i - 1]) <= std::abs(x[j - 1]); } );
  } 
  
  
 
  IntegerVector top_indices = indices[Rcpp::Range(0, 49)];
  
  return top_indices;
}

// [[Rcpp::export]]
NumericVector rmindex(NumericVector x, IntegerVector indices) {
  int numToRemove = indices.size();
  

  std::sort(indices.begin(), indices.end(), std::greater<int>());
  
  for (int i = 0; i < numToRemove; ++i) {
    int indexToRemove;
    if(indices[i] == 0){
      indexToRemove = 0;
    } else{
      indexToRemove = indices[i] - 1;
    }; 
    if (indexToRemove >= 0 && indexToRemove < x.size()) {
      x.erase(x.begin() + indexToRemove); // Erase element at the specified index
    }
  } 
  
  return x;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector arm_to_num(arma::vec x) {
  int n = x.n_elem;
  NumericVector result(n);
  
  for (int i = 0; i < n; ++i) {
    result[i] = x[i];
  }
  
  return result;
}
