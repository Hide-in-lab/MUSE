// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MUSE
List MUSE(arma::vec& data_mean_gamma, arma::vec& data_s2_gamma, arma::vec& data_mean_Gamma, arma::vec& data_s2_Gamma, int n_snps, int iter_times, Rcpp::Nullable<arma::mat> R_SNP);
RcppExport SEXP _MUSE_MUSE(SEXP data_mean_gammaSEXP, SEXP data_s2_gammaSEXP, SEXP data_mean_GammaSEXP, SEXP data_s2_GammaSEXP, SEXP n_snpsSEXP, SEXP iter_timesSEXP, SEXP R_SNPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type data_mean_gamma(data_mean_gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type data_s2_gamma(data_s2_gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type data_mean_Gamma(data_mean_GammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type data_s2_Gamma(data_s2_GammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_snps(n_snpsSEXP);
    Rcpp::traits::input_parameter< int >::type iter_times(iter_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type R_SNP(R_SNPSEXP);
    rcpp_result_gen = Rcpp::wrap(MUSE(data_mean_gamma, data_s2_gamma, data_mean_Gamma, data_s2_Gamma, n_snps, iter_times, R_SNP));
    return rcpp_result_gen;
END_RCPP
}
// rnorm_cpp
NumericVector rnorm_cpp(int n, double mean, double sd);
RcppExport SEXP _MUSE_rnorm_cpp(SEXP nSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm_cpp(n, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// dnorm_cpp
double dnorm_cpp(double x, double mean, double sd);
RcppExport SEXP _MUSE_dnorm_cpp(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(dnorm_cpp(x, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// mvrnorm_cpp
NumericMatrix mvrnorm_cpp(int n, NumericVector mean, NumericMatrix cov_matrix);
RcppExport SEXP _MUSE_mvrnorm_cpp(SEXP nSEXP, SEXP meanSEXP, SEXP cov_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cov_matrix(cov_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnorm_cpp(n, mean, cov_matrix));
    return rcpp_result_gen;
END_RCPP
}
// runi_cpp
double runi_cpp(double min_val, double max_val);
RcppExport SEXP _MUSE_runi_cpp(SEXP min_valSEXP, SEXP max_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type min_val(min_valSEXP);
    Rcpp::traits::input_parameter< double >::type max_val(max_valSEXP);
    rcpp_result_gen = Rcpp::wrap(runi_cpp(min_val, max_val));
    return rcpp_result_gen;
END_RCPP
}
// rinvgamma_cpp
NumericVector rinvgamma_cpp(int n, double shape, double scale);
RcppExport SEXP _MUSE_rinvgamma_cpp(SEXP nSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvgamma_cpp(n, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// rsample
CharacterVector rsample(CharacterVector characters, NumericVector probs, int sample_size);
RcppExport SEXP _MUSE_rsample(SEXP charactersSEXP, SEXP probsSEXP, SEXP sample_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type characters(charactersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(rsample(characters, probs, sample_size));
    return rcpp_result_gen;
END_RCPP
}
// rDirichlet
NumericMatrix rDirichlet(NumericVector alpha, int times);
RcppExport SEXP _MUSE_rDirichlet(SEXP alphaSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(rDirichlet(alpha, times));
    return rcpp_result_gen;
END_RCPP
}
// lm_cpp
List lm_cpp(NumericVector y, NumericMatrix X);
RcppExport SEXP _MUSE_lm_cpp(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(lm_cpp(y, X));
    return rcpp_result_gen;
END_RCPP
}
// cor_cpp
NumericMatrix cor_cpp(NumericMatrix input_matrix);
RcppExport SEXP _MUSE_cor_cpp(SEXP input_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type input_matrix(input_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_cpp(input_matrix));
    return rcpp_result_gen;
END_RCPP
}
// whichIn
IntegerVector whichIn(CharacterVector x, CharacterVector values);
RcppExport SEXP _MUSE_whichIn(SEXP xSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(whichIn(x, values));
    return rcpp_result_gen;
END_RCPP
}
// vtv
double vtv(NumericVector x, NumericVector y);
RcppExport SEXP _MUSE_vtv(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vtv(x, y));
    return rcpp_result_gen;
END_RCPP
}
// qt_cpp
NumericVector qt_cpp(NumericVector x, NumericVector probs);
RcppExport SEXP _MUSE_qt_cpp(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(qt_cpp(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// res_output
NumericVector res_output(NumericVector x);
RcppExport SEXP _MUSE_res_output(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(res_output(x));
    return rcpp_result_gen;
END_RCPP
}
// index_cpp
NumericVector index_cpp(NumericVector x, IntegerVector index);
RcppExport SEXP _MUSE_index_cpp(SEXP xSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(index_cpp(x, index));
    return rcpp_result_gen;
END_RCPP
}
// minor_genotype
NumericVector minor_genotype(NumericVector x);
RcppExport SEXP _MUSE_minor_genotype(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(minor_genotype(x));
    return rcpp_result_gen;
END_RCPP
}
// genotype_generation
List genotype_generation(List size_sample);
RcppExport SEXP _MUSE_genotype_generation(SEXP size_sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type size_sample(size_sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(genotype_generation(size_sample));
    return rcpp_result_gen;
END_RCPP
}
// order_abs
IntegerVector order_abs(NumericVector x, bool decreasing);
RcppExport SEXP _MUSE_order_abs(SEXP xSEXP, SEXP decreasingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type decreasing(decreasingSEXP);
    rcpp_result_gen = Rcpp::wrap(order_abs(x, decreasing));
    return rcpp_result_gen;
END_RCPP
}
// rmindex
NumericVector rmindex(NumericVector x, IntegerVector indices);
RcppExport SEXP _MUSE_rmindex(SEXP xSEXP, SEXP indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(rmindex(x, indices));
    return rcpp_result_gen;
END_RCPP
}
// arm_to_num
NumericVector arm_to_num(arma::vec x);
RcppExport SEXP _MUSE_arm_to_num(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(arm_to_num(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MUSE_MUSE", (DL_FUNC) &_MUSE_MUSE, 7},
    {"_MUSE_rnorm_cpp", (DL_FUNC) &_MUSE_rnorm_cpp, 3},
    {"_MUSE_dnorm_cpp", (DL_FUNC) &_MUSE_dnorm_cpp, 3},
    {"_MUSE_mvrnorm_cpp", (DL_FUNC) &_MUSE_mvrnorm_cpp, 3},
    {"_MUSE_runi_cpp", (DL_FUNC) &_MUSE_runi_cpp, 2},
    {"_MUSE_rinvgamma_cpp", (DL_FUNC) &_MUSE_rinvgamma_cpp, 3},
    {"_MUSE_rsample", (DL_FUNC) &_MUSE_rsample, 3},
    {"_MUSE_rDirichlet", (DL_FUNC) &_MUSE_rDirichlet, 2},
    {"_MUSE_lm_cpp", (DL_FUNC) &_MUSE_lm_cpp, 2},
    {"_MUSE_cor_cpp", (DL_FUNC) &_MUSE_cor_cpp, 1},
    {"_MUSE_whichIn", (DL_FUNC) &_MUSE_whichIn, 2},
    {"_MUSE_vtv", (DL_FUNC) &_MUSE_vtv, 2},
    {"_MUSE_qt_cpp", (DL_FUNC) &_MUSE_qt_cpp, 2},
    {"_MUSE_res_output", (DL_FUNC) &_MUSE_res_output, 1},
    {"_MUSE_index_cpp", (DL_FUNC) &_MUSE_index_cpp, 2},
    {"_MUSE_minor_genotype", (DL_FUNC) &_MUSE_minor_genotype, 1},
    {"_MUSE_genotype_generation", (DL_FUNC) &_MUSE_genotype_generation, 1},
    {"_MUSE_order_abs", (DL_FUNC) &_MUSE_order_abs, 2},
    {"_MUSE_rmindex", (DL_FUNC) &_MUSE_rmindex, 2},
    {"_MUSE_arm_to_num", (DL_FUNC) &_MUSE_arm_to_num, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MUSE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
