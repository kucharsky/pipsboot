// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// UPbrewer_rcpp
IntegerVector UPbrewer_rcpp(const NumericVector& pik);
RcppExport SEXP _pipsboot_UPbrewer_rcpp(SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(UPbrewer_rcpp(pik));
    return rcpp_result_gen;
END_RCPP
}
// UPpoisson_rcpp
IntegerVector UPpoisson_rcpp(const NumericVector& pik);
RcppExport SEXP _pipsboot_UPpoisson_rcpp(SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(UPpoisson_rcpp(pik));
    return rcpp_result_gen;
END_RCPP
}
// boot_05_pips_rcpp
IntegerMatrix boot_05_pips_rcpp(const NumericVector& x_sample, const double& x_pop_total, const int& n_boots, const Rcpp::Function& scheme);
RcppExport SEXP _pipsboot_boot_05_pips_rcpp(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP, SEXP schemeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< const double& >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_boots(n_bootsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type scheme(schemeSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_05_pips_rcpp(x_sample, x_pop_total, n_boots, scheme));
    return rcpp_result_gen;
END_RCPP
}
// boot_05_pips_rcpp_brewer
IntegerMatrix boot_05_pips_rcpp_brewer(const NumericVector& x_sample, const double& x_pop_total, const int& n_boots);
RcppExport SEXP _pipsboot_boot_05_pips_rcpp_brewer(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< const double& >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_05_pips_rcpp_brewer(x_sample, x_pop_total, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_AT2011_alg3_rcpp
IntegerMatrix boot_AT2011_alg3_rcpp(int N, int n, int n_boots);
RcppExport SEXP _pipsboot_boot_AT2011_alg3_rcpp(SEXP NSEXP, SEXP nSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_AT2011_alg3_rcpp(N, n, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_AT2011_alg45_rcpp
IntegerMatrix boot_AT2011_alg45_rcpp(const NumericVector& pi_sample, const int& n_boots, const Rcpp::Function& scheme);
RcppExport SEXP _pipsboot_boot_AT2011_alg45_rcpp(SEXP pi_sampleSEXP, SEXP n_bootsSEXP, SEXP schemeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_sample(pi_sampleSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_boots(n_bootsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type scheme(schemeSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_AT2011_alg45_rcpp(pi_sample, n_boots, scheme));
    return rcpp_result_gen;
END_RCPP
}
// boot_AT2014_rcpp
IntegerMatrix boot_AT2014_rcpp(const NumericVector& pi_sample, int n_boots);
RcppExport SEXP _pipsboot_boot_AT2014_rcpp(SEXP pi_sampleSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_sample(pi_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_AT2014_rcpp(pi_sample, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_Q_pips_rcpp
IntegerMatrix boot_Q_pips_rcpp(const NumericVector& x_sample, const double& x_pop_total, int n_boots);
RcppExport SEXP _pipsboot_boot_Q_pips_rcpp(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< const double& >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_Q_pips_rcpp(x_sample, x_pop_total, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_Qg_pips_rcpp
IntegerMatrix boot_Qg_pips_rcpp(const NumericVector& x_sample, const double& x_pop_total, const NumericVector& g_weights, int n_boots);
RcppExport SEXP _pipsboot_boot_Qg_pips_rcpp(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP g_weightsSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< const double& >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type g_weights(g_weightsSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_Qg_pips_rcpp(x_sample, x_pop_total, g_weights, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_SP_pips_rcpp
IntegerMatrix boot_SP_pips_rcpp(const NumericVector& x_sample, const NumericVector& pi_sample, const int& N, const int& n_boots, const Rcpp::Function& scheme);
RcppExport SEXP _pipsboot_boot_SP_pips_rcpp(SEXP x_sampleSEXP, SEXP pi_sampleSEXP, SEXP NSEXP, SEXP n_bootsSEXP, SEXP schemeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_sample(pi_sampleSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_boots(n_bootsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type scheme(schemeSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_SP_pips_rcpp(x_sample, pi_sample, N, n_boots, scheme));
    return rcpp_result_gen;
END_RCPP
}
// boot_Z_pips_rcpp
IntegerMatrix boot_Z_pips_rcpp(const NumericVector& n_replications, const NumericVector& pi_sample, int n_boots);
RcppExport SEXP _pipsboot_boot_Z_pips_rcpp(SEXP n_replicationsSEXP, SEXP pi_sampleSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type n_replications(n_replicationsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_sample(pi_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_Z_pips_rcpp(n_replications, pi_sample, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_xbal1_pips_rcpp
IntegerMatrix boot_xbal1_pips_rcpp(const NumericVector& x_sample, double x_pop_total, int n_boots, Function& scheme);
RcppExport SEXP _pipsboot_boot_xbal1_pips_rcpp(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP, SEXP schemeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< double >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    Rcpp::traits::input_parameter< Function& >::type scheme(schemeSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_xbal1_pips_rcpp(x_sample, x_pop_total, n_boots, scheme));
    return rcpp_result_gen;
END_RCPP
}
// boot_xbal1_pips_rcpp_brewer
IntegerMatrix boot_xbal1_pips_rcpp_brewer(const NumericVector& x_sample, double x_pop_total, int n_boots);
RcppExport SEXP _pipsboot_boot_xbal1_pips_rcpp_brewer(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< double >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_xbal1_pips_rcpp_brewer(x_sample, x_pop_total, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// boot_xbal2_pips_rcpp
IntegerMatrix boot_xbal2_pips_rcpp(const NumericVector& x_sample, double x_pop_total, int n_boots, Function& scheme);
RcppExport SEXP _pipsboot_boot_xbal2_pips_rcpp(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP, SEXP schemeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< double >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    Rcpp::traits::input_parameter< Function& >::type scheme(schemeSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_xbal2_pips_rcpp(x_sample, x_pop_total, n_boots, scheme));
    return rcpp_result_gen;
END_RCPP
}
// boot_xbal2_pips_rcpp_brewer
IntegerMatrix boot_xbal2_pips_rcpp_brewer(const NumericVector& x_sample, double x_pop_total, int n_boots);
RcppExport SEXP _pipsboot_boot_xbal2_pips_rcpp_brewer(SEXP x_sampleSEXP, SEXP x_pop_totalSEXP, SEXP n_bootsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    Rcpp::traits::input_parameter< double >::type x_pop_total(x_pop_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n_boots(n_bootsSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_xbal2_pips_rcpp_brewer(x_sample, x_pop_total, n_boots));
    return rcpp_result_gen;
END_RCPP
}
// halfdoubled_rcpp
IntegerVector halfdoubled_rcpp(const int& n);
RcppExport SEXP _pipsboot_halfdoubled_rcpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(halfdoubled_rcpp(n));
    return rcpp_result_gen;
END_RCPP
}
// inclusion_probabilities_rcpp
NumericVector inclusion_probabilities_rcpp(const NumericVector& values, const int& size);
RcppExport SEXP _pipsboot_inclusion_probabilities_rcpp(SEXP valuesSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< const int& >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(inclusion_probabilities_rcpp(values, size));
    return rcpp_result_gen;
END_RCPP
}
// inverse_hypergeometric_rcpp_scalar
int inverse_hypergeometric_rcpp_scalar(const int& N, const int& n);
RcppExport SEXP _pipsboot_inverse_hypergeometric_rcpp_scalar(SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(inverse_hypergeometric_rcpp_scalar(N, n));
    return rcpp_result_gen;
END_RCPP
}
// sampledata_bootrep_rcpp
NumericMatrix sampledata_bootrep_rcpp(const IntegerMatrix& bootstrap_samples, const NumericVector& x_sample);
RcppExport SEXP _pipsboot_sampledata_bootrep_rcpp(SEXP bootstrap_samplesSEXP, SEXP x_sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type bootstrap_samples(bootstrap_samplesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_sample(x_sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(sampledata_bootrep_rcpp(bootstrap_samples, x_sample));
    return rcpp_result_gen;
END_RCPP
}
// sampling_oneone_rcpp
IntegerVector sampling_oneone_rcpp(const int& n);
RcppExport SEXP _pipsboot_sampling_oneone_rcpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_oneone_rcpp(n));
    return rcpp_result_gen;
END_RCPP
}
// sampling_with_overreplacement_rcpp_scalar
IntegerVector sampling_with_overreplacement_rcpp_scalar(const int& N, const int& n);
RcppExport SEXP _pipsboot_sampling_with_overreplacement_rcpp_scalar(SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_with_overreplacement_rcpp_scalar(N, n));
    return rcpp_result_gen;
END_RCPP
}
// srswor_rcpp_01
IntegerVector srswor_rcpp_01(const int& N, const int& n);
RcppExport SEXP _pipsboot_srswor_rcpp_01(SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(srswor_rcpp_01(N, n));
    return rcpp_result_gen;
END_RCPP
}
// srswor_rcpp_pos
IntegerVector srswor_rcpp_pos(const int& N, const int& n);
RcppExport SEXP _pipsboot_srswor_rcpp_pos(SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(srswor_rcpp_pos(N, n));
    return rcpp_result_gen;
END_RCPP
}
// srswr_rcpp_zero
IntegerVector srswr_rcpp_zero(const int& N, const int& n);
RcppExport SEXP _pipsboot_srswr_rcpp_zero(SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(srswr_rcpp_zero(N, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pipsboot_UPbrewer_rcpp", (DL_FUNC) &_pipsboot_UPbrewer_rcpp, 1},
    {"_pipsboot_UPpoisson_rcpp", (DL_FUNC) &_pipsboot_UPpoisson_rcpp, 1},
    {"_pipsboot_boot_05_pips_rcpp", (DL_FUNC) &_pipsboot_boot_05_pips_rcpp, 4},
    {"_pipsboot_boot_05_pips_rcpp_brewer", (DL_FUNC) &_pipsboot_boot_05_pips_rcpp_brewer, 3},
    {"_pipsboot_boot_AT2011_alg3_rcpp", (DL_FUNC) &_pipsboot_boot_AT2011_alg3_rcpp, 3},
    {"_pipsboot_boot_AT2011_alg45_rcpp", (DL_FUNC) &_pipsboot_boot_AT2011_alg45_rcpp, 3},
    {"_pipsboot_boot_AT2014_rcpp", (DL_FUNC) &_pipsboot_boot_AT2014_rcpp, 2},
    {"_pipsboot_boot_Q_pips_rcpp", (DL_FUNC) &_pipsboot_boot_Q_pips_rcpp, 3},
    {"_pipsboot_boot_Qg_pips_rcpp", (DL_FUNC) &_pipsboot_boot_Qg_pips_rcpp, 4},
    {"_pipsboot_boot_SP_pips_rcpp", (DL_FUNC) &_pipsboot_boot_SP_pips_rcpp, 5},
    {"_pipsboot_boot_Z_pips_rcpp", (DL_FUNC) &_pipsboot_boot_Z_pips_rcpp, 3},
    {"_pipsboot_boot_xbal1_pips_rcpp", (DL_FUNC) &_pipsboot_boot_xbal1_pips_rcpp, 4},
    {"_pipsboot_boot_xbal1_pips_rcpp_brewer", (DL_FUNC) &_pipsboot_boot_xbal1_pips_rcpp_brewer, 3},
    {"_pipsboot_boot_xbal2_pips_rcpp", (DL_FUNC) &_pipsboot_boot_xbal2_pips_rcpp, 4},
    {"_pipsboot_boot_xbal2_pips_rcpp_brewer", (DL_FUNC) &_pipsboot_boot_xbal2_pips_rcpp_brewer, 3},
    {"_pipsboot_halfdoubled_rcpp", (DL_FUNC) &_pipsboot_halfdoubled_rcpp, 1},
    {"_pipsboot_inclusion_probabilities_rcpp", (DL_FUNC) &_pipsboot_inclusion_probabilities_rcpp, 2},
    {"_pipsboot_inverse_hypergeometric_rcpp_scalar", (DL_FUNC) &_pipsboot_inverse_hypergeometric_rcpp_scalar, 2},
    {"_pipsboot_sampledata_bootrep_rcpp", (DL_FUNC) &_pipsboot_sampledata_bootrep_rcpp, 2},
    {"_pipsboot_sampling_oneone_rcpp", (DL_FUNC) &_pipsboot_sampling_oneone_rcpp, 1},
    {"_pipsboot_sampling_with_overreplacement_rcpp_scalar", (DL_FUNC) &_pipsboot_sampling_with_overreplacement_rcpp_scalar, 2},
    {"_pipsboot_srswor_rcpp_01", (DL_FUNC) &_pipsboot_srswor_rcpp_01, 2},
    {"_pipsboot_srswor_rcpp_pos", (DL_FUNC) &_pipsboot_srswor_rcpp_pos, 2},
    {"_pipsboot_srswr_rcpp_zero", (DL_FUNC) &_pipsboot_srswr_rcpp_zero, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pipsboot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
