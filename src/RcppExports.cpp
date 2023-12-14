// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fastResidop
SEXP fastResidop(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _RUVIIIPRPS_fastResidop(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastResidop(A, B));
    return rcpp_result_gen;
END_RCPP
}
// fastResidop2
SEXP fastResidop2(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _RUVIIIPRPS_fastResidop2(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastResidop2(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matrixMult
SEXP matrixMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _RUVIIIPRPS_matrixMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matSubtraction
SEXP matSubtraction(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _RUVIIIPRPS_matSubtraction(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matSubtraction(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matTranspose
SEXP matTranspose(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _RUVIIIPRPS_matTranspose(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matTranspose(A));
    return rcpp_result_gen;
END_RCPP
}
// matInverse
SEXP matInverse(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _RUVIIIPRPS_matInverse(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matInverse(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RUVIIIPRPS_fastResidop", (DL_FUNC) &_RUVIIIPRPS_fastResidop, 2},
    {"_RUVIIIPRPS_fastResidop2", (DL_FUNC) &_RUVIIIPRPS_fastResidop2, 2},
    {"_RUVIIIPRPS_matrixMult", (DL_FUNC) &_RUVIIIPRPS_matrixMult, 2},
    {"_RUVIIIPRPS_matSubtraction", (DL_FUNC) &_RUVIIIPRPS_matSubtraction, 2},
    {"_RUVIIIPRPS_matTranspose", (DL_FUNC) &_RUVIIIPRPS_matTranspose, 1},
    {"_RUVIIIPRPS_matInverse", (DL_FUNC) &_RUVIIIPRPS_matInverse, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RUVIIIPRPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}