#include <Rcpp.h>
#include "../inst/include/FLasherEMSRR.h"

// Diagnostic only: verify whether the R-backed SRR callback preserves
// the derivative of recruitment with respect to SSB.
// [[Rcpp::export]]
Rcpp::NumericVector customSRR_AD_diagnostic(
    double ssb,
    std::string model_name,
    Rcpp::NumericVector params
) {
    using AD = CppAD::AD<double>;

    std::vector<AD> ax(1);
    ax[0] = ssb;
    CppAD::Independent(ax);

    std::vector<double> p(params.begin(), params.end());
    std::vector<AD> ay(1);
    ay[0] = customSRR<AD>(ax[0], p, model_name);

    CppAD::ADFun<double> f;
    f.Dependent(ax, ay);

    std::vector<double> x0(1, ssb);
    std::vector<double> value = f.Forward(0, x0);
    std::vector<double> jacobian = f.Jacobian(x0);

    return Rcpp::NumericVector::create(
        Rcpp::_["value"] = value[0],
        Rcpp::_["derivative"] = jacobian[0]
    );
}