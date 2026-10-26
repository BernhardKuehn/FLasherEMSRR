#include <Rcpp.h>
#include "../inst/include/FLasherEMSRR.h"
#include "../inst/include/custom_srr_atomic.h"

// Diagnostic: check the value and SSB derivative of an R-backed SRR.
// [[Rcpp::export]]
Rcpp::NumericVector customSRR_AD_diagnostic(
    double ssb,
    std::string model_name,
    Rcpp::NumericVector params
) {
  using AD = CppAD::AD<double>;
  
  // Must outlive every tape referencing its atomic operations.
  // Declared before f, so f is destroyed before atomic_scope.
  emsrr::AtomicScope atomic_scope;
  CppAD::ADFun<double> f;
  
  std::vector<AD> ax(1), ay(1);
  std::vector<double> p(params.begin(), params.end());
  bool recording = false;
  
  try {
    ax[0] = ssb;
    CppAD::Independent(ax);
    recording = true;
    
    ay[0] = customSRR<AD>(ax[0], p, model_name);
    
    f.Dependent(ax, ay);
    recording = false;
    
    std::vector<double> x0(1, ssb);
    std::vector<double> value = f.Forward(0, x0);
    std::vector<double> jacobian = f.Jacobian(x0);
    std::vector<double> sparse_jacobian = f.SparseJacobian(x0);
    
    return Rcpp::NumericVector::create(
      Rcpp::_["value"] = value[0],
                              Rcpp::_["derivative"] = jacobian[0],
                                                              Rcpp::_["sparse_derivative"] = sparse_jacobian[0]
    );
  } catch (...) {
    // Do not leave an active recording after a failed R callback.
    if (recording)
      AD::abort_recording();
    throw;
  }
}