#ifndef ADOUBLE_CONVERTER_H
#define ADOUBLE_CONVERTER_H

#include <Rcpp.h>
#include <cppad/cppad.hpp>

typedef CppAD::AD<double> adouble;

namespace Rcpp {
  namespace traits {
  
    // Specialize Exporter for converting SEXP to adouble in the traits namespace
    template <> 
    struct Exporter<adouble> {
        // The underlying value
        adouble value;
        
        // Constructor that extracts the adouble from SEXP
        Exporter(SEXP x) {
            if (!Rf_isReal(x) || Rf_length(x) != 1) {
                Rcpp::stop("Expected a numeric scalar for adouble conversion");
            }
            // Convert the R numeric scalar to double and then to adouble.
            value = adouble(REAL(x)[0]);
        }
        
        // Conversion operator that returns the adouble value.
        operator adouble() const {
            return value;
        }
    };
    
  } // namespace traits

    // Specialize wrap to convert an adouble back to SEXP.
    template<> 
    SEXP wrap(const adouble& a) {
        // Use CppAD::Value to extract the underlying double value.
        double v = CppAD::Value(a);
        return Rcpp::wrap(v);
    }
}

#endif // ADOUBLE_CONVERTER_H