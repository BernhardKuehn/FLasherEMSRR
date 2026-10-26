#ifndef FLASHER_EMSRR_CUSTOM_SRR_ATOMIC_H
#define FLASHER_EMSRR_CUSTOM_SRR_ATOMIC_H

#include <Rcpp.h>
#include "cppad/cppad.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace emsrr {

// Evaluate the same R expression as the existing customSRR callback.
// The model must be deterministic and must return one finite number.
inline double evaluate_r_srr(
    double ssb,
    const std::vector<double>& params,
    const std::string& expression
) {
    if (!std::isfinite(ssb) || ssb < 0.0)
        Rcpp::stop("customSRR: SSB must be finite and nonnegative");

    try {
        Rcpp::Environment env =
            Rcpp::Environment::global_env().new_child(true);

        env["ssb"] = ssb;
        env["params"] =
            Rcpp::NumericVector(params.begin(), params.end());

        Rcpp::Function parse =
            Rcpp::Environment::base_env()["parse"];
        Rcpp::ExpressionVector expressions =
            parse(Rcpp::_["text"] = expression);

        if (expressions.size() != 1)
            Rcpp::stop("Expected exactly one SRR expression");

        Rcpp::RObject result =
            Rcpp::Rcpp_eval(expressions[0], env);

        // Also supports the scalar FLQuant returned by your wrapper.
        const double value = Rcpp::as<double>(result);

        if (!std::isfinite(value) || value < 0.0)
            Rcpp::stop("Expected finite, nonnegative recruitment");

        return value;
    } catch (const std::exception& e) {
        Rcpp::stop("customSRR '%s', SSB=%.17g: %s",
                   expression.c_str(), ssb, e.what());
    }

    return 0.0; // Unreachable.
}

class RStockRecruit : public CppAD::atomic_base<double> {
private:
    const std::vector<double> params_;
    const std::string expression_;

    double value(double ssb) const {
        return evaluate_r_srr(ssb, params_, expression_);
    }

    double derivative(double ssb) const {
        const double h = 1e-5 * std::max(1.0, std::abs(ssb));
        double slope;

        if (ssb >= h) {
            slope = (value(ssb + h) - value(ssb - h)) / (2.0 * h);
        } else {
            // Second-order forward difference near the SSB=0 boundary.
            slope = (-3.0 * value(ssb) + 4.0 * value(ssb + h)
                     - value(ssb + 2.0 * h)) / (2.0 * h);
        }

        if (!std::isfinite(slope))
            Rcpp::stop("customSRR: non-finite numerical SSB derivative");

        return slope;
    }

public:
    RStockRecruit(
        const std::vector<double>& params,
        const std::string& expression
    ) :
        CppAD::atomic_base<double>(
            "RStockRecruit",
            CppAD::atomic_base<double>::bool_sparsity_enum
        ),
        params_(params),
        expression_(expression) {}

    bool forward(
        size_t p, size_t q,
        const CppAD::vector<bool>& vx,
        CppAD::vector<bool>& vy,
        const CppAD::vector<double>& tx,
        CppAD::vector<double>& ty
    ) override {
        // Support values and first derivatives, not higher orders.
        if (q > 1)
            Rcpp::stop("customSRR atomic supports first derivatives only");

        // One scalar input: SSB. Its Taylor coefficients are tx[0], tx[1].
        if (vx.size() != 0)
            vy[0] = vx[0];

        if (p == 0)
            ty[0] = value(tx[0]);

        if (q == 1)
            ty[1] = derivative(tx[0]) * tx[1];

        return true;
    }

    bool reverse(
        size_t q,
        const CppAD::vector<double>& tx,
        const CppAD::vector<double>& /*ty*/,
        CppAD::vector<double>& px,
        const CppAD::vector<double>& py
    ) override {
        // q=0 here supplies a first-order reverse derivative.
        if (q != 0)
            Rcpp::stop("customSRR atomic supports first derivatives only");

        px[0] = py[0] * derivative(tx[0]);
        return true;
    }

    bool for_sparse_jac(
        size_t q,
        const CppAD::vector<bool>& r,
        CppAD::vector<bool>& s,
        const CppAD::vector<double>& /*x*/
    ) override {
        // Structurally, recruitment may depend on SSB even if the
        // numerical slope happens to be zero at this particular point.
        for (size_t k = 0; k < q; ++k)
            s[k] = r[k];
        return true;
    }

    bool rev_sparse_jac(
        size_t q,
        const CppAD::vector<bool>& rt,
        CppAD::vector<bool>& st,
        const CppAD::vector<double>& /*x*/
    ) override {
        for (size_t k = 0; k < q; ++k)
            st[k] = rt[k];
        return true;
    }
};

// Ownership must extend beyond individual customSRR() calls.
class AtomicScope;

inline AtomicScope*& active_scope() {
    static AtomicScope* scope = nullptr;
    return scope;
}

class AtomicScope {
private:
    AtomicScope* previous_;
    std::vector<std::unique_ptr<RStockRecruit>> operations_;

public:
    AtomicScope() : previous_(active_scope()) {
        active_scope() = this;
    }

    ~AtomicScope() {
        active_scope() = previous_;
        // operations_ is destroyed after the projection's tapes.
    }

    AtomicScope(const AtomicScope&) = delete;
    AtomicScope& operator=(const AtomicScope&) = delete;

    RStockRecruit& make(
        const std::vector<double>& params,
        const std::string& expression
    ) {
        operations_.emplace_back(new RStockRecruit(params, expression));
        return *operations_.back();
    }
};

inline CppAD::AD<double> evaluate_ad_srr(
    const CppAD::AD<double>& ssb,
    const std::vector<double>& params,
    const std::string& expression
) {
    // No derivative connection is needed for an AD parameter.
    if (CppAD::Parameter(ssb))
        return evaluate_r_srr(CppAD::Value(ssb), params, expression);

    if (active_scope() == nullptr)
        Rcpp::stop("customSRR atomic needs an active emsrr::AtomicScope");

    RStockRecruit& operation = active_scope()->make(params, expression);

    CppAD::vector<CppAD::AD<double>> ax(1), ay(1);
    ax[0] = ssb;
    operation(ax, ay);
    return ay[0];
}

} // namespace emsrr
#endif