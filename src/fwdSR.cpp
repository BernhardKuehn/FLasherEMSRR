/*
 * Copyright 2014 FLR Team. Distributed under the GPL 2 or later
 * Maintainer: Finlay Scott, JRC
 * Stock–recruitment models for FLasherEMSRR.
 * Native models retain CppAD derivatives; other model expressions are
 * evaluated by the existing R-backed customSRR() callback.
 */

#include "../inst/include/fwdSR.h"
#include "../inst/include/adouble_converter.h" // used for converting SEXP types to adouble
#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// Native functions use the existing three-argument function-pointer type.
// They ignore model_name; customSRR uses it as an R expression.
// Declare here so the map below can refer to definitions later in this file.
template <typename T>
T ricker(const T srp, const std::vector<double> params,
         const std::string model_name);
template <typename T>
T bevholt(const T srp, const std::vector<double> params,
          const std::string model_name);
template <typename T>
T bevholtDa(const T srp, const std::vector<double> params,
            const std::string model_name);
template <typename T>
T constant(const T srp, const std::vector<double> params,
           const std::string model_name);
template <typename T>
T bevholtSS3(const T srp, const std::vector<double> params,
             const std::string model_name);
template <typename T>
T cushing(const T srp, const std::vector<double> params,
          const std::string model_name);
template <typename T>
T segreg(const T srp, const std::vector<double> params,
         const std::string model_name);
template <typename T>
T segregDa(const T srp, const std::vector<double> params,
           const std::string model_name);
template <typename T>
T survsrr(const T srp, const std::vector<double> params,
          const std::string model_name);
template <typename T>
T bevholtsig(const T srp, const std::vector<double> params,
             const std::string model_name);
template <typename T>
T mixedsrr(const T srp, const std::vector<double> params,
           const std::string model_name);

namespace {

// Prevent out-of-bounds access when an SRR is supplied with too few params.
void require_params(const std::vector<double>& p, std::size_t n,
                    const char* name) {
  if (p.size() < n) {
    Rcpp::stop("%s needs at least %d parameter(s), but got %d",
               name, static_cast<int>(n), static_cast<int>(p.size()));
  }
}
             
             // Used for dispatch ONLY. Do not remove whitespace from an R expression
             // passed to the callback: it could occur inside a string literal.
             std::string compact_model(std::string name) {
               name.erase(std::remove_if(name.begin(), name.end(),
                                         [](unsigned char c) { return std::isspace(c) != 0; }), name.end());
               if (name.compare(0, 4, "rec~") == 0)
                 name.erase(0, 4);
               return name;
             }
             
           } // namespace

/*-------------------------------------------------*/
// Templated class

/*! \brief Initialise the model map of SR functions. */
template <typename T>
void fwdSR_base<T>::init_model_map() {
  // Original FLasher names and aliases.
  map_model_name_to_function["ricker"] = &ricker<T>;
  map_model_name_to_function["Ricker"] = &ricker<T>;
  map_model_name_to_function["bevholt"] = &bevholt<T>;
  map_model_name_to_function["Bevholt"] = &bevholt<T>;
  map_model_name_to_function["bevholtDa"] = &bevholtDa<T>;
  map_model_name_to_function["bevholtda"] = &bevholtDa<T>;
  map_model_name_to_function["constant"] = &constant<T>;
  map_model_name_to_function["Constant"] = &constant<T>;
  map_model_name_to_function["mean"] = &constant<T>;
  map_model_name_to_function["Mean"] = &constant<T>;
  map_model_name_to_function["geomean"] = &constant<T>;
  map_model_name_to_function["Geomean"] = &constant<T>;
  map_model_name_to_function["bevholtSS3"] = &bevholtSS3<T>;
  map_model_name_to_function["bevholtss3"] = &bevholtSS3<T>;
  map_model_name_to_function["BevholtSS3"] = &bevholtSS3<T>;
  map_model_name_to_function["cushing"] = &cushing<T>;
  map_model_name_to_function["Cushing"] = &cushing<T>;
  map_model_name_to_function["segreg"] = &segreg<T>;
  map_model_name_to_function["Segreg"] = &segreg<T>;
  map_model_name_to_function["segregDa"] = &segregDa<T>;
  map_model_name_to_function["SegregDa"] = &segregDa<T>;
  map_model_name_to_function["survSRR"] = &survsrr<T>;
  map_model_name_to_function["survsrr"] = &survsrr<T>;
  map_model_name_to_function["bevholtsig"] = &bevholtsig<T>;
  map_model_name_to_function["Bevholtsig"] = &bevholtsig<T>;
  map_model_name_to_function["mixedsrr"] = &mixedsrr<T>;
  map_model_name_to_function["mixedSRR"] = &mixedsrr<T>;
  
  // This branch passes formula RHS text from R/coerce.R, not always
  // the canonical names above. Match known RHSs exactly.
  map_model_name_to_function["a"] = &constant<T>;
  map_model_name_to_function["a*ssb*exp(-b*ssb)"] = &ricker<T>;
  map_model_name_to_function["a*ssb/(b+ssb)"] = &bevholt<T>;
  map_model_name_to_function["a/(1+(b/ssb)^d)"] = &bevholtDa<T>;
  map_model_name_to_function["a*ssb^b"] = &cushing<T>;
  map_model_name_to_function["ifelse(ssb<=b,a*ssb,a*b)"] = &segreg<T>;
  map_model_name_to_function["ifelse(ssb^d<=b,a*ssb^d,a*b)"] =
    &segregDa<T>;
    
    map_model_name_to_function["customSRR"] = &customSRR<T>;
    map_model_name_to_function["CustomSRR"] = &customSRR<T>;
}

/*! \brief Empty constructor. */
template <typename T>
fwdSR_base<T>::fwdSR_base() : model(nullptr), deviances_mult(true) {
}

/*! \brief Main constructor for the fwdSR class
 *
 * Sets the function pointer to point at the correct SR function.
 * Assumes all dimensions (e.g. of deviances and parameters) have been checked in R (no checks made here).
 *
 * \param model_name The name of the SR function (must be present in the model map else it fails).
 * \param params_ip The SR parameters. The parameters are stored in the first dimension and can be disaggregated by time, area etc.
 * \param deviances_ip Residuals that can be applied to the predicted recruitment.
 * \param deviances_mult_ip Are the deviances multiplicative (true) or additive (false).
 */
template <typename T>
fwdSR_base<T>::fwdSR_base(const std::string model_name_ip,
                          const FLQuant params_ip,
                          const FLQuant deviances_ip,
                          const bool deviances_mult_ip)
  : model(nullptr), model_name(model_name_ip), params(params_ip),
    deviances(deviances_ip), deviances_mult(deviances_mult_ip) {
  
  init_model_map();
  
  const std::string key = compact_model(model_name_ip);
  if (key.empty())
    Rcpp::stop("SRR model name/formula RHS is empty");
  
  typename model_map_type::const_iterator model_pair_found =
    map_model_name_to_function.find(key);
  
  if (model_pair_found != map_model_name_to_function.end() &&
      key != "customSRR" && key != "CustomSRR") {
    model = model_pair_found->second;
    Rcpp::Rcout << "FLasherEMSRR SRR: native model [" << key << "]\n";
  } else {
    // A bare "customSRR" name is not an expression to evaluate.
    if (key == "customSRR" || key == "CustomSRR") {
      Rcpp::stop("Specify a custom SRR expression, e.g. "
                   "rec ~ gampredict(ssb, params), not just customSRR");
    }
    
    // Check the fallback is registered instead of assuming it exists.
    typename model_map_type::const_iterator fallback =
      map_model_name_to_function.find("customSRR");
    if (fallback == map_model_name_to_function.end() ||
        fallback->second == nullptr) {
      Rcpp::stop("SRR '%s' is not native and customSRR is unavailable",
                 model_name_ip.c_str());
    }
    
    model = fallback->second;
    Rcpp::Rcout << "FLasherEMSRR SRR: R custom expression ["
                << model_name_ip << "]\n";
  }
}

/*! \brief Intrusive wrap, primarily for tests. */
template <typename T>
fwdSR_base<T>::operator SEXP() const {
  return Rcpp::List::create(
    Rcpp::Named("params", params),
    Rcpp::Named("model_name", model_name),
    Rcpp::Named("deviances", deviances),
    Rcpp::Named("deviances_mult", deviances_mult));
}

/*! \brief Copy constructor. */
template <typename T>
fwdSR_base<T>::fwdSR_base(const fwdSR_base<T>& source) {
  model = source.model;
  model_name = source.model_name;
  params = source.params;
  deviances = source.deviances;
  deviances_mult = source.deviances_mult;
  map_model_name_to_function = source.map_model_name_to_function;
}

/*! \brief Assignment operator. */
template <typename T>
fwdSR_base<T>& fwdSR_base<T>::operator=(const fwdSR_base<T>& source) {
  if (this != &source) {
    model = source.model;
    model_name = source.model_name;
    params = source.params;
    deviances = source.deviances;
    deviances_mult = source.deviances_mult;
    map_model_name_to_function = source.map_model_name_to_function;
  }
  return *this;
}


/*! \name Get the SR parameters
 *
 * Given the year, unit, season, area and iter, returns the corresponding stock recruitment parameters.
 * \param year The year of the SR parameters to use.
 * \param unit The unit of the SR parameters to use.
 * \param season The season of the SR parameters to use.
 * \param area The area of the SR parameters to use.
 * \param iter The iter of the SR parameters to use.
 */
template <typename T>
std::vector<double> fwdSR_base<T>::get_params(
    unsigned int year, unsigned int unit, unsigned int season,
    unsigned int area, unsigned int iter) const {
  
  const int nparams = get_nparams();
  std::vector<double> model_params(nparams);
  
  // Preserve the existing dimension-recycling behaviour.
  if (year > params.get_nyear()) year = 1;
  if (unit > params.get_nunit()) unit = 1;
  if (season > params.get_nseason()) season = 1;
  if (area > params.get_narea()) area = 1;
  if (iter > params.get_niter()) iter = 1;
  
  for (int i = 1; i <= nparams; ++i)
    model_params[i - 1] =
      params(i, year, unit, season, area, iter);
  
  return model_params;
}

/*! \brief Evaluate the SR model for one SRP value. */
template <typename T>
T fwdSR_base<T>::eval_model(
    const T srp, int year, int unit, int season, int area, int iter,
    const std::string model_name_arg) const {
  
  std::vector<double> model_params =
    get_params(year, unit, season, area, iter);
  
  for (double param : model_params) {
    if (Rcpp::NumericVector::is_na(param)) {
      Rcpp::warning("An SR model param is NA. Setting rec to 0.");
      return T(0);
    }
  }
  if (model == nullptr)
    Rcpp::stop("SRR model function pointer is not initialised");
  
  return model(srp, model_params, model_name_arg);
}

/*! \brief Evaluate using a vector of five parameter indices. */
template <typename T>
T fwdSR_base<T>::eval_model(
    const T srp, const std::vector<unsigned int> indices,
    const std::string model_name_arg) const {
  
  if (indices.size() != 5)
    Rcpp::stop("In fwdSR::eval_model: params_indices must have length 5");
  
  return eval_model(srp, indices[0], indices[1], indices[2],
                    indices[3], indices[4], model_name_arg);
}

/*! \brief Predict recruitment
 *
 * Calculates the recruitment from an FLQuant of SRP, including the application of deviances.
 * The SRP can be a subset of the 'full' model SRP
 * (e.g. can be only one season out of all seasons, or several years out of all years).
 * It is therefore necessary to also pass in a vector of indices to specify the start
 * position of the SR params and deviances because we don't know the start position of the SRP argument relative to the whole operating model. 
 * i.e. we need to know the index of the params and deviances that correspond with the first value in the SRP vector.
 * The parameters and deviances dimensions are in line with the recruitment.
 * e.g. parameters in year 2, season 1 are used to calculated recruitment in year 2, season 1 given the SSB NOT applied to the SSB in year 2, season 1 to calculate recruitment the following year.
 * Internally, the method calls the eval() method to calculate deterministic recruitment then applies the deviances.
 *
 * \param srp The spawning reproductive potential that produces the recruitment.
 * \param initial_params_indices A vector of length 5 (year, unit, ... iter) to specify the start position of the indices of the SR params and deviances relative to the 'whole' operating model (starting at 1).
 */
template <typename T>
FLQuant_base<T> fwdSR_base<T>::predict_recruitment(const FLQuant_base<T> srp, const std::vector<unsigned int> initial_params_indices,
                                                   const std::string model_name){ 
  if (initial_params_indices.size() != 5){
    Rcpp::stop("In fwdSR::predict_recruitment. initial_params_indices must be of length 5.\n");
  }
  std::vector<unsigned int> srp_dim = srp.get_dim();
  if (srp_dim[0] != 1){
    Rcpp::stop("In fwdSR::predict_recruitment. srp must be of length 1 in the first dimension.\n");
  }
  // If deviances starting from the initial_params_indices are too small for the SRP it does not recycle = error
  // (e.g. if you pass in subset of SRP with 10 out 20 years but the initial_params_indices has year = 15, i.e. 5 years too short in deviances)
  // Iters are OK
  std::vector<unsigned int> res_dim = deviances.get_dim();
  // Distance between initial_params_indices and res_dim cannot be smaller than srp_dim
  for (unsigned int dim_counter = 1; dim_counter <= 4; ++dim_counter){
    if((res_dim[dim_counter] - initial_params_indices[dim_counter-1] + 1) < srp_dim[dim_counter]){
      Rcpp::stop("In fwdSR::predict_recruitment. Initial indices of deviances too small to cover the SRP\n");
    }
  }
  // Empty output object
  FLQuant_base<T> rec = srp;
  rec.fill(0.0);
  
  // SET sex ratio
  T sratio = 1.0;
  
  // GET unit dimnames, sorted
  std::vector<std::string> unit_names;
  unit_names = Rcpp::as<std::vector<std::string> >(deviances.get_dimnames()[2]);
  std::sort(unit_names.begin(), unit_names.end());
  std::vector<std::string> sex = { "F", "M" };
  
  // IF two unit(s) & unit_names %in% c('F','M'), sratio=0.5
  if (res_dim[2] == 2) {
    if(unit_names == sex) {
      sratio = 0.5;
    }
  } 
  
  // Going to have to loop over the dimensions and update the params and deviances indices - not nice
  std::vector<unsigned int> params_indices = initial_params_indices;
  //std::vector<unsigned int> deviances_indices = initial_deviances_indices;
  for (unsigned int year_counter = 1; year_counter <= srp_dim[1]; ++year_counter){
    params_indices[0] = initial_params_indices[0] + year_counter - 1;
    //deviances_indices[0] = initial_deviances_indices[0] + year_counter - 1;
    for (unsigned int unit_counter = 1; unit_counter <= srp_dim[2]; ++unit_counter){
      params_indices[1] = initial_params_indices[1] + unit_counter - 1;
      //deviances_indices[1] = initial_deviances_indices[1] + unit_counter - 1;
      for (unsigned int season_counter = 1; season_counter <= srp_dim[3]; ++season_counter){
        params_indices[2] = initial_params_indices[2] + season_counter - 1;
        //deviances_indices[2] = initial_deviances_indices[2] + season_counter - 1;
        for (unsigned int area_counter = 1; area_counter <= srp_dim[4]; ++area_counter){
          params_indices[3] = initial_params_indices[3] + area_counter - 1;
          //deviances_indices[3] = initial_deviances_indices[3] + area_counter - 1;
          for (unsigned int iter_counter = 1; iter_counter <= srp_dim[5]; ++iter_counter){
            params_indices[4] = initial_params_indices[4] + iter_counter - 1;
            //deviances_indices[4] = initial_deviances_indices[4] + iter_counter - 1;
            T rec_temp = eval_model(srp(1, year_counter, unit_counter, season_counter, area_counter, iter_counter),
                                    params_indices,
                                    model_name);
            rec_temp *= sratio;
            if (deviances_mult == true){
              rec_temp *= deviances(1, params_indices[0], params_indices[1], params_indices[2], params_indices[3], params_indices[4]);
            }
            else {
              rec_temp += deviances(1, params_indices[0], params_indices[1], params_indices[2], params_indices[3], params_indices[4]);
            }
            rec(1, year_counter, unit_counter, season_counter, area_counter, iter_counter) = rec_temp;
          }}}}}
  return rec;
}

template <typename T>
FLQuant_base<double> fwdSR_base<T>::get_params() const {
  return params;
}

template <typename T>
std::string fwdSR_base<T>::get_model_name() const {
  return model_name;
}

template <typename T>
int fwdSR_base<T>::get_nparams() const {
  return params.get_nquant();
}

template <typename T>
FLQuant_base<double> fwdSR_base<T>::get_deviances() const {
  return deviances;
}

template <typename T>
bool fwdSR_base<T>::get_deviances_mult() const {
  return deviances_mult;
}

template <typename T>
void fwdSR_base<T>::set_deviances(
    const FLQuant_base<double> new_deviances) {
  deviances = new_deviances;
}

template <typename T>
void fwdSR_base<T>::set_deviances_mult(
    const bool new_deviances_mult) {
  deviances_mult = new_deviances_mult;
}

/*! \brief Whether recruitment occurs for a unit in this timestep. */
template <typename T>
bool fwdSR_base<T>::does_recruitment_happen(
    unsigned int unit, unsigned int year, unsigned int season) const {
  std::vector<double> sr_params =
    get_params(year, unit, season, 1, 1);
  return !sr_params.empty() &&
    !Rcpp::NumericVector::is_na(sr_params[0]);
}

/*! \brief Whether recruitment has occurred by the given season. */
template <typename T>
bool fwdSR_base<T>::has_recruitment_happened(
    unsigned int unit, unsigned int year, unsigned int season) const {
  unsigned int nseasons = params.get_nseason();
  if (nseasons == 0 || season == 0 || season > nseasons)
    Rcpp::stop("Invalid season in has_recruitment_happened");
  
  for (unsigned int s = 1; s <= season; ++s) {
    if (does_recruitment_happen(unit, year, s))
      return true;
  }
  return false;
}

// Explicit instantiation of the class.
template class fwdSR_base<double>;
template class fwdSR_base<adouble>;

// //--------------------------------------------------------------------
// // SRR functions
// // The params in these functions do not have any disaggregation (e.g. by time or area)
// // They are only the params required to evaluate the SRR function at a specific point
// // The disaggregated parameter values are stored in the fwdSR_base class as an FLQuant_base object
// // The fwdSR.eval_model method sorts out time step of params etc and passes the correct set of params to these functions
// // These functions must have the same argument list so that it matches the typedef for the model pointer
template <typename T>
T customSRR(const T srp, const std::vector<double> params,
            const std::string model_name_arg) {
  if (model_name_arg.empty())
    Rcpp::stop("customSRR received an empty model expression");
  
  // The child environment prevents assignment to .GlobalEnv, but
  // inherits user-defined functions such as gampredict().
  Rcpp::Environment env =
    Rcpp::Environment::global_env().new_child(true);
  Rcpp::NumericVector r_params(params.begin(), params.end());
  env["ssb"] = srp;
  env["params"] = r_params;
  
  try {
    Rcpp::Function parseFunc =
      Rcpp::Environment::base_env()["parse"];
    SEXP expressions = parseFunc(
      Rcpp::_["text"] = model_name_arg);
    
    if (Rf_length(expressions) != 1)
      Rcpp::stop("customSRR requires exactly one R expression");
    
    SEXP result = Rcpp_eval(VECTOR_ELT(expressions, 0), env);
    const double value = Rcpp::as<double>(result);
    if (!std::isfinite(value))
      Rcpp::stop("customSRR returned a non-finite recruitment");
    
    // An R numeric result does NOT retain CppAD derivatives.
    return T(value);
  } catch (const std::exception& e) {
    Rcpp::stop("customSRR failed for '%s': %s",
               model_name_arg.c_str(), e.what());
  }
  return T(0); // Unreachable; makes return paths explicit.
}

template <typename T>
T ricker(const T srp, const std::vector<double> params,
         const std::string) {
  require_params(params, 2, "ricker");
  return params[0] * srp * exp(-params[1] * srp);
}

template <typename T>
T bevholt(const T srp, const std::vector<double> params,
          const std::string) {
  require_params(params, 2, "bevholt");
  if (params.size() > 2)
    return params[0] /
      (1.0 + pow(params[1] / srp, params[2]));
  return params[0] * srp / (params[1] + srp);
}

template <typename T>
T bevholtDa(const T srp, const std::vector<double> params,
            const std::string) {
  require_params(params, 3, "bevholtDa");
  return params[0] /
    (1.0 + pow(params[1] / srp, params[2]));
}

template <typename T>
T constant(const T, const std::vector<double> params,
           const std::string) {
  require_params(params, 1, "constant");
  return T(params[0]);
}

template <typename T>
T bevholtSS3(const T srp, const std::vector<double> params,
             const std::string) {
  require_params(params, 3, "bevholtSS3");
  const double s = params[0];
  const double R0 = params[1];
  const double v = params[2];
  return (4.0 * s * R0 * srp) /
    (v * (1.0 - s) + srp * (5.0 * s - 1.0));
}

template <typename T>
T cushing(const T srp, const std::vector<double> params,
          const std::string) {
  require_params(params, 2, "cushing");
  return params[0] * exp(log(srp) * params[1]);
}

template <typename T>
T segreg(const T srp, const std::vector<double> params,
         const std::string) {
  require_params(params, 2, "segreg");
  return srp <= params[1]
  ? params[0] * srp : T(params[0] * params[1]);
}

// CppAD must record both branches, not select one while taping.
template <>
adouble segreg<adouble>(const adouble srp,
                        const std::vector<double> params,
                        const std::string) {
  require_params(params, 2, "segreg");
  return CppAD::CondExpLe(
    srp, adouble(params[1]), params[0] * srp,
    adouble(params[0] * params[1]));
}

template <typename T>
T segregDa(const T srp, const std::vector<double> params,
           const std::string) {
  require_params(params, 3, "segregDa");
  const T powered = pow(srp, params[2]);
  return powered <= params[1]
  ? params[0] * powered : T(params[0] * params[1]);
}

template <>
adouble segregDa<adouble>(const adouble srp,
                          const std::vector<double> params,
                          const std::string) {
  require_params(params, 3, "segregDa");
  const adouble powered = pow(srp, params[2]);
  return CppAD::CondExpLe(
    powered, adouble(params[1]), params[0] * powered,
    adouble(params[0] * params[1]));
}

template <typename T>
T survsrr(const T ssf, const std::vector<double> params,
          const std::string) {
  require_params(params, 4, "survsrr");
  const double R0 = params[0];
  const double sfrac = params[1];
  const double beta = params[2];
  const double SB0 = params[3];
  const double sratio =
    params.size() > 4 ? params[4] : 0.5;
  const double z0 = log(1.0 / (SB0 / R0));
  const double zmax = z0 + sfrac * (0.0 - z0);
  const T zsurv = exp(
    (1.0 - pow(ssf / SB0, beta)) * (zmax - z0) + z0);
  return ssf * zsurv * sratio;
}

template <typename T>
T bevholtsig(const T srp, const std::vector<double> params,
             const std::string) {
  require_params(params, 3, "bevholtsig");
  return params[0] /
    (pow(params[1] / srp, params[2]) + 1.0);
}

template <typename T>
T mixedsrr(const T srp, const std::vector<double> params,
           const std::string) {
  require_params(params, 3, "mixedsrr");
  const double selector = params[2];
  if (!std::isfinite(selector) ||
      selector != std::floor(selector) ||
      selector < 1 || selector > 3)
    Rcpp::stop("mixedsrr selector must be 1, 2, or 3");
  
  // Do not call bevholt(params) here: its optional third parameter
  // would change the formula for a three-parameter mixedsrr.
  if (selector == 1)
    return params[0] * srp / (params[1] + srp);
  if (selector == 2)
    return params[0] * srp * exp(-params[1] * srp);
  return segreg(srp, params, "");
}

// Native models are used through pointers declared in this translation
// unit. Keep explicit instantiations for the externally declared callback.
template double customSRR<double>(
    const double, const std::vector<double>, const std::string);
template adouble customSRR<adouble>(
    const adouble, const std::vector<double>, const std::string);