/* 
 * Copyright 2014 FLR Team. Distributed under the GPL 2 or later
 * Maintainer: Finlay Scott, JRC
 */

#include "../inst/include/fwdSR.h"
#include "../inst/include/adouble_converter.h" // used for converting SEXP types to adouble
#include <Rcpp.h>
using namespace Rcpp;

/*-------------------------------------------------*/
// Templated class

/*! \brief Initialises the model map of SR functions
 *
 * The model map maps the names of SR functions to the addresses of the SR functions contained in the package.
 * This method initialises it.
 */
template <typename T>
void fwdSR_base<T>::init_model_map(){
    // Fill up the map
    map_model_name_to_function["customSRR"] = &customSRR;
    map_model_name_to_function["CustomSRR"] = &customSRR;
    return;
}

/*! \brief Empty constructor that creates empty members
 */
template <typename T>
fwdSR_base<T>::fwdSR_base(){
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
fwdSR_base<T>::fwdSR_base(const std::string model_name_ip, const FLQuant params_ip, const FLQuant deviances_ip, const bool deviances_mult_ip) {
    model_name = model_name_ip;
    // std::cout << model_name; // print for debugging
    params = params_ip;
    deviances = deviances_ip;
    deviances_mult = deviances_mult_ip;
    init_model_map();
    // Set the model pointer
    // always map to the customSRR-function as nothing else is defined
    typename model_map_type::const_iterator model_pair_found = map_model_name_to_function.find("customSRR");
    //std::cout << model_pair_found; // print for debugging
    if (model_pair_found != map_model_name_to_function.end()){
        model = model_pair_found->second; // pulls out value - the address of the SR function
    }
    else
        Rcpp::stop("SRR model not found\n");
} 

/*! \brief Intrusive 'wrap' for the fwdSR class.
 * Returns a List of stuff - used for tests, not really used for much else.
 */
template <typename T>
fwdSR_base<T>::operator SEXP() const{
    return Rcpp::List::create(Rcpp::Named("params", params),
                            Rcpp::Named("model_name", model_name),
                            Rcpp::Named("deviances",deviances),
                            Rcpp::Named("deviances_mult",deviances_mult));
}


/*! \brief Copy constructor for the fwdSR class.
 *
 * Ensures a deep copy.
 *
 * \param fwdSR_source The fwdSR to be copied.
 */
template <typename T>
fwdSR_base<T>::fwdSR_base(const fwdSR_base<T>& fwdSR_source){
    model = fwdSR_source.model; // Copy the pointer - we want it to point to the same place so copying should be fine.
    model_name = fwdSR_source.model_name;
    params = fwdSR_source.params;
    deviances = fwdSR_source.deviances;
    deviances_mult = fwdSR_source.deviances_mult;
}

/*! \brief Assignment operator for the fwdSR class.
 *
 * Ensures a deep copy.
 *
 * \param fwdSR_source The fwdSR to be copied.
 */
template <typename T>
fwdSR_base<T>& fwdSR_base<T>::operator = (const fwdSR_base<T>& fwdSR_source){
	if (this != &fwdSR_source){
        model = fwdSR_source.model; // Copy the pointer - we want it to point to the same place so copying should be fine.
        model_name = fwdSR_source.model_name;
        params = fwdSR_source.params;
        deviances = fwdSR_source.deviances;
        deviances_mult = fwdSR_source.deviances_mult;
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
std::vector<double> fwdSR_base<T>::get_params(unsigned int year, unsigned int unit, unsigned int season, unsigned int area, unsigned int iter) const{
    const int nparams = get_nparams();
    std::vector<double> model_params(nparams);
    // Parameters get recycled, i.e.  if requested year is bigger than years in the params FLQuant then we just take the first one.
    if (year > params.get_nyear()){
        year = 1;
    }
    if (unit > params.get_nunit()){
        unit = 1;
    }
    if (season > params.get_nseason()){
        season = 1;
    }
    if (area > params.get_narea()){
        area = 1;
    }
    if (iter > params.get_niter()){
        iter = 1;
    }
    for (int i = 1; i <= nparams; ++i){
        model_params[i-1] = params(i,year,unit,season,area,iter);
  // uncomment for debugging: 
  // Rprintf("Year: %i Unit: %i Season: %i Iter: %i Param %i: %f\n",  year, unit, season, iter, i, model_params[i-1]);
    }
    return model_params;
}

/*! \name Evaluate the SR model
 *
 * Produces a single value of recruitment given a single value of the SRP.
 * The SR parameters can vary with time, area etc. so it is necessary to also pass in the indices (starting from 1) for the parameters.
 * If the parameters are fixed in any dimension, then the index should be 1 for that dimension.
 * If the index passed in is greater than that dimension in the parameter FLQuant a value of 1 will be used instead.
 *
 */
//@{
/*!
 * \brief Evaluates the SR model
 * \param srp The spawning reproductive potential that produces the recruitment.
 * \param year The year of the SR parameters to use.
 * \param unit The unit of the SR parameters to use.
 * \param season The season of the SR parameters to use.
 * \param area The area of the SR parameters to use.
 * \param iter The iter of the SR parameters to use.
 */
template <typename T>
T fwdSR_base<T>::eval_model(const T srp, int year, int unit, int season, int area, int iter, const std::string model_name) const{
    // Get the parameters
    std::vector<double> model_params = get_params(year, unit, season, area, iter);
    
    // uncomment for debugging: 
    // Rprintf("Params: %f\n", model_params);
    // Check if any params are NA - if so, don't evaluate model, set rec to 0.0 for clean exit
    // Check first value to see if bad for clean exit
    //Rprintf("size model_params: %i \n", model_params.size());
    T rec = 0.0;
    bool paramNA = false;
    for (double param : model_params){
        if(Rcpp::NumericVector::is_na(param)){
            paramNA = true;
        }
    }
    if (paramNA == false){
      // evaluate with the model name passed as function argument
        rec = model(srp, model_params,model_name);
    }
    else {
        Rcpp::warning("An SR model param is NA. Setting rec to 0 else something bad will happen.\n");
    }
    return rec;
}
/*!
 * \brief Evaluates the SR model
 * \param srp The spawning reproductive potential that produces the recruitment.
 * \param params_indices The indices of the SR params (starting at 1).
 */
template <typename T>
T fwdSR_base<T>::eval_model(const T srp, const std::vector<unsigned int> params_indices, const std::string model_name) const{ 
    // Check length of params_indices
    if (params_indices.size() != 5){
        Rcpp::stop("In fwdSR::eval_model. params_indices must be of length 5.");
    }
    T rec = eval_model(srp, params_indices[0], params_indices[1], params_indices[2], params_indices[3], params_indices[4],
                       model_name);
    return rec;
}
//@}

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
FLQuant_base<double> fwdSR_base<T>::get_params() const{
    return params;
}

template <typename T>
std::string fwdSR_base<T>::get_model_name() const{
    return model_name;
}

/*! \brief Number of SR parameters.
 *
 * The number of SR parameters is the length of the first dimension of the parameters FLQuant.
 */
template <typename T>
int fwdSR_base<T>::get_nparams() const{
    return params.get_nquant();
}

/*! \brief Get the deviances.
 *
 * Returns the deviances.
 */
template <typename T>
FLQuant_base<double> fwdSR_base<T>::get_deviances() const{
    return deviances;
}

/*! \brief Get the deviances multiplier.
 *
 * Returns the deviances multiplier.
 */
template <typename T>
bool fwdSR_base<T>::get_deviances_mult() const{
    return deviances_mult;
}

template <typename T>
void fwdSR_base<T>::set_deviances(const FLQuant_base<double> new_deviances){
    deviances = new_deviances;
}

template <typename T>
void fwdSR_base<T>::set_deviances_mult(const bool new_deviances_mult){
    deviances_mult = new_deviances_mult;
}

//! Does recruitment happen for a unit in that timestep
/*!
  Each unit can recruit in a different season. Each unit can recruit only once per year.
  The first stock-recruitment parameter is checked.
  If it is NA, then recruitment does not happen.
  If the parameter is not NA, then recruitment happens.
  It is assumed that the timing pattern across iters is the same, e.g. if recruitment happens in season 1 for iter 1, it happens in season 1 for all iters.
  \param unit The unit we are checking for recruitment
  \param timestep The timestep we are checking for recruitment.
 */ 
template <typename T>
bool fwdSR_base<T>::does_recruitment_happen(unsigned int unit, unsigned int year, unsigned int season) const{
    unsigned int area = 1;
    // Just get first iter
    std::vector<double> sr_params = get_params(year, unit, season, area, 1); 
    bool did_spawning_happen = true;
    // Check the first parameter only
    if (Rcpp::NumericVector::is_na(sr_params[0])){
        did_spawning_happen = false;
    }
    return did_spawning_happen;
}

template <typename T>
bool fwdSR_base<T>::has_recruitment_happened(unsigned int unit, unsigned int year, unsigned int season) const{
    unsigned int area = 1;
    unsigned int has_spawning_happened = 0;

    // GET number of seasons
    unsigned int nseasons = params.get_nseason();

    // CREATE vector of length nseasons
    std::vector<int> spawns(nseasons);
    std::vector<double> sr_params = get_params(year, unit, 1, area, 1);
    spawns[0] = 0 + !Rcpp::NumericVector::is_na(sr_params[0]);

    // POPULATE with is.na(sr_params[0]), loop over 2:end
    for (unsigned int ns=1; ns < nseasons; ++ns){
      std::vector<double> sr_params = get_params(year, unit, ns+1, area, 1);
      spawns[ns] = spawns[ns-1] + !Rcpp::NumericVector::is_na(sr_params[0]);
    }
    
    has_spawning_happened = spawns[season-1] > 0;

    return has_spawning_happened;

}

// Explicit instantiation of class
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
T customSRR(const T srp, const std::vector<double> params,const std::string model_name){
    T rec;
  
  // uncommment for debugging
  // for(int i=0; i<params.size(); ++i){
  //   Rprintf("the value of params[%i] : %f \n", i, params[i]);
  // }
    // Use the global environment so that any function (such as predefined_function) is found.
    Environment env = Environment::global_env();
    
    // Insert the C++ arguments into the R environment.
    env.assign("ssb", srp);
    env.assign("params", params);
    
    // uncomment for debugging: 
    //double ssb_print = env["ssb"];
    //Rprintf("the value of SSB : %f \n", ssb_print);
    
    // Get R's built-in 'parse' function from the base namespace.
    Function parseFunc = Environment::base_env()["parse"];
    
    // Parse the provided formula call string into an R expression.
    SEXP exprList = parseFunc(_["text"] = model_name);
    
    // Evaluate the first expression from the parsed results in the global environment.
    SEXP result = Rcpp_eval(VECTOR_ELT(exprList, 0), env);
    // coerse to a double first
    double d = Rcpp::as<double>(result);
    
    // uncomment for debugging:
    //Rprintf("the value of Rec: %f \n", d);
    
    // than coersing it to the template format allowing to be cast to "adouble" correctly
    rec = T(d);
    
    // Convert and return the result as a double.
    return rec;
    //return static_cast<T>(REAL(rec)[0]);
    //return rec;
}
// 
// // create a 'custom' SRR that links to a R-call of a specified function
// 
// template <typename T>
// T customSRR(const T srp, const std::vector<double> params){
//   T rec;
//   // RInside R(srp, params, f_name);
//   // T rec;
//   // // call a specified R function from within C++
//   // // use params specified in the params slot as env. covars
//   // Rcpp::Function f(f_name)
//   // const T rec = f(ssb = srp,covars = params)
//   // R.parseEval(rec,ans);
//   // rec = double ans;
//   // return rec;
//   
//   // Rcpp::Function f(f_name);
//   
//   // Rcpp::NumericVector rec = f(srp,params);
//   // rec = f(srp,params);
//   return rec;
// }
// 
// 
// // Rcpp::NumericVector customSRR(const T srp, const std::vector<double> params, const char *f_name[]){
// //   // RInside R(srp, params, f_name);
// //   // call a specified R function from within C++
// //   // use params specified in the params slot as env. covars
// //   Rcpp::Function f(f_name);
// //   Rcpp::NumericVector rec = f(srp,params);
// //   return rec;
// // }
// 
// 
// template <typename T>
// T bevholt(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
// 
//     // rec = a * srp / (b + srp)
//     rec = params[0] * srp / (params[1] + srp);
// 
//     if (params.size() > 2) {
//       // rec = a / (1 + (b / srp) ^ d)
//       rec = params[0] / (1 + pow(params[1] / srp, params[2]));
//     }
//     
//     return rec;
// }
// 
// template <typename T>
// T bevholtDa(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
// 
//     // rec = a / (1 + (b / srp) ^ d)
//     rec = params[0] / (1 + pow(params[1] / srp, params[2]));
//     
//     return rec;
// }
// 
// template <typename T>
// T constant(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     // rec = a 
//     rec = params[0]; 
//     return rec;
// }
// 
// template <typename T>
// T bevholtSS3(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     // rec = (4 * s * R0 * ssb) / (v * (1 - s) + ssb * (5 * s - 1)) 
//     double s = params[0];
//     double R0 = params[1];
//     double v = params[2];
//     
//     // sratio is the recruits sex ratio, 1 if single sex model
//     //double sratio = 1.0;
//     //if (params.size() > 3) {
//     //  sratio = params[3];
//     //}
//     
//     // seasp is the prop of rec for the season, 1 if single rec
//     // TODO BUT srp needs to come from spwn season
//     //double seasp = 1.0;
//     //if (params.size() > 4) {
//     //  seasp = params[4];
//     //}
// 
//     // TODO ssbp is the prop of ssb active in a season, 1 if single rec
//     //double ssbp = 1.0;
//     //if (params.size() > 5) {
//     //  ssbp = params[5];
//     //}
//   
//     rec = (4.0 * s * R0 * srp) / (v * (1.0 - s) + srp * (5 * s - 1.0));
// 
//     return rec;
// }
// 
// template <typename T>
// T cushing(const T srp, const std::vector<double> params, char *f_name[]){
//   T rec;
//   // rec = a * srp ^ b
//   rec = params[0] * exp(log(srp)*params[1]);
//   return rec;
// }
// 
// template <typename T>
// T segreg(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     // rec = if(ssb < b) a * ssb else a * b
//     if(srp <= params[1]) {
//       rec = params[0] * srp;
//     }
//     else {
//       rec = params[0] * params[1];
//     }
//     return rec;
// }
// 
// template <typename T>
// T segregDa(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     // rec = if(ssb < b) a * ssb else a * b
//     if(pow(srp, params[2]) <= params[1]) {
//       rec = params[0] * pow(srp, params[2]);
//     }
//     else {
//       rec = params[0] * params[1];
//     }
//     return rec;
// }
// 
// 
// 
// template <typename T>
// T survsrr(const T ssf, const std::vector<double> params, char *f_name[]){
//     T rec;
//     
//     double R0 = params[0];
//     double sfrac = params[1];
//     double beta = params[2];
//     double SB0 = params[3];
// 
//     // sratio is the recruits sex ratio, default 0.5 assumes 2 sex model
//     double sratio = 0.5;
//     if (params.size() > 4) {
//       sratio = params[4];
//     }
// 
//     double z0 = log(1.0 / (SB0 / R0));
//     
//     double zmax = z0 + sfrac * (0.0 - z0);
//     
//     T zsurv = exp((1.0 - pow((ssf / SB0), beta)) * (zmax - z0) + z0);
// 
//     // Sex ratio at recruitment set at 1:1
//     rec = ssf * zsurv * sratio;
// 
//     return rec;
// }
// 
// template <typename T>
// T bevholtsig(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     // rec = a / ((b / srp) ^c + 1)
//     rec = params[0] / (pow((params[1] / srp), params[2]) + 1);
//     return rec;
// }
// 
// template <typename T>
// T mixedsrr(const T srp, const std::vector<double> params, char *f_name[]){
//     T rec;
//     
//     double a = params[0];
//     double b = params[1];
//     int mod = params[2];
// 
//     // 1 Bevholt, rec = a * srp / (b + srp)
//     if (mod == 1) {
//       rec = a * srp / (b + srp);
//     } else if (mod == 2) {
//     // 2 Ricker, rec = a * srp * exp (-b * srp)
//       rec = a * srp * exp(-b * srp);
//     // 3 Segreg, rec = if(ssb < b) a * ssb else a * b
//     } else if (mod == 3) {
//       if(srp <= b) {
//         rec = a * srp;
//       } else {
//         rec = a * b;
//       }
//     }
//     return rec;
// }
// 
// // Instantiate functions
template double customSRR(const double srp, const std::vector<double> params,const std::string model_name);
template adouble customSRR(const adouble srp, const std::vector<double> params,const std::string model_name);
// template double customSRR(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble customSRR(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double bevholt(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble bevholt(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double constant(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble constant(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double bevholtSS3(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble bevholtSS3(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double cushing(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble cushing(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double segreg(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble segreg(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double survsrr(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble survsrr(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double bevholtsig(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble bevholtsig(const adouble srp, const std::vector<double> params, char *f_name[]);
// template double mixedsrr(const double srp, const std::vector<double> params, char *f_name[]);
// template adouble mixedsrr(const adouble srp, const std::vector<double> params, char *f_name[]);
