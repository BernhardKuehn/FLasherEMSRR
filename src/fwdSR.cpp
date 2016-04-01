/* 
 * Copyright 2014 FLR Team. Distributed under the GPL 2 or later
 * Maintainer: Finlay Scott, JRC
 */

#include "../inst/include/fwdSR.h"

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
    map_model_name_to_function["ricker"] = &ricker;
    map_model_name_to_function["Ricker"] = &ricker;
    map_model_name_to_function["bevholt"] = &bevholt;
    map_model_name_to_function["Bevholt"] = &bevholt;
    map_model_name_to_function["constant"] = &constant;
    map_model_name_to_function["Constant"] = &constant;
    map_model_name_to_function["mean"] = &constant;
    map_model_name_to_function["Mean"] = &constant;
    map_model_name_to_function["geomean"] = &constant;
    map_model_name_to_function["Geomean"] = &constant;
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
 * Assumes all dimensions (e.g. of residuals and parameters) have been checked in R (no checks made here).
 *
 * \param model_name The name of the SR function (must be present in the model map else it fails).
 * \param params_ip The SR parameters. The parameters are stored in the first dimension and can be disaggregated by time, area etc.
 * \param residuals_ip Residuals that can be applied to the predicted recruitment.
 * \param residuals_mult_ip Are the residuals multiplicative (true) or additive (false).
 */
template <typename T>
fwdSR_base<T>::fwdSR_base(const std::string model_name, const FLQuant params_ip, const FLQuant residuals_ip, const bool residuals_mult_ip) {
    params = params_ip;
    residuals = residuals_ip;
    residuals_mult = residuals_mult_ip;
    init_model_map();

    // Set the model pointer
    typename model_map_type::const_iterator model_pair_found = map_model_name_to_function.find(model_name);
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
                            Rcpp::Named("residuals",residuals),
                            Rcpp::Named("residuals_mult",residuals_mult));
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
    params = fwdSR_source.params;
    residuals = fwdSR_source.residuals;
    residuals_mult = fwdSR_source.residuals_mult;
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
        params = fwdSR_source.params;
        residuals = fwdSR_source.residuals;
        residuals_mult = fwdSR_source.residuals_mult;
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
std::vector<double> fwdSR_base<T>::get_params(int year, int unit, int season, int area, int iter) const{
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
T fwdSR_base<T>::eval_model(const T srp, int year, int unit, int season, int area, int iter) const{
    // Get the parameters
    std::vector<double> model_params = get_params(year, unit, season, area, iter);
    // Finally, evaluate the function being pointed at
    T rec = model(srp, model_params);
    return rec;
}
/*!
 * \brief Evaluates the SR model
 * \param srp The spawning reproductive potential that produces the recruitment.
 * \param params_indices The indices of the SR params (starting at 1).
 */
template <typename T>
T fwdSR_base<T>::eval_model(const T srp, const std::vector<unsigned int> params_indices) const{ 
    // Check length of params_indices
    if (params_indices.size() != 5){
        Rcpp::stop("In fwdSR::eval_model. params_indices must be of length 5.");
    }
    T rec = eval_model(srp, params_indices[0], params_indices[1], params_indices[2], params_indices[3], params_indices[4]);
    return rec;
}
//@}

/*! \brief Predict recruitment
 *
 * Calculates the recruitment from an FLQuant of SRP, including the application of residuals.
 * The SRP can be a subset of the 'full' model SRP
 * (e.g. can be only one season out of all seasons, or several years out of all years).
 * It is therefore necessary to also pass in a vector of indices to specify the start
 * position of the SR params and residuals because we don't know the start position of the SRP argument relative to the whole operating model. 
 * i.e. we need to know the index of the params that corresponds with the first value in the SRP vector.
 * This is also true of the residuals.
 * Internally, the method calls the eval() method to calculate deterministic recruitment then applies the residuals.
 *
 * \param srp The spawning reproductive potential that produces the recruitment.
 * \param initial_params_indices A vector of length 5 (year, unit, ... iter) to specify the start position of the indices of the SR params relative to the 'whole' operating model (starting at 1).
 */
template <typename T>
FLQuant_base<T> fwdSR_base<T>::predict_recruitment(const FLQuant_base<T> srp, const std::vector<unsigned int> initial_params_indices, const std::vector<unsigned int> initial_residuals_indices){ 
    Rprintf("In predict_recruitment\n");
    if ((initial_params_indices.size() != 5) | (initial_residuals_indices.size() != 5)){
        Rcpp::stop("In fwdSR::predict_recruitment. initial_params_indices and initial_residuals_indices must be of length 5.\n");
    }
    std::vector<unsigned int> srp_dim = srp.get_dim();
    if (srp_dim[0] != 1){
        Rcpp::stop("In fwdSR::predict_recruitment. srp must be of length 1 in the first dimension.\n");
    }
    // Empty output object
    FLQuant_base<T> rec = srp;
    rec.fill(0.0);
    // Going to have to loop over the dimensions and update the params and residuals indices - not nice
    std::vector<unsigned int> params_indices = initial_params_indices;
    std::vector<unsigned int> residuals_indices = initial_residuals_indices;
    for (unsigned int year_counter = 1; year_counter <= srp_dim[1]; ++year_counter){
        params_indices[0] = initial_params_indices[0] + year_counter - 1;
        residuals_indices[0] = initial_residuals_indices[0] + year_counter - 1;
        for (unsigned int unit_counter = 1; unit_counter <= srp_dim[2]; ++unit_counter){
            params_indices[1] = initial_params_indices[1] + unit_counter - 1;
            residuals_indices[1] = initial_residuals_indices[1] + unit_counter - 1;
            for (unsigned int season_counter = 1; season_counter <= srp_dim[3]; ++season_counter){
                params_indices[2] = initial_params_indices[2] + season_counter - 1;
                residuals_indices[2] = initial_residuals_indices[2] + season_counter - 1;
                for (unsigned int area_counter = 1; area_counter <= srp_dim[4]; ++area_counter){
                    params_indices[3] = initial_params_indices[3] + area_counter - 1;
                    residuals_indices[3] = initial_residuals_indices[3] + area_counter - 1;
                    for (unsigned int iter_counter = 1; iter_counter <= srp_dim[5]; ++iter_counter){
                        params_indices[4] = initial_params_indices[4] + iter_counter - 1;
                        residuals_indices[4] = initial_residuals_indices[4] + iter_counter - 1;
                        T rec_temp = eval_model(srp(1, year_counter, unit_counter, season_counter, area_counter, iter_counter), params_indices);
                        if (residuals_mult == true){
                            rec_temp *= residuals(1, residuals_indices[0], residuals_indices[1], residuals_indices[2], residuals_indices[3], residuals_indices[4]);
                        }
                        else {
                            rec_temp += residuals(1, residuals_indices[0], residuals_indices[1], residuals_indices[2], residuals_indices[3], residuals_indices[4]);
                        }
                        rec(1, year_counter, unit_counter, season_counter, area_counter, iter_counter) = rec_temp;
                    }}}}}
    return rec;
}

/*! \brief Number of SR parameters.
 *
 * The number of SR parameters is the length of the first dimension of the parameters FLQuant.
 */
template <typename T>
int fwdSR_base<T>::get_nparams() const{
    return params.get_nquant();
}

/*! \brief Get the residuals.
 *
 * Returns the residuals.
 */
template <typename T>
FLQuant_base<double> fwdSR_base<T>::get_residuals() const{
    return residuals;
}

/*! \brief Get the residuals multiplier.
 *
 * Returns the residuals multiplier.
 */
template <typename T>
bool fwdSR_base<T>::get_residuals_mult() const{
    return residuals_mult;
}

// Explicit instantiation of class
template class fwdSR_base<double>;
template class fwdSR_base<adouble>;

//--------------------------------------------------------------------
// SRR functions
// The params in these functions do not have any disaggregation (e.g. by time or area)
// They are only the params required to evaluate the SRR function at a specific point
// The disaggregated parameter values are stored in the fwdSR_base class as an FLQuant_base object
// The fwdSR.eval_model method sorts out time step of params etc and passes the correct set of params to these functions
// These functions must have the same argument list so that it matches the typedef for the model pointer
template <typename T>
T ricker(const T srp, const std::vector<double> params){
    T rec;
    // rec = a * srp * exp (-b * srp)
    rec = params[0] * srp * exp(-params[1] * srp);
    return rec;
}

template <typename T>
T bevholt(const T srp, const std::vector<double> params){
    T rec;
    // rec = a * srp / (b + srp)
    rec = params[0] * srp / (params[1] + srp);
    return rec;
}

template <typename T>
T constant(const T srp, const std::vector<double> params){
    T rec;
    // rec = a 
    rec = params[0]; 
    return rec;
}

// Instantiate functions
template double ricker(const double ssb, const std::vector<double> params);
template adouble ricker(const adouble ssb, const std::vector<double> params);
template double bevholt(const double ssb, const std::vector<double> params);
template adouble bevholt(const adouble ssb, const std::vector<double> params);
template double constant(const double srp, const std::vector<double> params);
template adouble constant(const adouble srp, const std::vector<double> params);

