/* 
 * Copyright 2014 FLR Team. Distributed under the GPL 2 or later
 * Maintainer: Finlay Scott, JRC
 */

// For timing functions
#include <time.h>

#include "../inst/include/operating_model.h"

// Converting timestep to year and season and vice versa
// These are based on INDICES, not characters
template <typename T>
void year_season_to_timestep(const unsigned int year, const unsigned int season, const FLQuant_base<T>& flq, unsigned int& timestep){
    year_season_to_timestep(year, season, flq.get_nseason(), timestep);
}

template <typename T>
void timestep_to_year_season(const unsigned int timestep, const FLQuant_base<T>& flq, unsigned int& year, unsigned int& season){
    timestep_to_year_season(timestep, flq.get_nseason(), year, season);
}

void year_season_to_timestep(const unsigned int year, const unsigned int season, const unsigned int nseason, unsigned int& timestep){
    timestep = (year-1) * nseason + season;
}

void timestep_to_year_season(const unsigned int timestep, const unsigned int nseason, unsigned int& year, unsigned int& season){
    year =  (timestep-1) / nseason + 1; // integer divide - takes the floor
    season = (timestep-1) % nseason + 1;
}

// Instantiate
template void year_season_to_timestep(const unsigned int year, const unsigned int season, const FLQuant_base<double>& flq, unsigned int& timestep);
template void year_season_to_timestep(const unsigned int year, const unsigned int season, const FLQuant_base<adouble>& flq, unsigned int& timestep);

template void timestep_to_year_season(const unsigned int timestep, const FLQuant_base<double>& flq, unsigned int& year, unsigned int& season);
template void timestep_to_year_season(const unsigned int timestep, const FLQuant_base<adouble>& flq, unsigned int& year, unsigned int& season);

/*------------------------------------------------------------*/
// operatingModel class

/*! \brief Empty constructor that creates empty members
 */
operatingModel::operatingModel(){
    biols = fwdBiolsAD();
    fisheries = FLFisheriesAD();
    ctrl = fwdControl();
}

/*! \brief Main constructor for the operatingModel
 *
 * Checks dims of the arguments and everything is OK, creates an operatingModel.
 * The dimensions of the constituent member's FLQuants must all match along dimensions 1 - 5 else it throws and error.
 *
 * \param fisheries_in The fisheries.
 * \param biols_in The biological stocks.
 * \param ctrl_in The control object that controls the projections.
 */
operatingModel::operatingModel(const FLFisheriesAD fisheries_in, const fwdBiolsAD biols_in, const fwdControl ctrl_in){
    // Checking dims (1 - 5) of landings slots and biols are the same - but only for catches catching that biol?
    // Dim 6 must be 1 or n
    //Rcpp::IntegerVector catch_dim;
    //Rcpp::IntegerVector biol_dim;
    //const unsigned int nbiols = biols_in.get_nbiols();
    //const unsigned int nfisheries = fisheries_in.get_nfisheries();
    //const unsigned int ncatches = 0;
    //for (int unsigned biol_counter = 1; biol_counter <= nbiols; ++biol_counter){
    //    biol_dim = biols_in(biol_counter).n().get_dim();
    //    for (int unsigned fishery_counter = 1; fishery_counter <= nfisheries; ++fishery_counter){
    //        for (int unsigned catch_counter = 1; catch_counter <= ncatches; ++ catch_counter){
    //            catch_dim = fisheries_in(fishery_counter)(catch_counter).landings_n().get_dim(); 
    //            for (int dim_counter = 0; dim_counter < 5; ++dim_counter){
    //                if((biol_dim[dim_counter] != catch_dim[dim_counter])){
    //                    Rcpp::stop("In operatingModel constructor: Biol dims must be the same as Catch dims\n");
    //                }
    //            }
    //        }
    //    }
    //}
    // Check residuals of the SRR - need to have them for the projection years in the control object
    
    // Add ITER check for ctrl
    biols = biols_in;
    fisheries = fisheries_in;
    ctrl = ctrl_in;

    // Again add iteration check here
    // Iters in ctrl need to match those in biols
    // Shouldn't it recycle?
    if (ctrl.get_niter() != biols(1).n().get_niter()){
        Rcpp::stop("In operatingModel constructor. Iters in biol must equal those in fwdControl.\n");
    }
}

/*! \brief The copy constructor
 *
 * Ensures a deep copy and prevents members being pointed at by multiple instances
 * \param operatingModel_source The object to be copied.
 */
operatingModel::operatingModel(const operatingModel& operatingModel_source){
    biols = operatingModel_source.biols;
    fisheries = operatingModel_source.fisheries;
    ctrl = operatingModel_source.ctrl;
}

/*! \brief The assignment operator 
 *
 * Ensures a deep copy and prevents members being pointed at by multiple instances
 * \param operatingModel_source The object to be copied.
 */
operatingModel& operatingModel::operator = (const operatingModel& operatingModel_source){
	if (this != &operatingModel_source){
        biols = operatingModel_source.biols;
        fisheries = operatingModel_source.fisheries;
        ctrl = operatingModel_source.ctrl;
	}
	return *this;
}

/*! \brief Intrusive wrap
 *
 * Allows a direct return of an operatingModel to R.
 * It returns an R list of components.
 */
operatingModel::operator SEXP() const{
    //Rprintf("Wrapping operatingModel.\n");
    return Rcpp::List::create(
                            Rcpp::Named("biols", biols),
                            Rcpp::Named("fisheries", fisheries),
                            Rcpp::Named("ctrl", ctrl));
}

/*! \brief Calculates the spawning recruitment potential of a biol
 *
 * Calculated as SSB: N*mat*wt*exp(-Fprespwn - m*spwn) 
 * where the natural mortality m is assumed to be constant over the timestep
 * and Fprespwn is how much fishing mortality happened before spawning
 * \param biol_no The position of the biol in the biols list (starting at 1)
 * \param indices_min The minimum indices: year, unit etc (length 5)
 * \param indices_max The maximum indices: year, unit etc (length 5)
 */
FLQuantAD operatingModel::srp(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const{
    // Check indices_min and indices_max are of length 5
    if (indices_min.size() != 5 | indices_max.size() != 5){
        Rcpp::stop("In operatingModel srp subsetter. Indices not of length 5\n");
    }
    // Add age range to indices
    std::vector<unsigned int> dim = biols(biol_no).n().get_dim();
    std::vector<unsigned int> new_indices_min = indices_min;
    std::vector<unsigned int> new_indices_max = indices_max;
    new_indices_min.insert(new_indices_min.begin(), 1);
    new_indices_max.insert(new_indices_max.begin(), dim[0]);

    // Get Fprespawn = F * Fpropspawn for each F fishing the biol
    // 1. ID which F are fishing that B (timing is determined by ftime and spwn)
    // 2. Get F(F,C,B)
    // 3. Get Fpropspawn(F,C,B)
    // 4. sum (F * Fpropspawn)
    

    FLQuantAD f_pre_spwn;

    // Can we get a pointer / iterator to the biol and dereference it instead of typing biols(biol_no) all the time? Why? Faster?
    FLQuantAD srp = quant_sum(biols(biol_no).n(new_indices_min, new_indices_max) *
        biols(biol_no).wt(new_indices_min, new_indices_max) *
        biols(biol_no).mat(new_indices_min, new_indices_max) *
        // Add in f_pre_spwn here
        exp(-1.0 * biols(biol_no).m(new_indices_min, new_indices_max) * biols(biol_no).spwn(new_indices_min, new_indices_max)));

    return srp;
}

/*! \brief Calculates the proportion of fishing mortality that occurs before the stock spwns
 *
 * Given by the start and stop time of fishing (ftime member of FLFishery) and time of spawning (spwn member of fwdBiol)
 * \param fishery_no The position of the fishery in the fisheries list (starting at 1)
 * \param biol_no The position of the biol in the biols list (starting at 1)
 * \param indices_min The minimum indices: year, unit etc (length 5)
 * \param indices_max The maximum indices: year, unit etc (length 5)
 */
FLQuant operatingModel::f_prop_spwn(const int fishery_no, const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const{
    // Make an object to dump the result into
    FLQuant propf_out(1, indices_max[0]-indices_min[0]+1, indices_max[1]-indices_min[1]+1, indices_max[2]-indices_min[2]+1, indices_max[3]-indices_min[3]+1, indices_max[4]-indices_min[4]+1);
    // Need to calculate element by element as timing can change over years etc.
    double propf = 0.0;
    for (unsigned int year_count=indices_min[0]; year_count <= indices_max[0]; ++year_count){
        for (unsigned int unit_count=indices_min[1]; unit_count <= indices_max[1]; ++unit_count){
            for (unsigned int season_count=indices_min[2]; season_count <= indices_max[2]; ++season_count){
                for (unsigned int area_count=indices_min[3]; area_count <= indices_max[3]; ++area_count){
                    for (unsigned int iter_count=indices_min[4]; iter_count <= indices_max[4]; ++iter_count){
                        // Fix this depending on representation of fperiod
                        double fstart = fisheries(fishery_no).ftime()(1,year_count, unit_count, season_count, area_count, iter_count);
                        double fend = fisheries(fishery_no).ftime()(2,year_count, unit_count, season_count, area_count, iter_count);
                        double spwn = biols(biol_no).spwn()(1,year_count, unit_count, season_count, area_count, iter_count);
                        if (fend < spwn){
                            propf = 1.0;
                        }
                        else if (fstart > spwn){
                            propf = 0.0;
                        }
                        else {
                            propf = (spwn - fstart) / (fend - fstart);
                        }
                        // Dump it in the right place - ugliness abounds
                        propf_out(1, year_count - indices_min[0], unit_count - indices_min[1], season_count - indices_min[2], area_count - indices_min[3], iter_count - indices_min[4]); 
    }}}}}
    return propf_out;
}



///*!
// * \name Catchability
// * Calculate the catchability of a fishery / catch fishing a biol. 
// * Catchability = alpha * Biomass ^ -beta,
// * where alpha and beta are the catchability params.
// * It is assumed that the fishery / catch catches the biol (no check is made).
// */
////@{
///*!
// * \brief Catchability of a timestep (subset over year, unit, season, area, iter)
// * \param fishery_no the position of the fishery within the fisheries (starting at 1).
// * \param catch_no the position of the catch within the fishery (starting at 1).
// * \param biol_no the position of the biol within the biols (starting at 1).
// * \param indices_min The minimum indices year, unit etc (length 5)
// * \param indices_max The maximum indices year, unit etc (length 5)
// */
//FLQuantAD operatingModel::catch_q(const int fishery_no, const int catch_no, const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const{
//    if (indices_min.size() != 5 | indices_max.size() != 5){
//        Rcpp::stop("In operatingModel catch_q subsetter. Indices not of length 5\n");
//    }
//    // Make the empty FLQuant with one year and season
//    std::vector<unsigned int> new_dim {indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1};
//    FLQuantAD q(1, new_dim[0], new_dim[1], new_dim[2], new_dim[3], new_dim[4]);
//    FLQuantAD biomass = biols(biol_no).biomass(indices_min, indices_max);
//    // Need to add an extra value to the start of the indices when subsetting q_params
//    std::vector<unsigned int> indices_min_q = indices_min;
//    std::vector<unsigned int> indices_max_q = indices_max;
//    indices_min_q.insert(indices_min_q.begin(), 1);
//    indices_max_q.insert(indices_max_q.begin(), 1);
//    FLQuant q_params_a = fisheries(fishery_no, catch_no).catch_q_params(indices_min_q, indices_max_q);
//    indices_min_q[0] = 2;
//    indices_max_q[0] = 2;
//    FLQuant q_params_b = fisheries(fishery_no, catch_no).catch_q_params(indices_min_q, indices_max_q);
//    // The magic of lambda functions
//    std::transform(biomass.begin(), biomass.end(), q_params_b.begin(), biomass.begin(), [](adouble x, double y){return pow(x, -y);});
//    q = q_params_a * biomass;
//    return q;
//}
//
///*!
// * \brief Single value catchability for a specific quant, year, iter etc.
// * \param fishery_no the position of the fishery within the fisheries (starting at 1).
// * \param catch_no the position of the catch within the fishery (starting at 1).
// * \param biol_no the position of the biol within the biols (starting at 1).
// * \param year (starting at 1)
// * \param unit (starting at 1)
// * \param season (starting at 1)
// * \param area (starting at 1)
// * \param iter (starting at 1)
// */
//adouble operatingModel::catch_q(const int fishery_no, const int catch_no, const int biol_no, const unsigned int year, const unsigned int unit, const unsigned int season, const unsigned int area, const unsigned int iter) const{
//    // Call subsetter
//    std::vector<unsigned int> indices_min {year, unit, season, area, iter};
//    std::vector<unsigned int> indices_max {year, unit, season, area, iter};
//    FLQuantAD q = catch_q(fishery_no, catch_no, biol_no, indices_min, indices_max);
//    return q(1,1,1,1,1,1);
//}
//
///*!
// * \brief Catchability over all dimensions
// * \param fishery_no the position of the fishery within the fisheries (starting at 1).
// * \param catch_no the position of the catch within the fishery (starting at 1).
// * \param biol_no the position of the biol within the biols (starting at 1).
// */
//FLQuantAD operatingModel::catch_q(const int fishery_no, const int catch_no, const int biol_no) const{
//    FLQuantAD biomass = biols(biol_no).biomass();
//    std::vector<unsigned int> dim = biomass.get_dim();
//    std::vector<unsigned int> indices_min {1,1,1,1,1};
//    std::vector<unsigned int> indices_max {dim[1], dim[2], dim[3], dim[4], dim[5]};
//    FLQuantAD q = catch_q(fishery_no, catch_no, biol_no, indices_min, indices_max);
//    return q;
//}
////@}

/*! \name get_f
 * Calculate the instantaneous fishing mortality 
 * This method is the workhorse fishing mortality method that is called by other fishing mortality methods that do make checks.
 * F = effort * selectivity * catchability.
 * F = effort * selectivity * alpha * biomass ^ -beta
 */
//@{
/*! \brief Calculate the instantaneous fishing mortality of a single biol from a single fishery / catch over a subset of dimensions.
 * It is assumed that the fishery / catch actually fishes the biol (no check is made).
 * \param fishery_no the position of the fishery within the fisheries (starting at 1).
 * \param catch_no the position of the catch within the fishery (starting at 1).
 * \param biol_no the position of the biol within the biols (starting at 1).
 * \param indices_min The minimum indices quant, year, unit etc (length 6)
 * \param indices_max The maximum indices quant, year, unit etc (length 6)
*/
FLQuantAD operatingModel::get_f(const int fishery_no, const int catch_no, const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
    if (indices_min.size() != 6 | indices_max.size() != 6){
        Rcpp::stop("In operatingModel get_f subsetter. Indices not of length 6\n");
    }
    // Lop off the first value from the indices to get indices without quant - needed for effort and catch_q
    std::vector<unsigned int> indices_min5(indices_min.begin()+1, indices_min.end());
    std::vector<unsigned int> indices_max5(indices_max.begin()+1, indices_max.end());
    FLQuantAD effort = fisheries(fishery_no).effort(indices_min5, indices_max5);
    FLQuantAD sel = fisheries(fishery_no, catch_no).catch_sel()(indices_min, indices_max);
    FLQuantAD biomass = biols(biol_no).biomass(indices_min5, indices_max5);
    FLQuantAD fout(indices_max[0]-indices_min[0]+1, indices_max[1]-indices_min[1]+1, indices_max[2]-indices_min[2]+1, indices_max[3]-indices_min[3]+1, indices_max[4]-indices_min[4]+1, indices_max[5] - indices_min[5]+1);
    // Ugly code - q_params, effort and sel are all different dims
    std::vector<unsigned int> dim = fout.get_dim();
    for (int iter_count = 1; iter_count <= dim[5]; ++iter_count){
        for (int area_count = 1; area_count <= dim[4]; ++area_count){
            for (int season_count = 1; season_count <= dim[3]; ++season_count){
                for (int unit_count = 1; unit_count <= dim[2]; ++unit_count){
                    for (int year_count = 1; year_count <= dim[1]; ++year_count){
                        // Get alpha * biomass ^ -beta * effort (not age structured)
                        std::vector<double> q_params = fisheries(fishery_no, catch_no).catch_q_params(year_count + indices_min[1] - 1, unit_count + indices_min[2] - 1, season_count + indices_min[3] - 1, area_count + indices_min[4] - 1, iter_count + indices_min[5] - 1);
                        adouble temp_qe = q_params[0] * pow(biomass(1, year_count, unit_count, season_count, area_count, iter_count), -1.0 * q_params[1]) * effort(1, year_count, unit_count, season_count, area_count, iter_count);
                        for (int quant_count = 1; quant_count <= dim[0]; ++quant_count){
                            // * sel (which is age structured)
                            fout(quant_count, year_count, unit_count, season_count, area_count, iter_count) =  temp_qe * sel(quant_count, year_count, unit_count, season_count, area_count, iter_count);
    }}}}}}
    return fout;
}

/*! \brief Calculate the instantaneous fishing mortality of a single biol from a single fishery / catch over all dimensions.
 * It is assumed that the fishery / catch actually fishes the biol (no check is made).
 * \param fishery_no the position of the fishery within the fisheries (starting at 1).
 * \param catch_no the position of the catch within the fishery (starting at 1).
 * \param biol_no the position of the biol within the biols (starting at 1).
 */
FLQuantAD operatingModel::get_f(const int fishery_no, const int catch_no, const int biol_no) const {
    // Just call the subset method with full indices
    std::vector<unsigned int> indices_max = biols(biol_no).n().get_dim();
    std::vector<unsigned int> indices_min(6,1);
    FLQuantAD f = get_f(fishery_no, catch_no, biol_no, indices_min, indices_max);
    return f;
}

//@}

/*! \brief Total instantaneous fishing mortality on a biol over a subset of dimensions
 * \param biol_no the position of the biol within the biols (starting at 1).
 * \param indices_min minimum indices for subsetting (quant - iter, vector of length 6)
 * \param indices_max maximum indices for subsetting (quant - iter, vector of length 6)
 */
FLQuantAD operatingModel::get_f(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const{
    unsigned int fishery_no;
    unsigned int catch_no;
    // We need to know the Fishery / Catches that catch the biol
    const Rcpp::IntegerMatrix FC =  ctrl.get_FC(biol_no);
    // What happens if no-one is fishing that biol? FC.nrow() == 0 so loop never triggers
    FLQuantAD total_f(indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1, indices_max[5] - indices_min[5] + 1); 
    total_f.fill(0.0);
    for (unsigned int f_counter=0; f_counter < FC.nrow(); ++f_counter){
        total_f = total_f + get_f(FC(f_counter,0), FC(f_counter,1), biol_no, indices_min, indices_max);
    }
    return total_f;
}

/*! \brief Total instantaneous fishing mortality on a biol over all dimensions
 *
 * \param biol_no the position of the biol within the biols (starting at 1).
 */
FLQuantAD operatingModel::get_f(const int biol_no) const {
    std::vector<unsigned int> indices_max = biols(biol_no).n().get_dim();
    std::vector<unsigned int> indices_min(6,1);
    FLQuantAD f = get_f(biol_no, indices_min, indices_max);
    return f;
}
//@}

///*! \name Calculate the partial instantaneous fishing mortality on a single biol from a single fishery / catch
// * Checks are made to see if the fishery / catch actually catches the biol.
// * If the fishery / catch does not actually fish that biol, then the partial fishing mortality will be 0.
// * biol_no is included as an argument in case the fishery / catch catches more than one biol.
// */
////@{
///*! \brief Instantaneous partial fishing mortality over a subset of dimensions
// * \param fishery_no the position of the fishery within the fisheries (starting at 1).
// * \param catch_no the position of the catch within the fishery (starting at 1).
// * \param biol_no the position of the biol within the biols (starting at 1).
// * indices_min minimum subset indices (length 6)
// * indices_max maximum subset indices (length 6)
// */
//FLQuantAD operatingModel::partial_f(const int fishery_no, const int catch_no, const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    if (indices_min.size() != 6 | indices_max.size() != 6){
//        Rcpp::stop("In operatingModel partial_f subsetter. Indices not of length 6\n");
//    }
//    FLQuantAD partial_f;
//    std::vector<int> B = ctrl.get_B(fishery_no, catch_no);
//    // Do F / C catch anything? If not, something has probably gone wrong but we'll leave it at the moment.
//    if (B.size() == 0){
//        partial_f = biols(biol_no).n(); // Just get FLQuant of the correct dims
//        partial_f.fill(0.0);
//    }
//    std::vector<int>::iterator B_iterator = find(B.begin(), B.end(), biol_no);
//    if(B_iterator != B.end()){
//        partial_f = get_f(fishery_no, catch_no, biol_no, indices_min, indices_max);
//    }
//    // biol_no is not caught by F/C - return 0
//    else {
//        partial_f = biols(biol_no).n(); // Just get FLQuant of the correct dims
//        partial_f.fill(0.0);
//    }
//    return partial_f;
//}
///*! \brief Instantaneous partial fishing mortality over all dimensions
// * \param fishery_no the position of the fishery within the fisheries (starting at 1).
// * \param catch_no the position of the catch within the fishery (starting at 1).
// * \param biol_no the position of the biol within the biols (starting at 1).
// */
//FLQuantAD operatingModel::partial_f(const int fishery_no, const int catch_no, const int biol_no) const {
//    // Just call the subset method with full indices
//    //Rprintf("In partial_f FLQ method\n");
//    //Rcpp::IntegerVector raw_dims = biols(biol_no).n().get_dim();
//    //std::vector<unsigned int> indices_max = Rcpp::as<std::vector<unsigned int>>(raw_dims);
//    std::vector<unsigned int> indices_max = biols(biol_no).n().get_dim();
//    std::vector<unsigned int> indices_min(6,1);
//    FLQuantAD f = partial_f(fishery_no, catch_no, biol_no, indices_min, indices_max);
//    return f;
//}
////@}



////! Calculate the total mortality on a biol
///*!
//    Calculate the total instantaneous mortality on a single biol.
//    This is the sum of the fishing mortality from all fishery / catches that fish it and the natural mortality.
//    The total mortality is calculated over all dimensions (quant, year, etc. ).
//    \param biol_no the position of the biol within the biols (starting at 1).
// */
//FLQuantAD operatingModel::z(const int biol_no) const {
//    FLQuantAD z = total_f(biol_no) + biols(biol_no).m();
//    return z;
//}
//// The indices_min and max version
//FLQuantAD operatingModel::z(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    FLQuantAD z = total_f(biol_no, indices_min, indices_max) + biols(biol_no).m(indices_min, indices_max);
//    return z;
//}
//
//
//// The timestep that fmult affects to calculate the target value
//// e.g. if Biomass is target, then adjust the fmult in the previous timestep
//// Add more target types to it
////int operatingModel::get_target_effort_timestep(const int target_no){
////    fwdControlTargetType target_type = ctrl.get_target_type(target_no);
////    unsigned int target_year = ctrl.get_target_year(target_no);
////    unsigned int target_season = ctrl.get_target_season(target_no);
////    unsigned int target_timestep = 0;
////    year_season_to_timestep(target_year, target_season, biols(1).n(), target_timestep);
////    // Biomass timesteps
////    if((target_type == target_ssb) ||
////       (target_type == target_biomass)){
////        target_timestep = target_timestep - 1;
////    }
////    return target_timestep;
////}
//
////! Project the Fisheries in the operatingModel by a single timestep
///*!
//    Projects the Fisheries in the operatingModel by a single timestep.
//    All catches, landings and discards in the Fisheries are updated for that timestep based on effort in that timestep
//    The Baranov catch equation is used to calculate catches.
//    This assumes that the instantaneous rate of fishing and natural mortalities are constant over time and age and occur simultaneously.
//    \param timestep The time step for the projection.
// */
//void operatingModel::project_fisheries(const int timestep){
//    // C = (pF / Z) * (1 - exp(-Z)) * N
//    Rprintf("In operatingModel::project_fisheries\n");
//    // What timesteps / years / seasons are we dealing with?
//    unsigned int year = 0;
//    unsigned int season = 0;
//    std::vector<unsigned int> catch_dim = fisheries(1,1).landings_n().get_dim();
//    timestep_to_year_season(timestep, catch_dim[3], year, season);
//    // timestep checks
//    if ((year > catch_dim[1]) | (season > catch_dim[3])){
//        Rcpp::stop("In operatingModel::project_fisheries. timestep outside of range");
//    }
//
//    // CAREFUL WITH NUMBER OF ITERS - do we allow different iterations across objects?
//    // Maybe move this below into loop?
//    unsigned int niter = catch_dim[5];
//    // Not yet set up for units and areas
//    unsigned int unit = 1;
//    unsigned int area = 1;
//
//    // Get total_fs for all biols
//    // Possible different age ranges - hence the indices_min and max
//    FLQuant7AD total_fs;
//    for (int biol_counter=1; biol_counter <= biols.get_nbiols(); ++biol_counter){
//        std::vector<unsigned int> biol_dim = biols(1).n().get_dim();
//        std::vector<unsigned int> indices_min{1, year, unit, season, area, 1};
//        std::vector<unsigned int> indices_max{biol_dim[0], year, unit, season, area, niter};
//        total_fs(total_f(biol_counter, indices_min, indices_max));  
//    }
//    
//    std::vector<int> biols_fished;
//    unsigned int biol_no;
//    FLQuantAD landings;
//    FLQuantAD discards;
//    FLQuantAD discards_ratio_temp;
//    FLQuantAD z_temp;
//    for (int fishery_counter=1; fishery_counter <= fisheries.get_nfisheries(); ++fishery_counter){
//        for (int catch_counter=1; catch_counter <= fisheries(fishery_counter).get_ncatches(); ++catch_counter){
//            // Set up indices for subsetting the timestep
//            // Catch must have same dims as the biol it catches
//            std::vector<unsigned int> catch_dim = fisheries(fishery_counter, catch_counter).landings_n().get_dim();
//            std::vector<unsigned int> indices_min{1, year, unit, season, area, 1};
//            std::vector<unsigned int> indices_max{catch_dim[0], year, unit, season, area, niter};
//            FLQuantAD catches(indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1, indices_max[5] - indices_min[5] + 1); 
//            catches.fill(0.0);
//            // Get B indices for the F/C
//            biols_fished = ctrl.get_B(fishery_counter, catch_counter);
//            if (biols_fished.size() == 0){
//                Rcpp::stop("In project timestep, fishery / catch does not fish anything. Probably a setup error. Stopping\n");
//            }
//            // Get catches from each biol that is fished - in case a Catch fishes more than 1 Biol 
//            for (int biol_counter=1; biol_counter <= biols_fished.size(); ++biol_counter){
//                //Rprintf("biol being fished: %i\n", biols_fished[biol_counter-1]);
//                // C = (pF / Z) * (1 - exp(-Z)) * N
//                biol_no = biols_fished[biol_counter-1];
//                //Rprintf("fishery: %i; catch: %i; biol: %i\n", fishery_counter, catch_counter, biol_no);
//                //Inefficient! Calculating f again - both partial and total_f use the same method which is now being called multiple times
//                z_temp = total_fs(biol_no) + biols(biol_no).m(indices_min, indices_max);
//                // This line is expensive - extra calculation of partial_f
//                catches = catches +
//                    (partial_f(fishery_counter, catch_counter, biol_no, indices_min, indices_max) / z_temp) *
//                    (1.0 - exp(-1.0 * z_temp)) *
//                    biols(biol_no).n(indices_min, indices_max); 
//            } 
//            // Calc landings and discards from catches and discards.ratio
//            discards_ratio_temp = fisheries(fishery_counter, catch_counter).discards_ratio(indices_min, indices_max); 
//            landings = catches * (1.0 - discards_ratio_temp);
//            discards = catches * discards_ratio_temp;
//            // Dish out catches to landings and discards
//            // Only into the timestep - we don't want to overwrite other timesteps
//            for (unsigned int quant_count = 1; quant_count <= fisheries(fishery_counter, catch_counter).landings_n().get_nquant(); ++quant_count){
//                // Again - careful with iters
//                for (unsigned int iter_count = 1; iter_count <= fisheries(fishery_counter, catch_counter).landings_n().get_niter(); ++iter_count){
//                    fisheries(fishery_counter, catch_counter).landings_n()(quant_count,year,unit,season,area,iter_count) = landings(quant_count,1,unit,1,area,iter_count); 
//                    fisheries(fishery_counter, catch_counter).discards_n()(quant_count,year,unit,season,area,iter_count) = discards(quant_count,1,unit,1,area,iter_count); 
//                }
//            }
//        }
//    }
//    return;
//}
//
////! Project the Biols in the operatingModel by a single timestep
///*!
//    Projects the Biols in the operatingModel by a single timestep.
//    All abundances in the Biols are updated in the timestep based on the effort in the previous timestep.
//    This assumes that the instantaneous rate of fishing and natural mortalities are constant over time and age and occur simultaneously.
//    If a Biol is caught by multiple Catches, the fishing mortalities happen at the same time (and at the same time as the natural mortality) in the timestep.
//    \param timestep The time step for the projection.
// */
///*
//void operatingModel::project_biols(const int timestep){
//    // N2 = N1 * (exp -Z)
//    Rprintf("In operatingModel::project_biols\n");
//    if (timestep < 2){
//        Rcpp::stop("In operatingModel::project_biols. timestep must be at least 2.");
//    }
//
//    unsigned int year = 0;
//    unsigned int season = 0;
//    unsigned int prev_year = 0;
//    unsigned int prev_season = 0;
//    std::vector<unsigned int> biol_dim = biols(1).n().get_dim(); // All biols have same year and season range, so pick first one
//    timestep_to_year_season(timestep, biol_dim[3], year, season);
//    timestep_to_year_season(timestep-1, biol_dim[3], prev_year, prev_season);
//    // timestep checks
//    if ((year > biol_dim[1]) | (season > biol_dim[3])){
//        Rcpp::stop("In operatingModel::project_biols. Timestep outside of range");
//    }
//    // CAREFUL WITH NUMBER OF ITERS - do we allow different iterations across objects? - No
//    // In which case we will need to store niter at an operatingModel level
//    unsigned int niter = biol_dim[5];
//    // Not yet set up for units and areas
//    unsigned int unit = 1;
//    unsigned int area = 1;
//
//    // Update biol in timestep based on effort in previous timestep
//    unsigned int ssb_year = 0;
//    unsigned int ssb_season = 0;
//    // In what age do we put next timestep abundance? If beginning of year, everything is one year older
//    int age_shift = 0;
//    if (season == 1){
//        age_shift = 1;
//    }
//
//    for (unsigned int biol_counter=1; biol_counter <= biols.get_nbiols(); ++biol_counter){
//        // Indices for subsetting the previous timestep
//        std::vector<unsigned int> biol_dim = biols(biol_counter).n().get_dim();
//        std::vector<unsigned int> indices_min{1, prev_year, unit, prev_season, area, 1};
//        std::vector<unsigned int> indices_max{biol_dim[0], prev_year, unit, prev_season, area, niter};
//
//        // SSB and Recruitment
//        // Get the year and season of the SSB that will result in recruitment in the timestep
//        unsigned int ssb_timestep = timestep - biols(biol_counter).srr.get_timelag();
//        if (ssb_timestep < 1){
//            Rcpp::stop("project_timestep: ssb timestep outside of range");
//        }
//        timestep_to_year_season(ssb_timestep, biol_dim[3], ssb_year, ssb_season);
//        // Update min and max_indices for ssb and get it
//        std::vector<unsigned int> ssb_indices_min = {ssb_year, unit, ssb_season, area, 1};
//        std::vector<unsigned int> ssb_indices_max = {ssb_year, unit, ssb_season, area, niter};
//        FLQuantAD ssb_flq = ssb(biol_counter, ssb_indices_min, ssb_indices_max);
//
//        // Get abundances in next timestep from rec+1 to pg
//        FLQuantAD z_temp = z(biol_counter, indices_min, indices_max);
//        FLQuantAD next_n = biols(biol_counter).n(indices_min, indices_max) * exp(-1.0 * z_temp);
//
//        // Update biol in next timestep using next n = current n * exp(-z) taking care of age shifts and plus groups
//        for (int iter_count = 1; iter_count <= niter; ++iter_count){
//            // Recruitment
//            // rec is calculated in a scalar way - i.e. not passing it a vector of SSBs, so have to do one iter at a time
//            adouble ssb_temp = ssb_flq(1,1,1,1,1,iter_count);
//            //Rprintf("year: %i, season: %i, ssb_year: %i, ssb_season: %i, ssb_temp: %f\n", year, season, ssb_year, ssb_season, Value(ssb_temp));
//            adouble rec_temp = biols(biol_counter).srr.eval_model(ssb_temp, year, 1, season, 1, iter_count);
//            //Rprintf("iter_count %i; ssb_temp %f; rec_temp %f\n", iter_count, Value(ssb_temp), Value(rec_temp));
//            // Apply the residuals to rec_temp
//            // Use of if statement is OK because for each taping it will only branch the same way (depending on residuals_mult)
//            if (biols(biol_counter).srr.get_residuals_mult()){
//                //rec_temp = rec_temp * exp(biols(biol_counter).srr.get_residuals()(1,year,1,season,1,iter_count));
//                rec_temp = rec_temp * biols(biol_counter).srr.get_residuals()(1,year,1,season,1,iter_count);
//            }
//            else{
//                rec_temp = rec_temp + biols(biol_counter).srr.get_residuals()(1,year,1,season,1,iter_count);
//            }
//
//            // Update ages rec+1 upto final age 
//            for (int quant_count = 2; quant_count <= biol_dim[0]; ++quant_count){
//                biols(biol_counter).n()(quant_count, year, 1, season, 1, iter_count) = next_n(quant_count - age_shift, 1,unit,1,area,iter_count);
//            }
//            // Fix Plus Group and add in recruitment
//            // If next_season is start of the year assume last age is a plusgroup and abundance in first age group is ONLY recruitment
//            if (season == 1){
//                biols(biol_counter).n()(biol_dim[0], year, 1, season, 1, iter_count) = biols(biol_counter).n()(biol_dim[0], year, 1, season, 1, iter_count) + next_n(biol_dim[0],1,unit,1,area,iter_count);
//                // Rec in first age is only the calculated recruitment
//                biols(biol_counter).n()(1,year,1,season,1,iter_count) = rec_temp;
//            }
//            // Otherwise add in recruitment to first age group
//            else {
//                // Add rec to the existing age 1
//                biols(biol_counter).n()(1, year, 1, season, 1, iter_count) = next_n(1, 1, unit, 1, area, iter_count) + rec_temp;
//            }
//        }
//    }
//
//
//    return;
//}
//*/
//
//
//// Baranov catch equation: assumes that instantaneous rate of fishing and natural mortalities are constant over time and age
//// natural mortality and fishing mortality occur simultaneously.
//// If a biol is caught my multiple catches, the Fs happen at the same time (and at the same time as M) in the timestep
////! Project the operatingModel by a single timestep
///*!
//    Projects the operatingModel by a single timestep.
//    All catches, landings and discards in the Fisheries are updated for that timestep.
//    All abundnces in the Biols are updated for the following timestep.
//    The Baranov catch equation is used to calculate catches.
//    This assumes that the instantaneous rate of fishing and natural mortalities are constant over time and age and occur simultaneously.
//    If a Biol is caught by multiple Catches, the fishing mortalities happen at the same time (and at the same time as the natural mortality) in the timestep.
//    \param timestep The time step for the projection.
// */
///*
//void operatingModel::project_timestep(const int timestep){
//    // C = (pF / Z) * (1 - exp(-Z)) * N
//    // N2 = N1 * (exp -Z)
//    Rprintf("In project\n");
//
//    // Timing starts - CPU time
//    clock_t start, end, start2, end2;
//    start2 = clock();
//
//    // What timesteps / years / seasons are we dealing with?
//    unsigned int year = 0;
//    unsigned int season = 0;
//    unsigned int next_year = 0;
//    unsigned int next_season = 0;
//    timestep_to_year_season(timestep, biols(1).n(), year, season);
//    timestep_to_year_season(timestep+1, biols(1).n(), next_year, next_season); 
//
//    // Bit of faffing dims are not unsigned in Rcpp::IntegerVector
//    //Rcpp::IntegerVector biol_raw_dim = biols(1).n().get_dim();
//    //std::vector<unsigned int> biol_dim = Rcpp::as<std::vector<unsigned int>>(biol_raw_dim);
//    std::vector<unsigned int> biol_dim = biols(1).n().get_dim();
//    // timestep checks
//    if ((year > biol_dim[1]) | (season > biol_dim[3])){
//        Rcpp::stop("project_timestep: timestep outside of range");
//    }
//
//    // CAREFUL WITH NUMBER OF ITERS - do we allow different iterations across objects? - No
//    // In which case we will need to store niter at an operatingModel level
//    unsigned int niter = biol_dim[5];
//    // Not yet set up for units and areas
//    unsigned int unit = 1;
//    unsigned int area = 1;
//
//    // Set up indices for subsetting the timestep
//    std::vector<unsigned int> indices_min{1, year, unit, season, area, 1};
//    std::vector<unsigned int> indices_max{biol_dim[0], year, unit, season, area, niter};
//
//    start = clock();
//    // Get total_fs for all biols
//    FLQuant7AD total_fs;
//    for (int biol_counter=1; biol_counter <= biols.get_nbiols(); ++biol_counter){
//        total_fs(total_f(biol_counter, indices_min, indices_max));  
//    }
//
//    std::vector<int> biols_fished;
//    unsigned int biol_no;
//    // Make some empty objects for temp storage
//    std::vector<unsigned int> subset_dims = {indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1, indices_max[5] - indices_min[5] + 1}; 
//    FLQuantAD catches(subset_dims[0], subset_dims[1], subset_dims[2], subset_dims[3], subset_dims[4], subset_dims[5]); 
//    FLQuantAD landings;
//    FLQuantAD discards;
//    FLQuantAD discards_ratio_temp;
//    FLQuantAD z_temp;
//    for (int fishery_counter=1; fishery_counter <= fisheries.get_nfisheries(); ++fishery_counter){
//        for (int catch_counter=1; catch_counter <= fisheries(fishery_counter).get_ncatches(); ++catch_counter){
//            catches.fill(0.0);
//            // Get B indices for the F/C
//            biols_fished = ctrl.get_B(fishery_counter, catch_counter);
//            if (biols_fished.size() == 0){
//                Rcpp::stop("In project timestep, fishery / catch does not fish anything. Probably a setup error. Stopping\n");
//            }
//            // Get catches from each biol that is fished - in case a Catch fishes more than 1 Biol 
//            for (int biol_counter=1; biol_counter <= biols_fished.size(); ++biol_counter){
//                //Rprintf("biol being fished: %i\n", biols_fished[biol_counter-1]);
//                // C = (pF / Z) * (1 - exp(-Z)) * N
//                biol_no = biols_fished[biol_counter-1];
//                //Rprintf("fishery: %i; catch: %i; biol: %i\n", fishery_counter, catch_counter, biol_no);
//                //In efficient! Calculating f again - both partial and total_f use the same method
//                z_temp = total_fs(biol_no) + biols(biol_no).m(indices_min, indices_max);
//                // This line is expensive - extra calculation of partial_f
//                catches = catches +
//                    (partial_f(fishery_counter, catch_counter, biol_no, indices_min, indices_max) / z_temp) *
//                    (1.0 - exp(-1.0 * z_temp)) *
//                    biols(biol_no).n(indices_min, indices_max); 
//            } 
//            // Calc landings and discards from catches and discards.ratio
//            discards_ratio_temp = fisheries(fishery_counter, catch_counter).discards_ratio(indices_min, indices_max); 
//            landings = catches * (1.0 - discards_ratio_temp);
//            discards = catches * discards_ratio_temp;
//            // Dish out catches to landings and discards
//            // Only into the timestep - we don't want to overwrite other timesteps
//            for (unsigned int quant_count = 1; quant_count <= fisheries(fishery_counter, catch_counter).landings_n().get_nquant(); ++quant_count){
//                // Again - careful with iters
//                for (unsigned int iter_count = 1; iter_count <= fisheries(fishery_counter, catch_counter).landings_n().get_niter(); ++iter_count){
//                    fisheries(fishery_counter, catch_counter).landings_n()(quant_count,year,unit,season,area,iter_count) = landings(quant_count,1,unit,1,area,iter_count); 
//                    fisheries(fishery_counter, catch_counter).discards_n()(quant_count,year,unit,season,area,iter_count) = discards(quant_count,1,unit,1,area,iter_count); 
//                }
//            }
//        }
//    }
//    end = clock();
//    //Rprintf("After updating catches: %f\n", (end - start) / (double)(CLOCKS_PER_SEC));
//    //Rprintf("After updating catches\n");
//    
//    start = clock();
//    //Rprintf("Updating biols\n");
//    // Update biol in next timestep only if within time range
//    if ((next_year <= biol_dim[1]) & (next_season <= biol_dim[3])){
//        unsigned int ssb_year = 0;
//        unsigned int ssb_season = 0;
//        std::vector<unsigned int> ssb_indices_min;
//        std::vector<unsigned int> ssb_indices_max;
//        adouble rec_temp = 0.0;
//        adouble ssb_temp = 0.0;
//        FLQuantAD next_n;
//        // In what age do we put next timestep abundance? If beginning of year, everything is one year older
//        int age_shift = 0;
//        if (next_season == 1){
//            age_shift = 1;
//        }
//        for (unsigned int biol_counter=1; biol_counter <= biols.get_nbiols(); ++biol_counter){
//            // SSB and Recruitment
//            // Get the year and season of the SSB that will result in recruitment in the next timestep
//            // Add one to timestep because we are getting the recruitment in timestep+1
//            unsigned int ssb_timestep = timestep - biols(biol_counter).srr.get_timelag() + 1;
//            if (ssb_timestep < 1){
//                Rcpp::stop("project_timestep: ssb timestep outside of range");
//            }
//            timestep_to_year_season(ssb_timestep, biol_dim[3], ssb_year, ssb_season);
//            // Update min and max_indices for ssb and get it
//            ssb_indices_min = {ssb_year, unit, ssb_season, area, 1};
//            ssb_indices_max = {ssb_year, unit, ssb_season, area, niter};
//            FLQuantAD ssb_flq= ssb(biol_counter, ssb_indices_min, ssb_indices_max);
//            // Get abundances in next timestep from rec+1 to pg
//            z_temp = total_fs(biol_counter) + biols(biol_counter).m(indices_min, indices_max);
//            next_n = biols(biol_counter).n(indices_min, indices_max) * exp(-1.0 * z_temp);
//
//            // Update biol in next timestep using next n = current n * exp(-z) taking care of age shifts and plus groups
//            for (int iter_count = 1; iter_count <= niter; ++iter_count){
//                // Recruitment
//                // rec is calculated in a scalar way - i.e. not passing it a vector of SSBs, so have to do one iter at a time
//                ssb_temp = ssb_flq(1,1,1,1,1,iter_count);
//                //Rprintf("year: %i, season: %i, ssb_year: %i, ssb_season: %i, ssb_temp: %f\n", year, season, ssb_year, ssb_season, Value(ssb_temp));
//                rec_temp = biols(biol_counter).srr.eval_model(ssb_temp, next_year, 1, next_season, 1, iter_count);
//                //Rprintf("iter_count %i; ssb_temp %f; rec_temp %f\n", iter_count, Value(ssb_temp), Value(rec_temp));
//                // Apply the residuals to rec_temp
//                // Use of if statement is OK because for each taping it will only branch the same way (depending on residuals_mult)
//                if (biols(biol_counter).srr.get_residuals_mult()){
//                    rec_temp = rec_temp * biols(biol_counter).srr.get_residuals()(1,next_year,1,next_season,1,iter_count);
//                }
//                else{
//                    rec_temp = rec_temp + biols(biol_counter).srr.get_residuals()(1,next_year,1,next_season,1,iter_count);
//                }
//
//                // Update ages rec+1 upto final age 
//                for (int quant_count = 2; quant_count <= biol_dim[0]; ++quant_count){
//                    biols(biol_counter).n()(quant_count, next_year, 1, next_season, 1, iter_count) = next_n(quant_count - age_shift, 1,unit,1,area,iter_count);
//                }
//                // Fix Plus Group and add in recruitment
//                // If next_season is start of the year assume last age is a plusgroup and abundance in first age group is ONLY recruitment
//                if (next_season == 1){
//                    biols(biol_counter).n()(biol_dim[0], next_year, 1, next_season, 1, iter_count) = biols(biol_counter).n()(biol_dim[0], next_year, 1, next_season, 1, iter_count) + next_n(biol_dim[0],1,unit,1,area,iter_count);
//                    // Rec in first age is only the calculated recruitment
//                    biols(biol_counter).n()(1,next_year,1,next_season,1,iter_count) = rec_temp;
//                }
//                // Otherwise add in recruitment to first age group
//                else {
//                    // Add rec to the existing age 1
//                    biols(biol_counter).n()(1, next_year, 1, next_season, 1, iter_count) = next_n(1, 1, unit, 1, area, iter_count) + rec_temp;
//                }
//            }
//        }
//    }
//    end = clock();
//    //Rprintf("After updating biols: %f\n", (end - start) / (double)(CLOCKS_PER_SEC));
//
//
//    end2 = clock();
//    Rprintf("project_timestep CPU time: %f\n", (end2 - start2) / (double)(CLOCKS_PER_SEC));
//     
//    Rprintf("Leaving project_timestep\n");
//    
//    return; 
//}
//*/
//
//
//// Returns the indices of the age range, starts at 0
//std::vector<unsigned int> operatingModel::get_target_age_range_indices(const unsigned int target_no, const unsigned int sim_target_no, const unsigned int biol_no) const {
//    std::vector<unsigned int> age_range = ctrl.get_age_range(target_no, sim_target_no);
//    std::vector<unsigned int> age_range_indices(2);
//    // Get age dimnames - then find position of ages from control in them
//    // Convert the age names to a vector of strings
//    std::vector<std::string> age_names = Rcpp::as<std::vector<std::string> >(biols(biol_no).n().get_dimnames()[0]);
//    // Use find() to match names - precheck in R that they exist - if not find, returns the last
//    std::vector<std::string>::iterator age_min_iterator = find(age_names.begin(), age_names.end(), std::to_string(age_range[0]));
//    if(age_min_iterator != age_names.end()){
//        age_range_indices[0] = std::distance(age_names.begin(), age_min_iterator);
//    }
//    else {
//        Rcpp::stop("minAge in control not found in dimnames of FLBiol\n");
//    }
//    std::vector<std::string>::iterator age_max_iterator = find(age_names.begin(), age_names.end(), std::to_string(age_range[1]));
//    if(age_max_iterator != age_names.end()){
//        age_range_indices[1] = std::distance(age_names.begin(), age_max_iterator);
//    }
//    else {
//        Rcpp::stop("maxAge in control not found in dimnames of FLBiol\n");
//    }
//    return age_range_indices;
//}
//
///*! \brief Get the current target value in the operating model
// *
// * Interrogates the control object to get the target type and time, then
// * returns the current value in the operating Model.
// * Any FLQuant range can be specified, i.e. to get values over a range of years and seasons
// * A target can be relative to another fishery / catch / biol.
// * Therefore it is possible to the evaluate the target of the relative fishery / catch / biol.
// * \param target_no References the target column in the control dataframe.
// * \param sim_target_no Which simultaneous target.
// * \param indices_min Vector of minimum FLQuant dims for subsetting.
// * \param indices_max Vector of maximum FLQuant dims for subsetting.
// * \param relative_target Get the relative target (evaluates the target for the relative fishery / catch / biol)
// */
//FLQuantAD operatingModel::eval_target(const unsigned int target_no, const unsigned int sim_target_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max, const bool relative_target) const {
//    Rprintf("In operatingModel::eval_target. target_no: %i sim_target_no: %i\n", target_no, sim_target_no);
//    FLQuantAD out;
//    fwdControlTargetType target_type = ctrl.get_target_type(target_no, sim_target_no);
//    auto niter = ctrl.get_niter();
//    // Do we want to evaluate the relative target?
//    auto fishery_no = 0;
//    auto catch_no = 0;
//    auto biol_no = 0;
//    if (relative_target){
//        Rprintf("Getting values from relative target quantity\n");
//        fishery_no = ctrl.get_target_int_col(target_no, sim_target_no, "relFishery");
//        catch_no = ctrl.get_target_int_col(target_no, sim_target_no, "relCatch");
//        biol_no = ctrl.get_target_int_col(target_no, sim_target_no, "relBiol");
//    }
//    else {
//        fishery_no = ctrl.get_target_int_col(target_no, sim_target_no, "fishery");
//        catch_no = ctrl.get_target_int_col(target_no, sim_target_no, "catch");
//        biol_no = ctrl.get_target_int_col(target_no, sim_target_no, "biol");
//    }
//
//    //// Check patterns of F, C and B to see if they make sense
//    //// You can have 
//    //if ((Rcpp::IntegerVector::is_na(catch_no) ^ Rcpp::IntegerVector::is_na(fishery_no)) | (!(Rcpp::IntegerVector::is_na(biol_no) ^ Rcpp::IntegerVector::is_na(fishery_no)))){
//    //    Rcpp::stop("In operatingModel::eval_target. Either biol_no or (catch_no and fishery_no) or all three must not be NA\n");
//    //}
//
//    // If we have a catch_no, we must also have a fishery_no: XOR! 
//    if(Rcpp::IntegerVector::is_na(catch_no) ^ Rcpp::IntegerVector::is_na(fishery_no)){
//        Rcpp::stop("In operatingModel::eval_target. If you specify a catch_no, you must also specify a fishery_no (relative or not)\n");
//    }
//
//    // Get the output depending on target type
//    switch(target_type){
//        //case target_effort:
//        //break;
//        case target_fbar: {
//            Rprintf("target_fbar\n");
//            // Indices only 5D when passed in - needs age range
//            if (Rcpp::IntegerVector::is_na(biol_no)) {
//                Rcpp::stop("In operatingModel::eval_target. Asking for fbar when biol_na is NA. Not yet implemented.\n");
//            }
//            std::vector<unsigned int> age_range_indices = get_target_age_range_indices(target_no, sim_target_no, biol_no); // indices start at 0
//            std::vector<unsigned int> age_indices_min = indices_min;
//            std::vector<unsigned int> age_indices_max = indices_max;
//            age_indices_min.insert(age_indices_min.begin(), age_range_indices[0] + 1); // +1 because accessor starts at 1 but age indices 0 - sorry
//            age_indices_max.insert(age_indices_max.begin(), age_range_indices[1] + 1);
//            // Is it total fbar on a biol, or fbar of an FLCatch
//            if (!Rcpp::IntegerVector::is_na(biol_no) & !Rcpp::IntegerVector::is_na(catch_no)){
//                Rprintf("fbar is from FLCatch %i in FLFishery %i on biol %i\n", catch_no, fishery_no, biol_no);
//                out = fbar(fishery_no, catch_no, biol_no, age_indices_min, age_indices_max);
//            }
//            else {
//                Rprintf("catch_no is NA. fbar is total fbar from biol %i\n", biol_no);
//                Rprintf("ages_indices_min\n");
//                for (auto i=0; i<6; ++i){
//                    Rprintf("age_indices_min %i\n", age_indices_min[i]);
//                }
//                Rprintf("ages_indices_max\n");
//                for (auto i=0; i<6; ++i){
//                    Rprintf("age_indices_max %i\n", age_indices_max[i]);
//                }
//                out = fbar(biol_no, age_indices_min, age_indices_max);
//            }
//            break;
//        }
//        case target_catch: {
//            Rprintf("target_catch\n");
//            // Is it total catch of a biol, or catch of an FLCatch
//            if (Rcpp::IntegerVector::is_na(biol_no) & !Rcpp::IntegerVector::is_na(catch_no)){
//                Rprintf("biol_no is NA, catch is from FLCatch %i in FLFishery %i\n", catch_no, fishery_no);
//                out = fisheries(fishery_no, catch_no).catches(indices_min, indices_max);
//            }
//            else if (!Rcpp::IntegerVector::is_na(biol_no) & Rcpp::IntegerVector::is_na(catch_no)){
//                Rprintf("catch is total catch from biol %i\n", biol_no);
//                out =  catches(biol_no, indices_min, indices_max);
//            }
//            else {
//                Rcpp::stop("In operatingModel::eval_target. Asking for catch from a particular catch and biol. It's a special case that is not yet implemented. Can you ask for total catch from just the biol instead?\n");
//            }
//            break;
//        }
//        case target_landings: {
//            Rprintf("target_landings\n");
//            if (Rcpp::IntegerVector::is_na(biol_no) & !Rcpp::IntegerVector::is_na(catch_no)){
//                out = fisheries(fishery_no, catch_no).landings(indices_min, indices_max);
//            }
//            else if (!Rcpp::IntegerVector::is_na(biol_no) & Rcpp::IntegerVector::is_na(catch_no)){
//                Rprintf("landings are total landings from biol %i\n", biol_no);
//                out = landings(biol_no, indices_min, indices_max);
//            }
//            else {
//                Rcpp::stop("In operatingModel::eval_target. Asking for landings from a particular catch and biol. It's a special case that is not yet implemented. Can you ask for total landings from just the biol instead?\n");
//            }
//            break;
//        }
//        case target_discards: {
//            Rprintf("target_discards\n");
//            if (Rcpp::IntegerVector::is_na(biol_no) & !Rcpp::IntegerVector::is_na(catch_no)){
//                out = fisheries(fishery_no, catch_no).discards(indices_min, indices_max);
//            }
//            else if (!Rcpp::IntegerVector::is_na(biol_no) & Rcpp::IntegerVector::is_na(catch_no)){
//                Rprintf("discards are total discards from biol %i\n", biol_no);
//                out = discards(biol_no, indices_min, indices_max);
//            }
//            else {
//                Rcpp::stop("In operatingModel::eval_target. Asking for discards from a particular catch and biol. It's a special case that is not yet implemented. Can you ask for total discards from just the biol instead?\n");
//            }
//            break;
//        }
//        case target_srp: {
//            // Spawning reproductive potential
//            Rprintf("target_srp\n");
//            if (Rcpp::IntegerVector::is_na(biol_no)){
//                Rcpp::stop("In operatingModel eval_target. Trying to evaluate SRP target when biol no. is NA. Problem with control object?\n");
//            }
//            else {
//                out = ssb(biol_no, indices_min, indices_max);
//            }
//        break;
//    }
//        case target_biomass: {
//            Rprintf("target_biomass\n");
//            if (Rcpp::IntegerVector::is_na(biol_no)){
//                Rcpp::stop("In operatingModel eval_target. Trying to evaluate biomass target when biol no. is NA. Problem with control object?\n");
//            }
//            else {
//                out = biols(biol_no).biomass(indices_min, indices_max);
//            }
//            break;
//        }
//        default:
//            Rcpp::stop("target_type not found in switch statement - giving up\n");
//        break;
//    }
//    return out;
//}
//
//
///*! \name Get the current target values in the operating model
// */
////@{
///*! \brief Get the current target values in the operating model for all simultaneous targets.
// *
// * Returns a vector of the current simultaneous target values.
// * If the target is relative, the method calculates the current relative state.
// * The values can be compared to the desired target values from get_target_value(). 
// * \param target_no References the target column in the control dataframe. Starts at 1.
// */
//std::vector<adouble> operatingModel::get_target_value_hat(const int target_no) const{
//    Rprintf("In get_target_value_hat target_no\n");
//    auto nsim_target = ctrl.get_nsim_target(target_no);
//    std::vector<adouble> value;
//    for (auto sim_target_count = 1; sim_target_count <= nsim_target; ++sim_target_count){
//        auto sim_target_value = get_target_value_hat(target_no, sim_target_count);
//        value.insert(value.end(), sim_target_value.begin(), sim_target_value.end());
//    }
//    return value;
//} 
//
///*! \brief Get the current target values in the operating model for all simultaneous targets.
// *
// * Returns a vector of the current simultaneous target values.
// * If the target is relative, the method calculates the current relative state.
// * The values can be compared to the desired target values from get_target_value(). 
// * \param target_no References the target column in the control dataframe. Starts at 1.
// * \param sim_target_no References the simultaneous target in the target set. Starts at 1.
// */
//std::vector<adouble> operatingModel::get_target_value_hat(const int target_no, const int sim_target_no) const{
//    Rprintf("In get_target_value_hat sim_target_no\n");
//    auto niter = ctrl.get_niter();
//    std::vector<adouble> value;
//    // Are we dealing with absolute or relative values?
//    // Each target set can have a mix
//    Rcpp::IntegerVector target_rel_year = ctrl.get_target_int_col(target_no, "relYear");
//    Rcpp::IntegerVector target_rel_season = ctrl.get_target_int_col(target_no, "relSeason");
//    //Rprintf("sim_target_no: %i\n", sim_target_no);
//    unsigned int year = ctrl.get_target_int_col(target_no, sim_target_no, "year");
//    unsigned int season = ctrl.get_target_int_col(target_no, sim_target_no, "season");
//    // Indices for subsetting the target values
//    std::vector<unsigned int> indices_min = {year,1,season,1,1};
//    std::vector<unsigned int> indices_max = {year,1,season,1,niter};
//    // Get the current absolute values, i.e. not relative
//    FLQuantAD sim_target_value = eval_target(target_no, sim_target_no, indices_min, indices_max);
//    //Rprintf("sim_target_value: %f\n", Value(sim_target_value(1,1,1,1,1,1)));
//    // Test for a relative target value - is the relative year or season NA?
//    unsigned int rel_year = target_rel_year[sim_target_no-1];
//    unsigned int rel_season = target_rel_season[sim_target_no-1];
//    bool target_rel_year_na = Rcpp::IntegerVector::is_na(rel_year);
//    bool target_rel_season_na = Rcpp::IntegerVector::is_na(rel_season);
//    // Both are either NA, or neither are, if one or other is NA then something has gone wrong (XOR!)
//    if ((target_rel_year_na ^ target_rel_season_na)){
//        Rcpp::stop("in operatingModel::get_target_value sim_target_no. Only one of rel_year or rel_season is NA. Must be neither or both.\n");
//    }
//    // If target is relative we have to get the current actual value (done above) and the actual value it is relative to
//    // Then calculate the relative difference (to be compared to the target)
//    if (!target_rel_year_na){
//        Rprintf("Relative target\n");
//        // Get the value we are relative to
//        std::vector<unsigned int> rel_indices_min = {rel_year,1,rel_season,1,1};
//        std::vector<unsigned int> rel_indices_max = {rel_year,1,rel_season,1,niter};
//        FLQuantAD sim_target_rel_value = eval_target(target_no, sim_target_no, rel_indices_min, rel_indices_max, true);
//        // current proportion: value = value / rel_value
//        std::transform(sim_target_rel_value.begin(), sim_target_rel_value.end(), sim_target_value.begin(), sim_target_value.begin(),
//                    [](adouble x, adouble y){return y / x;});
//    }
//    // Clumsy - if relative we already copied the result back into sim_target_value - could copy direct into value and else{} the next line
//    value.insert(value.end(), sim_target_value.begin(), sim_target_value.end());
//    Rprintf("Leaving get_target_value_hat sim_target\n");
//    return value;
//} 
////@}
//
//
///*! \name Get the desired target values in the operating model
// */
////@{
///*! \brief Get the desired target values in the operating model for all simultaneous targets in a target set.
// * If the target values are relative, the returned values are the proportions from the control object, not the actual values.
// * If values are based on max / min some calculation is required.
// * Returns a vector of the simultaneous target values of a target set. 
// * \param target_no References the target column in the control dataframe. Starts at 1.
// */
//std::vector<double> operatingModel::get_target_value(const int target_no) const{
//    Rprintf("In get_target_value all sim targets\n");
//    auto nsim_target = ctrl.get_nsim_target(target_no);
//    std::vector<double> value;
//    for (auto sim_target_count = 1; sim_target_count <= nsim_target; ++sim_target_count){
//        auto sim_target_value = get_target_value(target_no, sim_target_count);
//        value.insert(value.end(), sim_target_value.begin(), sim_target_value.end());
//    }
//    return value;
//}
///*! \brief Get the desired target values in the operating model for a single simultaneous target
// * If the target values are relative, the returned values are the proportions from the control object, not the actual values.
// * If values are based on max / min some calculation is required.
// * Returns a vector of values of the simultaneous target from the target set. 
// * \param target_no References the target column in the control dataframe. Starts at 1.
// * \param sim_target_no References the target column in the control dataframe. Starts at 1.
// */
//std::vector<double> operatingModel::get_target_value(const int target_no, const int sim_target_no) const{
//    // Pull out values for all iterations for the sim targets from the control object
//    auto niter = ctrl.get_niter();
//    std::vector<double> value(niter);
//    // Are we dealing with a min / max value?
//    // If so we need to get the current state of the operating model to compare with
//    double max_col = ctrl.get_target_num_col(target_no, sim_target_no, "max");
//    double min_col = ctrl.get_target_num_col(target_no, sim_target_no, "min");
//    auto max_na = Rcpp::NumericVector::is_na(max_col);
//    auto min_na = Rcpp::NumericVector::is_na(min_col);
//    if (!max_na | !min_na){
//        // Need to make <double> from <adouble>
//        std::vector<adouble> current_value = get_target_value_hat(target_no, sim_target_no);
//        std::transform(current_value.begin(), current_value.end(), value.begin(), [](adouble x){return Value(x);});
//        if(!max_na){
//            Rprintf("Max target\n");
//            // value is min of current_value and max_value - put result in value
//            std::vector<double> max_value = ctrl.get_target_value(target_no, sim_target_no, 3);
//            std::transform(max_value.begin(), max_value.end(), value.begin(), value.begin(), [](double x, double y) {return std::min(x, y);});
//        }
//        if(!min_na){
//            Rprintf("Min target\n");
//            // value is max of current_value and min_value - put result in correct position in value
//            std::vector<double> min_value = ctrl.get_target_value(target_no, sim_target_no, 1);
//            std::transform(min_value.begin(), min_value.end(), value.begin(), value.begin(), [](double x, double y) {return std::max(x, y);});
//        }
//    }
//    // If not min or max, just get the values from the control object
//    else {
//        value = ctrl.get_target_value(target_no, sim_target_no, 2);
//    }
//    return value;
//}
////@}
//
//
///*
//void operatingModel::run(const double indep_min, const double indep_max){
//    Rprintf("In run\n");
//
//    auto niter = ctrl.get_niter(); // number of iters taken from control object
//    Rprintf("niter: %i\n", niter);
//    auto ntarget = ctrl.get_ntarget();
//    Rprintf("%i targets to solve\n", ntarget);
//    auto neffort = fisheries.get_nfisheries(); // how many efforts, i.e. number of FLFishery objects in OM
//    //Rprintf("%i fisheries to solve effort for\n", neffort);
//
//    auto effort_mult_initial = 1.0; 
//    //auto effort_mult_initial = 2.0; 
//    // Set up effort multipliers - do all iters at once, but keep timesteps, areas, units separate
//    std::vector<double> effort_mult(neffort * niter, effort_mult_initial);
//    std::vector<adouble> effort_mult_ad(neffort * niter, effort_mult_initial);
//
//    // Update N at start of minimum target timestep (not effort timestep)
//    // Assumes it's in the first target set
//    Rcpp::IntegerVector target_timestep = ctrl.get_target_int_col(1, "timestep");
//    // Get min of this
//    auto min_target_timestep = *std::min_element(std::begin(target_timestep), std::end(target_timestep));
//    //Rprintf("min target timestep:%i\n", min_target_timestep);
//    
//    // Assuming that the target timesteps are contiguous,
//    // we first update the biology in the first target timestep.
//    // This ensures that we have abundance numbers in the first timestep of the projection
//    project_biols(min_target_timestep);
//
//    // Loop over targets and solve all simultaneous targets in that target set
//    // e.g. With 2 fisheries with 2 efforts, we can set 2 catch targets to be solved at the same time
//    // Indexing of targets starts at 1
//    for (auto target_count = 1; target_count <= ntarget; ++target_count){
//        Rprintf("\nProcessing target: %i\n", target_count);
//        auto nsim_targets = ctrl.get_nsim_target(target_count);
//        Rprintf("Number of simultaneous targets: %i\n", nsim_targets);
//
//        // Timestep in which we find effort has to be the same for all simultaneous targets in a target set
//        // Cannot just read from the control object because the abundance (biomass, SSB etc.) is determined
//        // by effort in previous timestep
//        // Get time step of first sim target.
//        unsigned int target_effort_timestep = ctrl.get_target_effort_timestep(target_count, 1);
//        unsigned int target_effort_year = 0;
//        unsigned int target_effort_season = 0;
//        timestep_to_year_season(target_effort_timestep, biols(1).n().get_nseason(), target_effort_year, target_effort_season);
//        Rprintf("target_effort_timestep: %i\n", target_effort_timestep);
//        Rprintf("target_effort_year: %i\n", target_effort_year);
//        Rprintf("target_effort_season: %i\n", target_effort_season);
//
//        //Rprintf("Getting desired target values\n");
//        std::vector<double> target_value = get_target_value(target_count);
//
//        // Turn tape on
//        CppAD::Independent(effort_mult_ad);
//        //Rprintf("Turned on tape\n");
//
//        // Update fisheries.effort() with effort multiplier in the effort timestep (area and unit effectively ignored)
//        //Rprintf("Updating effort with multipler\n");
//        for (int fisheries_count = 1; fisheries_count <= fisheries.get_nfisheries(); ++fisheries_count){
//            for (int iter_count = 1; iter_count <= niter; ++ iter_count){
//                fisheries(fisheries_count).effort()(1, target_effort_year, 1, target_effort_season, 1, iter_count) = 
//                    fisheries(fisheries_count).effort()(1, target_effort_year, 1, target_effort_season, 1, iter_count) * 
//                    effort_mult_ad[(fisheries_count - 1) * niter + iter_count - 1];
//            }
//        }
//
//        //Rprintf("About to project\n");
//        // Project fisheries in the target effort timestep
//        // (landings and discards are functions of effort in the effort timestep)
//        project_fisheries(target_effort_timestep); 
//        // Project biology in the target effort timestep plus 1
//        // (biology abundances are functions of effort in the previous timestep)
//        // Only update if there is room
//        if ((target_effort_timestep+1) <= (biols(1).n().get_dim()[1] * biols(1).n().get_dim()[3])){
//            project_biols(target_effort_timestep+1); 
//        }
//        //Rprintf("Back from project\n");
//
//        // Calc error
//        //Rprintf("Getting new state of operating model\n");
//        std::vector<adouble> target_value_hat = get_target_value_hat(target_count); 
//        //Rprintf("Back from target_value_hat\n");
//        //Rprintf("Length of target_value_hat: %i\n", target_value.size());
//        //Rprintf("Length of target_value: %i\n", target_value.size());
//        // Check they are the same length? 
//        std::vector<adouble> error(target_value_hat.size());
//        Rprintf("Calculating error\n");
//        std::transform(target_value.begin(), target_value.end(), target_value_hat.begin(), error.begin(),
//                [](double x, adouble y){return x - y;});
//                //[](double x, adouble y){return (x - y) * (x - y);}); // squared error - doesn't seem so effective
//
//        Rprintf("target_value: %f\n", target_value[0]);
//        Rprintf("target_value_hat: %f\n", Value(target_value_hat[0]));
//        Rprintf("error: %f\n", Value(error[0]));
//        // Stop recording
//        //Rprintf("Turning off tape\n");
//        CppAD::ADFun<double> fun(effort_mult_ad, error);
//        //Rprintf("Turned off tape\n");
//
//        // Solve the target
//        // Reset initial solver value - solver changes the values
//        std::fill(effort_mult.begin(), effort_mult.end(), effort_mult_initial);
//
//        //Rprintf("Calling NR\n");
//        // indep_min and max should be arguments to run and passable from R
//        auto nr_out = newton_raphson(effort_mult, fun, niter, nsim_targets, indep_min, indep_max);
//        //Rprintf("NR done\n");
//
//        // Check nr_out - if not all 1 then something has gone wrong - flag up warning
//
//        
//        Rprintf("effort_mult: %f\n", effort_mult[0]);
//
//        Rprintf("Updating effort with solved effort mult\n");
//        for (int fisheries_count = 1; fisheries_count <= fisheries.get_nfisheries(); ++fisheries_count){
//            for (int iter_count = 1; iter_count <= niter; ++ iter_count){
//                fisheries(fisheries_count).effort()(1, target_effort_year, 1, target_effort_season, 1, iter_count) = 
//                   fisheries(fisheries_count).effort()(1, target_effort_year, 1, target_effort_season, 1, iter_count) * 
//                    effort_mult[(fisheries_count - 1) * niter + iter_count - 1] / effort_mult_initial;
//            }
//        }
//
//        //Rprintf("Projecting again\n");
//        //project_timestep(target_effort_timestep);
//        project_fisheries(target_effort_timestep); 
//        // If space, update biols too
//        if ((target_effort_timestep+1) <= (biols(1).n().get_dim()[1] * biols(1).n().get_dim()[3])){
//            project_biols(target_effort_timestep+1); 
//        }
//        //Rprintf("Done projecting again\n");
//    }
//    Rprintf("Leaving run\n");
//}
//*/
//
//
//
////---------------Target methods ----------------------------
//
//// SSB calculations - Actual SSB that results in recruitment
//
//
//
///*!
// * \name Spawning Stock Biomass of
// *
// * SSB = Wt * mat * N * (-z * spwn)
// *
// * where spwn is the proportion through the timestep where spawning occurs.
// * As F and natural mortality happen simultaneously, the total mortality before spawning is z * spwn.
// *
// */
////@{
///*! \brief SSB of a single biol over a subset of dimensions
// * \param biol_no The index position of the biological stocks in the fwdBiols member (starting at 1).
// * \indices_min minimum indices for subsetting (year - iter, an integer vector of length 5).
// * \indices_max maximum indices for subsetting (year - iter, an integer vector of length 5).
// */
//FLQuantAD operatingModel::ssb(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const{
//    if(indices_min.size() != 5 | indices_max.size() != 5){
//        Rcpp::stop("In ssb() subset method. indices_min and max not of length 5");
//    }
//    // Add quant dim to indices
//    std::vector<unsigned int> quant_indices_min = indices_min;
//    quant_indices_min.insert(quant_indices_min.begin(), 1);
//    std::vector<unsigned int> dims = biols(biol_no).n().get_dim();
//    std::vector<unsigned int> quant_indices_max = indices_max;
//    quant_indices_max.insert(quant_indices_max.begin(), dims[0]);
//    FLQuantAD ssb = quant_sum(biols(biol_no).wt(quant_indices_min, quant_indices_max) * biols(biol_no).fec(quant_indices_min, quant_indices_max) * biols(biol_no).n(quant_indices_min, quant_indices_max) * exp(-1.0*((total_f(biol_no, quant_indices_min, quant_indices_max) + biols(biol_no).m(quant_indices_min, quant_indices_max)) * biols(biol_no).spwn(quant_indices_min, quant_indices_max))));
//    return ssb;
//}
//
///*! \brief SSB of a single biol over a subset of dimensions
// * \param biol_no The index position of the biological stocks in the fwdBiols member (starting at 1).
// * \param total_f Total instantaneous fishing mortality on the biol (used for efficiency) 
// * \indices_min Minimum indices for subsetting.
// * \indices_max Maximum indices for subsetting.
// */
////FLQuantAD operatingModel::ssb(const int biol_no, &FLQuantAD total_f, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
////    // Check dims of total_f match those of others (might not need to)
////    FLQuantAD dummy;
////    return dummy;
////}
//
///*! \brief SSB of a single biol over all dimensions
// * \param biol_no The index position of the biological stocks in the fwdBiols member (starting at 1).
// */
//FLQuantAD operatingModel::ssb(const int biol_no) const {
//    std::vector<unsigned int> indices_min {1,1,1,1,1};
//    //std::vector<unsigned int> indices_max = Rcpp::as<std::vector<unsigned int>>(biols(biol_no).n().get_dim());
//    std::vector<unsigned int> indices_max = biols(biol_no).n().get_dim();
//    indices_max.erase(indices_max.begin(), indices_max.begin()+1);
//    return ssb(biol_no, indices_min, indices_max);
//}
////@}
//
///*
//// Return all iterations but single timestep as an FLQuant
//FLQuantAD operatingModel::ssb(const int timestep, const int unit, const int area, const int biol_no) const {
//    FLQuantAD full_ssb = ssb(biol_no);
//    int year = 0;
//    int season = 0;
//    timestep_to_year_season(timestep, full_ssb, year, season);
//    FLQuantAD out = full_ssb(1,1,year,year,unit,unit,season,season,area,area,1,full_ssb.get_niter());
//    return out;
//}
//
//// Return a single value given timestep, unit, area and iter
//adouble operatingModel::ssb(const int timestep, const int unit, const int area, const int iter, const int biol_no) const {
//    FLQuantAD full_ssb = ssb(biol_no);
//    int year = 0;
//    int season = 0;
//    timestep_to_year_season(timestep, full_ssb, year, season);
//    adouble out = full_ssb(1,year,unit,season,area,iter);
//    return out;
//}
//
//// Return a single value given year, season, unit, area and iter
//adouble operatingModel::ssb(const int year, const int unit, const int season, const int area, const int iter, const int biol_no) const {
//    FLQuantAD full_ssb = ssb(biol_no);
//    adouble out = full_ssb(1,year,unit,season,area,iter);
//    return out;
//}
//*/
//
//
//// indices_min and indices_max - these are indices starting at 1, not dimnames
//// i.e. the first dimension does not hold actual ages, just indices
//FLQuantAD operatingModel::fbar(const int fishery_no, const int catch_no, const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    // Get the F, then mean over the first dimension
//    FLQuantAD f = partial_f(fishery_no, catch_no, biol_no, indices_min, indices_max); 
//    FLQuantAD fbar = quant_mean(f);
//    return fbar;
//}
//
//FLQuantAD operatingModel::fbar(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    // Get the F, then mean over the first dimension
//    FLQuantAD f = total_f(biol_no, indices_min, indices_max); 
//    FLQuantAD fbar = quant_mean(f);
//    return fbar;
//}
//
//
///*! \brief Subset the total landings from a single biol 
// * If the catch that catches the biol only catches from that biol, the landings are taken directly from the landings in the FLCatch object, i.e. they are not recalculated. Therefore, if effort or something has changed, this method will not reflect those changes.
// * If a biol is fished by a catch that also catches from another biol (e.g. a catch fishing on two sub stocks)
// * the landings from each biol have to be recalculated as there is no way of splitting the total catch
// * into the catches from each biol.
// * The sum of the (recalculated) landings may not equal the total landings already stored in the object.
// *
// * \param biol_no Position of the chosen biol in the biols list
// * \param indices_min minimum indices for subsetting (year - iter, integer vector of length 5)
// * \param indices_max maximum indices for subsetting (year - iter, integer vector of length 5)
// */
//FLQuantAD operatingModel::landings(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    if(indices_min.size() != 5 | indices_max.size() != 5){
//        Rcpp::stop("In operatingModel landings on a biol subset method. indices_min and max must be of length 5\n");
//    }
//    // Indices for the full FLQ, i.e. including first dimension
//    std::vector<unsigned int> dim = biols(biol_no).n().get_dim();
//    std::vector<unsigned int> quant_indices_min = indices_min;
//    std::vector<unsigned int> quant_indices_max = indices_max;
//    quant_indices_min.insert(quant_indices_min.begin(), 1);
//    quant_indices_max.insert(quant_indices_max.begin(), dim[0]);
//    unsigned int fishery_no;
//    unsigned int catch_no;
//    // Empty quant for storage
//    FLQuantAD total_landings(1, indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1); 
//    total_landings.fill(0.0);
//    // Get the Fishery / Catches that catch the biol
//    const Rcpp::IntegerMatrix FC =  ctrl.get_FC(biol_no);
//    // Do any of these FCs catch another biol 
//    int nbiols = 0;
//    std::vector<int> biols_fished;
//    // Loop over the FCs that catch from that biol 
//    for (int FC_counter = 0; FC_counter < FC.nrow(); ++FC_counter){
//        // What biols are also fished by that FC
//        biols_fished = ctrl.get_B(FC(FC_counter, 0), FC(FC_counter, 1)); 
//        // If that FC only catches that biol, add catches into total - no need to partition
//        if(biols_fished.size() == 1) { 
//            total_landings = total_landings + fisheries(FC(FC_counter,0), FC(FC_counter,1)).landings(indices_min, indices_max);
//        }
//        // Get partition of the total landings from the FC that comes from the biol
//        // Sadly this involves the baranov equation again
//        else {
//            // Slow as total_f calls partial_f so we end up calling it again
//            FLQuantAD z = total_f(biol_no, quant_indices_min, quant_indices_max) + biols(biol_no).m(quant_indices_min, quant_indices_max);
//            FLQuantAD catch_n =  (partial_f(FC(FC_counter,0), FC(FC_counter,1), biol_no, quant_indices_min, quant_indices_max) / z) * (1.0 - exp(-1.0 * z)) * biols(biol_no).n(quant_indices_min, quant_indices_max);
//            FLQuantAD landings_n = catch_n * (1.0 - fisheries(FC(FC_counter,0), FC(FC_counter,1)).discards_ratio(quant_indices_min, quant_indices_max));
//            FLQuantAD landings = quant_sum(landings_n * fisheries(FC(FC_counter,0), FC(FC_counter,1)).landings_wt(quant_indices_min, quant_indices_max));
//            total_landings = total_landings + landings;
//        }
//    }
//    return total_landings;
//}
//
///*! \brief Subset the total discards from a single biol 
// * If the catch that catches the biol only catches from that biol, the discards are taken directly from the discards in the FLCatch object, i.e. they are not recalculated. Therefore, if effort or something has changed, this method will not reflect those changes.
// * If a biol is fished by a catch that also catches from another biol (e.g. a catch fishing on two sub stocks)
// * the discards from each biol have to be recalculated as there is no way of splitting the total catch
// * into the catches from each biol.
// * The sum of the (recalculated) discards may not equal the total discards already stored in the object.
// *
// * \param biol_no Position of the chosen biol in the biols list
// * \param indices_min minimum indices for subsetting (year - iter, integer vector of length 5)
// * \param indices_max maximum indices for subsetting (year - iter, integer vector of length 5)
// */
//FLQuantAD operatingModel::discards(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    if(indices_min.size() != 5 | indices_max.size() != 5){
//        Rcpp::stop("In operatingModel discards on a biol subset method. indices_min and max must be of length 5\n");
//    }
//    // Indices for the full FLQ, i.e. including first dimension
//    std::vector<unsigned int> dim = biols(biol_no).n().get_dim();
//    std::vector<unsigned int> quant_indices_min = indices_min;
//    std::vector<unsigned int> quant_indices_max = indices_max;
//    quant_indices_min.insert(quant_indices_min.begin(), 1);
//    quant_indices_max.insert(quant_indices_max.begin(), dim[0]);
//    unsigned int fishery_no;
//    unsigned int catch_no;
//    // Empty quant for storage
//    FLQuantAD total_discards(1, indices_max[0] - indices_min[0] + 1, indices_max[1] - indices_min[1] + 1, indices_max[2] - indices_min[2] + 1, indices_max[3] - indices_min[3] + 1, indices_max[4] - indices_min[4] + 1); 
//    total_discards.fill(0.0);
//    // Get the Fishery / Catches that catch the biol
//    const Rcpp::IntegerMatrix FC =  ctrl.get_FC(biol_no);
//    // Do any of these FCs catch another biol 
//    int nbiols = 0;
//    std::vector<int> biols_fished;
//    // Loop over the FCs that catch from that biol 
//    for (int FC_counter = 0; FC_counter < FC.nrow(); ++FC_counter){
//        // What biols are also fished by that FC
//        biols_fished = ctrl.get_B(FC(FC_counter, 0), FC(FC_counter, 1)); 
//        // If that FC only catches that biol, add catches into total - no need to partition
//        if(biols_fished.size() == 1) { 
//            total_discards = total_discards + fisheries(FC(FC_counter,0), FC(FC_counter,1)).discards(indices_min, indices_max);
//        }
//        // Get partition of the total discards from the FC that comes from the biol
//        // Sadly this involves the baranov equation again
//        else {
//            // Slow as total_f calls partial_f so we end up calling it again
//            FLQuantAD z = total_f(biol_no, quant_indices_min, quant_indices_max) + biols(biol_no).m(quant_indices_min, quant_indices_max);
//            FLQuantAD catch_n =  (partial_f(FC(FC_counter,0), FC(FC_counter,1), biol_no, quant_indices_min, quant_indices_max) / z) * (1.0 - exp(-1.0 * z)) * biols(biol_no).n(quant_indices_min, quant_indices_max);
//            FLQuantAD discards_n = catch_n * (fisheries(FC(FC_counter,0), FC(FC_counter,1)).discards_ratio(quant_indices_min, quant_indices_max));
//            FLQuantAD discards = quant_sum(discards_n * fisheries(FC(FC_counter,0), FC(FC_counter,1)).discards_wt(quant_indices_min, quant_indices_max));
//            total_discards = total_discards + discards;
//        }
//    }
//    return total_discards;
//}
//
///*! \brief Subset the total catches from a single biol 
// * If the catch that catches the biol only catches from that biol, the catches are taken directly from the catches in the FLCatch object, i.e. they are not recalculated. Therefore, if effort or something has changed, this method will not reflect those changes.
// * If a biol is fished by a catch that also catches from another biol (i.e. a catch fishing on two sub stocks)
// * the catches from each biol have to be recalculated as there is no way of splitting the total catch
// * into the catches from each biol.
// * The sum of the (recalculated) catches may not equal the total catch already stored in the object.
// *
// * \param biol_no Position of the chosen biol in the biols list
// * \param indices_min minimum indices for subsetting (year - iter, integer vector of length 5)
// * \param indices_max maximum indices for subsetting (year - iter, integer vector of length 5)
// */
//FLQuantAD operatingModel::catches(const int biol_no, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max) const {
//    if(indices_min.size() != 5 | indices_max.size() != 5){
//        Rcpp::stop("In operatingModel catches on a biol subset method. indices_min and max must be of length 5\n");
//    }
//    FLQuantAD total_catches = landings(biol_no, indices_min, indices_max) + discards(biol_no, indices_min, indices_max);
//    return total_catches;
//}
//
//
////-------------------------------------------------
//
//
////// Define some pointers to function 
////typedef std::vector<adouble> (*funcPtrAD)(const std::vector<adouble>& x);
////typedef std::vector<double> (*funcPtr)(const std::vector<double>& x);
////
////
////// Returning a whole ADFun might be expensive - return as pointer?
////CppAD::ADFun<double> tape_my_func(std::vector<double>& xin, funcPtrAD fun){
////    // How to easily make an adouble vector?
////    std::vector<adouble> x(xin.size());
////    std::copy(xin.begin(), xin.end(), x.begin());
////    // Tape on
////    CppAD::Independent(x);
////    std::vector<adouble> y = fun(x);
////    CppAD::ADFun<double> f(x, y);
////    return f;
////}
////
////// [[Rcpp::export]]
////Rcpp::List solve_my_func(std::vector<double> xin, Rcpp::XPtr<funcPtrAD> xptr){
////    funcPtrAD func = *xptr; // the typedef
////    CppAD::ADFun<double> fun = tape_my_func(xin, func);
////    std::vector<int> success = newton_raphson(xin, fun, 1, 2, -1e9, 1e9, 50, 1e-12);
////    return Rcpp::List::create(Rcpp::Named("indep") = xin,
////                              Rcpp::Named("success") = success);
////}
////
