/* 
 * Copyright 2014 FLR Team. Distributed under the GPL 2 or later
 * Maintainer: Finlay Scott, JRC
 */


#include "../../inst/include/fwdBiol.h"

/*-------------------------------------------------------*/

// [[Rcpp::export]]
fwdBiol test_fwdBiol_sexp_constructor(SEXP flb_sexp){
	fwdBiol fwdb(flb_sexp);
	return fwdb;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolAD_sexp_constructor(SEXP flb_sexp){
	fwdBiolAD fwdb(flb_sexp);
	return fwdb;
}

// [[Rcpp::export]]
fwdBiol test_fwdBiol_as_wrap(fwdBiol fwdb){
	return fwdb;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolAD_as_wrap(fwdBiolAD fwdb){
	return fwdb;
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiol_fwdSR_constructor(SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdSR fwsr(model_name, params, residuals, residuals_mult);
    fwdBiol fwb(flb_sexp, fwsr);
	return Rcpp::List::create(Rcpp::Named("fwb", fwb),
				Rcpp::Named("srr",fwb.get_srr()));
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiolAD_fwdSRAD_constructor(SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdSRAD fwsr(model_name, params, residuals, residuals_mult);
    fwdBiolAD fwb(flb_sexp, fwsr);
	return Rcpp::List::create(Rcpp::Named("fwb", fwb),
				Rcpp::Named("srr",fwb.get_srr()));
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiol_FLSR_bits_constructor(SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdBiol fwb(flb_sexp, model_name, params, residuals, residuals_mult);
	return Rcpp::List::create(Rcpp::Named("fwb", fwb),
				Rcpp::Named("srr",fwb.get_srr()));
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiolAD_FLSR_bits_constructor(SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdBiolAD fwb(flb_sexp, model_name, params, residuals, residuals_mult);
	return Rcpp::List::create(Rcpp::Named("fwb", fwb),
				Rcpp::Named("srr",fwb.get_srr()));
}

// [[Rcpp::export]]
fwdBiol test_fwdBiol_copy_constructor(fwdBiol fwdb1){
	fwdBiol fwdb2(fwdb1); // uses copy constructor
    return fwdb2;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolAD_copy_constructor(fwdBiolAD fwdb1){
	fwdBiolAD fwdb2(fwdb1); // uses copy constructor
    return fwdb2;
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiol_copy_constructor2(fwdBiol fwdb1, int quant, int year, int unit, int season, int area, int iter, double value){
	fwdBiol fwdb2(fwdb1); // uses copy constructor
	fwdb2.n()(quant,year,unit,season,area,iter) = value;
	return Rcpp::List::create(Rcpp::Named("fwdb1", fwdb1),
				Rcpp::Named("fwdb2",fwdb2));
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiolAD_copy_constructor2(fwdBiolAD fwdb1, int quant, int year, int unit, int season, int area, int iter, double value){
    adouble value_ad = value;
	fwdBiolAD fwdb2(fwdb1); // uses copy constructor
	fwdb2.n()(quant,year,unit,season,area,iter) = value_ad;
	return Rcpp::List::create(Rcpp::Named("fwdb1", fwdb1),
				Rcpp::Named("fwdb2",fwdb2));
}

// [[Rcpp::export]]
fwdBiol test_fwdBiol_assignment_operator(fwdBiol fwdb1){
	fwdBiol fwdb2 = fwdb1; 
    return fwdb2;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolAD_assignment_operator(fwdBiolAD fwdb1){
	fwdBiolAD fwdb2 = fwdb1; 
    return fwdb2;
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiol_assignment_operator2(fwdBiol fwdb1, int quant, int year, int unit, int season, int area, int iter, double value){
	fwdBiol fwdb2 = fwdb1; 
	fwdb2.n()(quant,year,unit,season,area,iter) = value;
	return Rcpp::List::create(Rcpp::Named("fwdb1", fwdb1),
				Rcpp::Named("fwdb2",fwdb2));
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiolAD_assignment_operator2(fwdBiolAD fwdb1, int quant, int year, int unit, int season, int area, int iter, double value){
    adouble value_ad = value;
	fwdBiolAD fwdb2 = fwdb1; 
	fwdb2.n()(quant,year,unit,season,area,iter) = value_ad;
	return Rcpp::List::create(Rcpp::Named("fwdb1", fwdb1),
				Rcpp::Named("fwdb2",fwdb2));
}

// [[Rcpp::export]]
Rcpp::NumericVector test_fwdBiol_const_get_accessors(const fwdBiol fwdb,int quant, int year, int unit, int season, int area, int iter){
    Rcpp::NumericVector out(6);
    out[0] = fwdb.n()(quant, year, unit, season, area, iter);
    out[1] = fwdb.m()(quant, year, unit, season, area, iter);
    out[2] = fwdb.wt()(quant, year, unit, season, area, iter);
    out[3] = fwdb.fec()(quant, year, unit, season, area, iter);
    out[4] = fwdb.spwn()(quant, year, unit, season, area, iter);
    out[5] = fwdb.mat()(quant, year, unit, season, area, iter);
    return out;
}

// [[Rcpp::export]]
Rcpp::List test_fwdBiol_const_get_accessors_subset(const fwdBiol fwdb, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max){
    return Rcpp::List::create(
        Rcpp::Named("n", fwdb.n(indices_min, indices_max)),
        Rcpp::Named("m", fwdb.m(indices_min, indices_max)),
        Rcpp::Named("wt", fwdb.wt(indices_min, indices_max)),
        Rcpp::Named("fec", fwdb.fec(indices_min, indices_max)),
        Rcpp::Named("spwn", fwdb.spwn(indices_min, indices_max)),
        Rcpp::Named("mat", fwdb.mat(indices_min, indices_max)));
}


// [[Rcpp::export]]
Rcpp::NumericVector test_fwdBiol_get_accessors(fwdBiol fwdb,int quant, int year, int unit, int season, int area, int iter){
    Rcpp::NumericVector out(6);
    out[0] = fwdb.n()(quant, year, unit, season, area, iter);
    out[1] = fwdb.m()(quant, year, unit, season, area, iter);
    out[2] = fwdb.wt()(quant, year, unit, season, area, iter);
    out[3] = fwdb.fec()(quant, year, unit, season, area, iter);
    out[4] = fwdb.spwn()(quant, year, unit, season, area, iter);
    out[5] = fwdb.mat()(quant, year, unit, season, area, iter);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector test_fwdBiolAD_const_get_accessors(const fwdBiolAD fwdb,int quant, int year, int unit, int season, int area, int iter){
    Rcpp::NumericVector out(6);
    //out[0] = fwdb.n()(quant, year, unit, season, area, iter).value();
    out[0] = Value(fwdb.n()(quant, year, unit, season, area, iter));
    out[1] = fwdb.m()(quant, year, unit, season, area, iter);
    out[2] = fwdb.wt()(quant, year, unit, season, area, iter);
    out[3] = fwdb.fec()(quant, year, unit, season, area, iter);
    out[4] = fwdb.spwn()(quant, year, unit, season, area, iter);
    out[5] = fwdb.mat()(quant, year, unit, season, area, iter);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector test_fwdBiolAD_get_accessors(fwdBiolAD fwdb,int quant, int year, int unit, int season, int area, int iter){
    Rcpp::NumericVector out(6);
    //out[0] = fwdb.n()(quant, year, unit, season, area, iter).value();
    out[0] = Value(fwdb.n()(quant, year, unit, season, area, iter));
    out[1] = fwdb.m()(quant, year, unit, season, area, iter);
    out[2] = fwdb.wt()(quant, year, unit, season, area, iter);
    out[3] = fwdb.fec()(quant, year, unit, season, area, iter);
    out[4] = fwdb.spwn()(quant, year, unit, season, area, iter);
    out[5] = fwdb.mat()(quant, year, unit, season, area, iter);
    return out;
}

// [[Rcpp::export]]
fwdBiol test_fwdBiol_set_accessors(fwdBiol fwdb, int quant, int year, int unit, int season, int area, int iter, Rcpp::NumericVector values){
    fwdb.n()(quant, year, unit, season, area, iter) = values[0];
    fwdb.m()(quant, year, unit, season, area, iter) = values[1];
    fwdb.wt()(quant, year, unit, season, area, iter) = values[2];
    fwdb.fec()(quant, year, unit, season, area, iter) = values[3];
    fwdb.spwn()(quant, year, unit, season, area, iter) = values[4];
    fwdb.mat()(quant, year, unit, season, area, iter) = values[5];
    return fwdb;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolAD_set_accessors(fwdBiolAD fwdb, int quant, int year, int unit, int season, int area, int iter, Rcpp::NumericVector values){
    adouble value_ad = values[0];
    fwdb.n()(quant, year, unit, season, area, iter) = value_ad;
    fwdb.m()(quant, year, unit, season, area, iter) = values[1];
    fwdb.wt()(quant, year, unit, season, area, iter) = values[2];
    fwdb.fec()(quant, year, unit, season, area, iter) = values[3];
    fwdb.spwn()(quant, year, unit, season, area, iter) = values[4];
    fwdb.mat()(quant, year, unit, season, area, iter) = values[5];
    return fwdb;
}

// [[Rcpp::export]]
FLQuantAD fwdBiolAD_biomass_FLQ(fwdBiolAD fwdb){
    return fwdb.biomass();
}

// [[Rcpp::export]]
FLQuantAD fwdBiolAD_biomass_subset(fwdBiolAD fwdb, const std::vector<unsigned int> indices_min, const std::vector<unsigned int> indices_max){
    return fwdb.biomass(indices_min, indices_max);
}



/*-------------------------------------------------------*/

// fwdBiols constructors

// List of FLBiol and SR components
// [[Rcpp::export]]
fwdBiols test_fwdBiols_list_constructor(const Rcpp::List biols_in){
    fwdBiols biols(biols_in);
    return biols;
}

// [[Rcpp::export]]
fwdBiolsAD test_fwdBiolsAD_list_constructor(const Rcpp::List biols_in){
    fwdBiolsAD biols(biols_in);
    return biols;
}

// fwdBiol constructor
// [[Rcpp::export]]
fwdBiolsAD test_fwdBiolsAD_fwdBiolAD_constructor(SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdBiolAD fwb(flb_sexp, model_name, params, residuals, residuals_mult);
    fwdBiolsAD fwbs(fwb);
    return fwbs;
}

// fwdBiolAD () get operator
// [[Rcpp::export]]
fwdBiolAD test_fwdBiolsAD_const_get_single_index_accessor(const Rcpp::List fwbs_list, const int element){
	const fwdBiolsAD fwbs(fwbs_list);
    fwdBiolAD fwb = fwbs(element);
    return fwb;
}

// [[Rcpp::export]]
fwdBiolAD test_fwdBiolsAD_get_single_index_accessor(const Rcpp::List fwbs_list, const int element){
	fwdBiolsAD fwbs(fwbs_list);
    fwdBiolAD fwb = fwbs(element);
    return fwb;
}

// [[Rcpp::export]]
double test_fwdBiolsAD_const_get_value_accessor(const Rcpp::List fwbs_list, const int biol_no, int quant, int year, int unit, int season, int area, int iter){
	const fwdBiolsAD fwbs(fwbs_list);
    adouble ad_value_out = fwbs(biol_no).n()(quant,year,unit,season,area,iter);
    return Value(ad_value_out);
}

// [[Rcpp::export]]
double test_fwdBiolsAD_get_value_accessor(const Rcpp::List fwbs_list, const int biol_no, int quant, int year, int unit, int season, int area, int iter){
	fwdBiolsAD fwbs(fwbs_list);
    adouble ad_value_out = fwbs(biol_no).n()(quant,year,unit,season,area,iter);
    return Value(ad_value_out);
}

// [[Rcpp::export]]
fwdBiolsAD test_fwdBiolsAD_set_single_index_accessor(const Rcpp::List fwbs_list, const int element, SEXP flb_sexp, const std::string model_name, const FLQuant params, const FLQuant residuals, const bool residuals_mult){
    fwdBiolAD fwb(flb_sexp, model_name, params, residuals, residuals_mult);
	fwdBiolsAD fwbs(fwbs_list);
    fwbs(element) = fwb;
    return fwbs;
}

// [[Rcpp::export]]
fwdBiolsAD test_fwdBiolsAD_set_value_accessor(const Rcpp::List fwbs_list, const int biol_no, int quant, int year, int unit, int season, int area, int iter, double value){
	fwdBiolsAD fwbs(fwbs_list);
    adouble ad_value = value;
    fwbs(biol_no).n()(quant,year,unit,season,area,iter) = ad_value;
    return fwbs;
}

// copy constructor
// [[Rcpp::export]]
Rcpp::List test_fwdBiolsAD_copy_constructor(const Rcpp::List fwbs_list, const int biol_no, const std::vector<int> dims, const double value){
    fwdBiolsAD fwbs1(fwbs_list);
    fwdBiolsAD fwbs2(fwbs1);
    fwbs1(biol_no).n()(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5]) = value;
	return Rcpp::List::create(Rcpp::Named("fwbs1", fwbs1),
				Rcpp::Named("fwbs2", fwbs2));
}

// assignment operator
// [[Rcpp::export]]
Rcpp::List test_fwdBiolsAD_assignment_operator(const Rcpp::List fwbs_list, const int biol_no, const std::vector<int> dims, const double value){
    fwdBiolsAD fwbs1(fwbs_list);
    fwdBiolsAD fwbs2;
    fwbs2 = fwbs1;
    fwbs1(biol_no).n()(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5]) = value;
	return Rcpp::List::create(Rcpp::Named("fwbs1", fwbs1),
				Rcpp::Named("fwbs2", fwbs2));
}

// iterators
// [[Rcpp::export]]
Rcpp::List test_fwdBiolsAD_const_iterator(const Rcpp::List biols_in){
    Rcpp::List out;
    const fwdBiolsAD biols(biols_in);
    for (const auto biol : biols){
        out.push_back(biol.n());
    }
    return out;
}

// [[Rcpp::export]]
fwdBiolsAD test_fwdBiolsAD_iterator(const Rcpp::List biols_in, int quant, int year, int unit, int season, int area, int iter, double value){
    fwdBiolsAD biols(biols_in);
    for (auto& biol : biols){
        biol.n()(quant, year, unit, season, area, iter) = value;
    }
    return biols;
}

