//
// Created by Jakob Goldmann on 10.04.25.
//

#ifndef BISAM_RCPP_H
#define BISAM_RCPP_H

// [[Rcpp::export]]
Rcpp::List b_ism_rcpp(
    arma::mat data,
    int i_index,
    int t_index,
    int y_index,
    long long Ndraw,
    long long Nburn,
    std::string b_prior,
    double lambda_b,
    double c0,
    double C0,
    double va,
    double vb,
    double tau,
    bool geweke,
    bool use_phiinit,
    bool const_val,
    bool ife,
    bool tfe,
    bool iis,
    bool sis
);

#endif //BISAM_RCPP_H
