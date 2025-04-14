//
// Created by Jakob Goldmann on 08.04.25.
//
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "biasm_model.h"
#include "bisam_types.h"

// TODO HARD CODED VALUES!!!

namespace bisam {
    // [[Rcpp::plugins(cpp17)]]
    // [[Rcpp::depends(RcppArmadillo)]]

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
        bool use_phiinit,
        bool const_val,
        bool ife,
        bool tfe,
        bool iis,
        bool sis
    ) {
        // Call the C++ function
        BisamResult result = estimate_model(
            data,
            i_index,
            t_index,
            y_index,
            Ndraw,
            Nburn,
            b_prior,
            lambda_b,
            c0,
            C0,
            va,
            vb,
            tau,
            use_phiinit,
            const_val,
            ife,
            tfe,
            iis,
            sis,
            ComputationStrategy::SPLIT_SEQUENTIAL
        );

        // Convert the output struct to an R list
        Rcpp::List ret;
        ret["b_store"]        = result.beta_samples;
        ret["g_store"]        = result.gamma_samples;
        ret["s2_store"]       = result.sigma2_samples;
        ret["w_store"]        = result.indicator_samples;
        ret["w_store_means"]  = result.indicator_means;
        ret["b_store_means"]  = result.beta_means;
        ret["s2_store_means"] = result.sigma2_means;

        return ret;
    }
}
