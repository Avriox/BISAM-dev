//
// Created by jakob on 4/13/25.
//
#ifndef BIASM_MODEL_H
#define BIASM_MODEL_H
#include <string>

#include "bisam_types.h"
#include <RcppArmadillo.h>
#include "../include/utils.h"

namespace bisam {
    BisamResult estimate_model(
        arma::mat &data,
        int i_index,
        int t_index,
        int y_index,
        long long num_draws,
        long long num_burnin,
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
        bool sis,
        ComputationStrategy strategy
    );
}


#endif //BIASM_MODEL_H
