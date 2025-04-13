//
// Created by jakob on 4/13/25.
//

#ifndef MODELSELECTION_STRATEGY_H
#define MODELSELECTION_STRATEGY_H
#include <vector>
#include <RcppArmadillo.h>
#include "bisam_types.h"
#include "mombf_bridge.h"

namespace bisam {
    arma::Col<int> model_selection_with_strategy(const arma::vec &y, const arma::mat &x, int niter, int thinning,
                                                 int burnin, arma::Col<int> &deltaini_input, bool center, bool scale,
                                                 bool XtXprecomp, double phi, double tau, double priorSkew,
                                                 double prDeltap, arma::vec thinit,
                                                 InitType initpar_type,
                                                 ComputationStrategy strategy, int n = 3);

    DataPartition partition_data(

        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        int num_partitions
    );

    arma::Col<int> combine_partition_results( // Previously combine_results, more specific
        const std::vector<arma::Col<int> > &results,
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    );
}

#endif //MODELSELECTION_STRATEGY_H
