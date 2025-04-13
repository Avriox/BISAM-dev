//
// Created by jakob on 4/13/25.
//
#include "../include/modelselection_strategy.h"

namespace bisam {
    arma::Col<int> model_selection_with_strategy(const arma::vec &y, const arma::mat &x, int niter, int thinning,
                                                 int burnin, arma::Col<int> &deltaini_input, bool center, bool scale,
                                                 bool XtXprecomp, double phi, double tau, double priorSkew,
                                                 double prDeltap, arma::vec thinit,
                                                 InitType initpar_type,
                                                 ComputationStrategy strategy, int n) {
        switch (strategy) {
            case ComputationStrategy::STANDARD:
                // Simple passthrough to modelSelection
                return modelSelection(y,
                                      x,
                                      niter,
                                      thinning,
                                      burnin,
                                      deltaini_input,
                                      center,
                                      scale,
                                      XtXprecomp,
                                      phi,
                                      tau,
                                      priorSkew,
                                      prDeltap,
                                      thinit,
                                      initpar_type);

            case ComputationStrategy::SPLIT_SEQUENTIAL: {
                // Prepare the split data
                DataPartition split_data = partition_data(y, x, deltaini_input, thinit, n);

                // Process each part sequentially
                std::vector<arma::Col<int> > results(n);
                for (int part = 0; part < n; part++) {
                    results[part] = modelSelection(
                        split_data.y_parts[part],
                        split_data.common_x,
                        niter,
                        thinning,
                        burnin,
                        split_data.delta_init_parts[part],
                        center,
                        scale,
                        XtXprecomp,
                        phi,
                        tau,
                        priorSkew,
                        prDeltap,
                        split_data.theta_init_parts[part],
                        initpar_type);
                }

                // Combine and return results
                return combine_partition_results(results, split_data.start_columns, split_data.end_columns,
                                                 x.n_cols);
            }

            case ComputationStrategy::SPLIT_PARALLEL: {
                // Prepare the split data (same as sequential)
                DataPartition split_data = partition_data(y, x, deltaini_input, thinit, n);

                // Process each part in parallel (placeholder for now)
                std::vector<arma::Col<int> > results(n);

                // Here you would add parallel execution code
                // For example using OpenMP:
                // #pragma omp parallel for
                for (int part = 0; part < n; part++) {
                    results[part] = modelSelection(
                        split_data.y_parts[part],
                        split_data.common_x,
                        niter,
                        thinning,
                        burnin,
                        split_data.delta_init_parts[part],
                        center,
                        scale,
                        XtXprecomp,
                        phi,
                        tau,
                        priorSkew,
                        prDeltap,
                        split_data.theta_init_parts[part],
                        initpar_type);
                }

                // Combine and return results (same as sequential)
                return combine_partition_results(results, split_data.start_columns, split_data.end_columns,
                                                 x.n_cols);
            }

            default:
                throw std::runtime_error("Unknown model selection strategy");
        }
    }

    DataPartition partition_data( // Previously split_data, clearer name

        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        int num_partitions
    ) {
        DataPartition data;

        size_t n_rows      = y.size();
        size_t n_cols      = x.n_cols;
        size_t thinit_size = theta_init.size(); // Check the actual size of thinit

        // Calculate sizes for each part
        size_t rows_per_part  = n_rows / num_partitions;
        size_t cols_per_part  = n_cols / num_partitions;
        size_t rows_remainder = n_rows % num_partitions;
        size_t cols_remainder = n_cols % num_partitions;

        // Extract the same submatrix of x for all parts (as per requirement)
        // We'll use the first partition's dimensions
        size_t x_start_row = 0;
        size_t x_end_row   = rows_per_part + (rows_remainder > 0 ? 1 : 0) - 1;
        size_t x_start_col = 0;
        size_t x_end_col   = cols_per_part + (cols_remainder > 0 ? 1 : 0) - 1;

        // Extract the submatrix for x (same for all parts)
        data.common_x = x.submat(x_start_row, x_start_col, x_end_row, x_end_col);

        // Prepare vectors for all parts
        data.y_parts.resize(num_partitions);
        data.delta_init_parts.resize(num_partitions);
        data.theta_init_parts.resize(num_partitions);
        data.start_columns.resize(num_partitions);
        data.end_columns.resize(num_partitions);

        // Check if thinit needs to be split
        bool split_thinit = (thinit_size == n_cols);

        // Split y and deltaini_input
        for (int part = 0; part < num_partitions; part++) {
            // Calculate row range for this part
            size_t start_row = part * rows_per_part + std::min(static_cast<size_t>(part), rows_remainder);
            size_t end_row   = (part + 1) * rows_per_part + std::min(static_cast<size_t>(part + 1), rows_remainder) - 1;

            // Calculate column range for this part
            size_t start_col = part * cols_per_part + std::min(static_cast<size_t>(part), cols_remainder);
            size_t end_col   = (part + 1) * cols_per_part + std::min(static_cast<size_t>(part + 1), cols_remainder) - 1;

            // Extract the corresponding part of y
            data.y_parts[part] = y.subvec(start_row, end_row);

            // Extract the corresponding part of deltaini_input
            data.delta_init_parts[part] = delta_initial.subvec(start_col, end_col);

            // Handle thinit appropriately based on its size
            if (split_thinit) {
                // If thinit has the same length as x.n_cols, split it accordingly
                data.theta_init_parts[part] = theta_init.subvec(start_col, end_col);
            } else if (thinit_size > 0) {
                // If thinit is not empty but doesn't match x.n_cols, use the full vector
                data.theta_init_parts[part] = theta_init;
            } else {
                // If thinit is empty, create an empty vector
                data.theta_init_parts[part] = arma::vec();
            }

            // Store the column indices for later reconstruction
            data.start_columns[part] = start_col;
            data.end_columns[part]   = end_col;
        }

        return data;
    }

    arma::Col<int> combine_partition_results( // Previously combine_results, more specific
        const std::vector<arma::Col<int> > &results,
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    ) {
        arma::Col<int> combined_result(total_columns, arma::fill::zeros);

        for (size_t part = 0; part < results.size(); part++) {
            for (size_t j = 0; j <= (end_columns[part] - start_columns[part]); j++) {
                combined_result(start_columns[part] + j) = results[part](j);
            }
        }

        return combined_result;
    }
}
