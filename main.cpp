#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <RcppArmadillo.h>

#include "include/biasm_model.h"
#include "r-testing-files/simulation_data/simulation_datasets.h"
#include "validation.h"
#include "result_storage.h"

// ========================================================================
// CONFIGURATION SECTION
// ========================================================================
const bool ENABLE_VALIDATION_OUTPUT = false; // Set to false to disable all validation printing
const int RUNS_PER_DATASET          = 10;     // Number of times to run each dataset for timing
const bool APPEND_TO_EXISTING_FILES  = true; // Set to true to append to existing result files

// MCMC Settings
const int MCMC_ITERATIONS = 5000; // Total MCMC iterations
const int MCMC_BURNIN     = 500;  // Burn-in period

// Results Storage Configuration
const std::string EXPERIMENT_NAME = "cumulative_gains_study"; // File name base
const std::string RUN_NAME = "cpp_split_newton_parallel_thopt_o3-native-loop_flto";   // This run's identifier - CHANGE THIS FOR EACH RUN

// Select which datasets to test (comment/uncomment as needed)
std::vector<std::string> test_datasets = {
    // "tiny",    // Very small for quick testing
    // "small_a", // Small for detailed analysis
    "small_b", // Another small dataset
    // "med_a",     // Medium size
    // "large_a"    // Large size - takes longer
};

int main() {
    // Initialize results storage
    bisam::ResultsStorage results_storage(EXPERIMENT_NAME, APPEND_TO_EXISTING_FILES);

    // Get all available datasets
    auto datasets = get_all_datasets();

    if (!ENABLE_VALIDATION_OUTPUT) {
        std::cout << "BISAM Performance Testing" << std::endl;
        std::cout << "Experiment: " << EXPERIMENT_NAME << std::endl;
        std::cout << "Run: " << RUN_NAME << std::endl;
        std::cout << "Testing " << test_datasets.size() << " datasets with " << RUNS_PER_DATASET << " runs each" << std::endl;
        std::cout << "MCMC: " << MCMC_ITERATIONS << " iterations, " << MCMC_BURNIN << " burn-in" << std::endl;
        std::cout << "Append mode: " << (APPEND_TO_EXISTING_FILES ? "ON" : "OFF") << std::endl;
    }

    for (const auto &dataset: datasets) {
        // Skip datasets not in test list
        if (std::find(test_datasets.begin(), test_datasets.end(), dataset.name) == test_datasets.end()) {
            continue;
        }

        // Print dataset information
        if (ENABLE_VALIDATION_OUTPUT) {
            BisamValidator::print_dataset_info(dataset);
        } else {
            std::cout << "\n" << dataset.name << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx << "): ";
        }

        // Variables to store validation results from first run
        ValidationResults first_run_validation;
        bool validation_done = false;

        // Run the dataset multiple times for timing
        for (int run = 0; run < RUNS_PER_DATASET; run++) {
            if (ENABLE_VALIDATION_OUTPUT) {
                std::cout << "\nRun " << (run + 1) << "/" << RUNS_PER_DATASET << " for dataset " << dataset.name << std::endl;
            }

            // Start timing for this individual run
            auto start_time = std::chrono::high_resolution_clock::now();

            // Run BISAM with optimized settings
            bisam::BisamResult result = bisam::estimate_model(
                dataset.data,
                0, 1, 2, // column indices: unit, time, y
                MCMC_ITERATIONS,
                MCMC_BURNIN,
                "g",      // prior type
                100.0,    // prior scale
                0.001, 0.001, // tolerances
                1.0, 1.0, 1.0, // prior parameters
                true,     // include indicators
                false, false, false, // no fixed effects for performance testing
                true, true,   // IIS, SIS
                2, 1, 2, 0,   // method parameters
                100000, 1, 1, // MCMC parameters
                0.01, 0.01,   // alpha, lambda
                1, 0.5, {1, 1}, // additional parameters
                bisam::ComputationStrategy::SPLIT_PARALLEL // Use optimized strategy
            );

            // End timing for this individual run
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double execution_time_ms = duration.count() / 1000.0; // Convert to milliseconds

            // Store timing result
            results_storage.add_timing_result(dataset.name, RUN_NAME, run + 1, execution_time_ms,
                                            dataset.n, dataset.t, dataset.nx);

            // Store estimation results
            results_storage.add_estimation_results(dataset.name, RUN_NAME, run + 1, result, dataset);

            // Always do validation but control printing
            ValidationResults validation = BisamValidator::validate_results(result, dataset, 0.5, ENABLE_VALIDATION_OUTPUT);
            results_storage.add_validation_summary(dataset.name, RUN_NAME, run + 1, validation);

            // Handle output based on setting
            if (ENABLE_VALIDATION_OUTPUT) {
                if (run == 0) {
                    first_run_validation = validation;
                    BisamValidator::print_validation_summary(first_run_validation, dataset.name);
                    validation_done = true;
                }
            } else {
                // Just show timing - no validation output
                if (run == 0) {
                    std::cout << std::fixed << std::setprecision(1) << execution_time_ms << "ms";
                    first_run_validation = validation;
                    validation_done = true;
                } else {
                    std::cout << ", " << std::fixed << std::setprecision(1) << execution_time_ms << "ms";
                }
            }


        }

        // Print brief summary for this dataset
        if (!ENABLE_VALIDATION_OUTPUT && validation_done) {
            std::cout << " | Beta: " << (first_run_validation.beta_max_error < 0.5 ? "OK" : "FAIL")
                      << " | Break: " << (first_run_validation.break_detected_correctly ? "OK" : "FAIL") << std::endl;
        }
    }

    // Save all results to files
    std::cout << "\nSaving results..." << std::endl;
    results_storage.save_results();

    std::cout << "Completed run: " << RUN_NAME << std::endl;

    return 0;
}

// ========================================================================
// USAGE INSTRUCTIONS
// ========================================================================

/*
SIMPLE USAGE:

1. Set EXPERIMENT_NAME to the study name (this determines file names)
2. Set RUN_NAME to identify this specific run/version
3. Set APPEND_TO_EXISTING_FILES = true to add to existing files
4. Choose your datasets in test_datasets
5. Compile and run

For cumulative gains study:
- Keep EXPERIMENT_NAME = "cumulative_gains_study"
- Change RUN_NAME for each version: "R_original", "cpp_basic", "cpp_optimized", etc.
- Set APPEND_TO_EXISTING_FILES = true
- Results go to cumulative_gains_study_timing.csv with RUN_NAME stored in each row

For root finding comparison:
- Set EXPERIMENT_NAME = "root_finding_comparison"
- Change RUN_NAME for each algorithm: "jenkins_traub", "newton_raphson", etc.
- Recompile with different algorithms and run
- Results append to root_finding_comparison_timing.csv

For parallelization study:
- Set EXPERIMENT_NAME = "parallelization_study"
- Change RUN_NAME for each core count: "1_core", "4_cores", "8_cores", etc.
- Results append to parallelization_study_timing.csv
*/