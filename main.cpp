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
const int RUNS_PER_DATASET          = 1;     // Number of times to run each dataset for timing
const bool APPEND_TO_EXISTING_FILES  = true; // Set to true to append to existing result files

// MCMC Settings
const int MCMC_ITERATIONS = 5000; // Total MCMC iterations
const int MCMC_BURNIN     = 500;  // Burn-in period

// Results Storage Configuration
const std::string EXPERIMENT_NAME = "timelength"; // File name base

// Select which datasets to test (comment/uncomment as needed)
std::vector<std::string> test_datasets = {
    // "tiny",    // Very small for quick testing
    // "small_b", // Small for detailed analysis
    // "small_b", // Another small dataset
    // "med_a",     // Medium size
    // "large_a"    // Large size - takes longer
    // "large_b",
    // "rootfind_stepsize_050",
    // "rootfind_stepsize_075",
    // "rootfind_stepsize_100",
    // "rootfind_stepsize_150",
    // "rootfind_stepsize_300",
    // "rootfind_stepsize_500",

        // "rootfind_timelength_t010",
        // "rootfind_timelength_t020",
        // "rootfind_timelength_t030",
        // "rootfind_timelength_t040",
        // "rootfind_timelength_t050",
        // "rootfind_timelength_t060",
        "rootfind_timelength_t070",
        // "rootfind_timelength_t080",
        // "rootfind_timelength_t090",
        // "rootfind_timelength_t100",
        // "rootfind_timelength_t110",
        // "rootfind_timelength_t120",
};

int main() {
    // Initialize results storage
    bisam::ResultsStorage results_storage(EXPERIMENT_NAME, APPEND_TO_EXISTING_FILES);

    // Get all available datasets
    auto datasets = get_all_datasets();

    std::cout << "BISAM Performance Testing" << std::endl;
    std::cout << "Experiment: " << EXPERIMENT_NAME << std::endl;
    std::cout << "Testing thread counts 1-30" << std::endl;
    std::cout << "Testing " << test_datasets.size() << " datasets with " << RUNS_PER_DATASET << " runs each" << std::endl;
    std::cout << "MCMC: " << MCMC_ITERATIONS << " iterations, " << MCMC_BURNIN << " burn-in" << std::endl;
    std::cout << "Append mode: " << (APPEND_TO_EXISTING_FILES ? "ON" : "OFF") << std::endl;

    // Loop through thread counts from 1 to 30
        // Update RUN_NAME for current thread count



             for (const auto &dataset: datasets) {
            // Skip datasets not in test list
            if (std::find(test_datasets.begin(), test_datasets.end(), dataset.name) == test_datasets.end()) {
                continue;
            }

                 std::string RUN_NAME = dataset.name;

            // Print dataset information
            if (ENABLE_VALIDATION_OUTPUT) {
                BisamValidator::print_dataset_info(dataset);
            } else {
                std::cout << dataset.name << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx << "): ";
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

                // Run BISAM with optimized settings and current thread count
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
                    bisam::ComputationStrategy::SPLIT_SEQUENTIAL, // Use optimized strategy
                    8  // Use current thread count
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
                 results_storage.save_results();
                 results_storage.clear_memory(); // Clear accumulated results to prevent memory issues
        }

        // Save results after each thread count iteration to ensure data is persisted




    // Final save (redundant but ensures everything is saved)
    std::cout << "\nFinal save of all results..." << std::endl;
    results_storage.save_results();

    std::cout << "Completed multithreading performance test (1-30 cores)" << std::endl;

    return 0;
}