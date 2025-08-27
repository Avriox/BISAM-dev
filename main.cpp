#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <RcppArmadillo.h>

#include "include/biasm_model.h"
#include "r-testing-files/simulation_data/simulation_datasets.h"
#include "validation_include/validation.h"
#include "validation_include/result_storage.h"

// ========================================================================
// CONFIGURATION SECTION
// ========================================================================
const bool ENABLE_VALIDATION_OUTPUT = false; // Set to false to disable all validation printing
const int RUNS_PER_DATASET          = 1;     // Number of times to run each dataset for timing
const bool APPEND_TO_EXISTING_FILES = true;  // Set to true to append to existing result files

// MCMC Settings
const int MCMC_ITERATIONS = 2000; // Total MCMC iterations
const int MCMC_BURNIN     = 500;  // Burn-in period

// Results Storage Configuration
const std::string EXPERIMENT_NAME = "step_size_investigation"; // File name base

// Select which datasets to test (comment/uncomment as needed)
std::vector<std::string> test_datasets = {
    // "tiny", // Very small for quick testing
    // "small_b", // Small for detailed analysis
    // "small_b", // Another small dataset
    // "med_a",     // Medium size
    // "large_a"    // Large size - takes longer
    // "large_b",
    // "rootfind_stepsize_050",
    "rootfind_stepsize_300_001"
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
    // "rootfind_timelength_t070",
    // "rootfind_timelength_t080",
    // "rootfind_timelength_t090",
    // "rootfind_timelength_t100",
    // "rootfind_timelength_t110",
    // "rootfind_timelength_t120",
};

int main(int argc, char *argv[]) {
    // Check if dataset name provided as command line argument
    std::vector<std::string> active_test_datasets;
    int start_run               = 1;         // Default start from run 1
    std::string run_name_prefix = "jenkins"; // Default prefix

    // Debug: Print all arguments
    std::cout << "Arguments received: " << argc << std::endl;
    for (int i = 0; i < argc; i++) {
        std::cout << "  argv[" << i << "] = '" << argv[i] << "'" << std::endl;
    }

    if (argc > 1) {
        // Use dataset from command line argument
        active_test_datasets = {argv[1]};
        std::cout << "Running single dataset from command line: " << argv[1] << std::endl;

        // Check if start run number provided as second argument
        if (argc > 2) {
            start_run = std::atoi(argv[2]);
            if (start_run < 1 || start_run > RUNS_PER_DATASET) {
                std::cerr << "ERROR: Start run must be between 1 and " << RUNS_PER_DATASET << std::endl;
                return 5;
            }
            std::cout << "Starting from run " << start_run << std::endl;
        }

        // Check if run name prefix provided as third argument
        if (argc > 3) {
            run_name_prefix = argv[3];
            std::cout << "Using run name prefix from command line: '" << run_name_prefix << "'" << std::endl;
        } else {
            std::cout << "Using default run name prefix: '" << run_name_prefix << "'" << std::endl;
        }
    } else {
        // Use hardcoded list (existing behavior)
        active_test_datasets = test_datasets;
    }

    try {
        // Initialize results storage
        bisam::ResultsStorage results_storage(EXPERIMENT_NAME, APPEND_TO_EXISTING_FILES);

        // Get all available datasets
        auto datasets = get_all_datasets();

        std::cout << "BISAM Performance Testing" << std::endl;
        std::cout << "Experiment: " << EXPERIMENT_NAME << std::endl;
        std::cout << "Testing thread counts 1-30" << std::endl;
        std::cout << "Testing " << active_test_datasets.size() << " datasets with " << RUNS_PER_DATASET << " runs each"
                << std::endl;
        std::cout << "MCMC: " << MCMC_ITERATIONS << " iterations, " << MCMC_BURNIN << " burn-in" << std::endl;
        std::cout << "Append mode: " << (APPEND_TO_EXISTING_FILES ? "ON" : "OFF") << std::endl;

        bool found_dataset = false;
        int completed_runs = 0;

        // Loop through thread counts from 1 to 30
        // Update RUN_NAME for current thread count

        for (const auto &dataset: datasets) {
            // Skip datasets not in test list
            if (std::find(active_test_datasets.begin(), active_test_datasets.end(), dataset.name) ==
                active_test_datasets.end()) {
                continue;
            }

            found_dataset = true;

            // Extract step size from dataset name (e.g., "rootfind_stepsize_050" -> "050")
            std::string step_size    = "000"; // Default fallback
            std::string dataset_name = dataset.name;
            size_t pos               = dataset_name.find_last_of("_");
            if (pos != std::string::npos && pos + 1 < dataset_name.length()) {
                step_size = dataset_name.substr(pos + 1);
            }

            std::string RUN_NAME = run_name_prefix + "_" + step_size;

            // Debug: Print RUN_NAME construction
            std::cout << "DEBUG: Constructing RUN_NAME:" << std::endl;
            std::cout << "  dataset.name: '" << dataset.name << "'" << std::endl;
            std::cout << "  step_size: '" << step_size << "'" << std::endl;
            std::cout << "  run_name_prefix: '" << run_name_prefix << "'" << std::endl;
            std::cout << "  final RUN_NAME: '" << RUN_NAME << "'" << std::endl;

            // Print dataset information
            if (ENABLE_VALIDATION_OUTPUT) {
                BisamValidator::print_dataset_info(dataset);
            } else {
                std::cout << dataset.name << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx <<
                        "): ";
            }

            // Variables to store validation results from first run
            ValidationResults first_run_validation;
            bool validation_done = false;

            // Run the dataset multiple times for timing
            for (int run = start_run - 1; run < RUNS_PER_DATASET; run++) {
                // start_run is 1-based, loop is 0-based
                if (ENABLE_VALIDATION_OUTPUT) {
                    std::cout << "\nRun " << (run + 1) << "/" << RUNS_PER_DATASET << " for dataset " << dataset.name <<
                            std::endl;
                }

                try {
                    // Start timing for this individual run
                    auto start_time = std::chrono::high_resolution_clock::now();

                    // Run BISAM with optimized settings and current thread count
                    bisam::BisamResult result = bisam::estimate_model(
                        dataset.data,
                        0, 1, 2, // column indices: unit, time, y
                        MCMC_ITERATIONS,
                        MCMC_BURNIN,
                        "g",           // prior type
                        100.0,         // prior scale
                        0.001, 0.001,  // tolerances
                        1.0, 1.0, 1.0, // prior parameters
                        true,          // include indicators
                        false,
                        false,
                        false, // no fixed effects for performance testing
                        true,
                        true,                                         // IIS, SIS
                        2, 1, 2, 0,                                   // method parameters
                        100000, 1, 1,                                 // MCMC parameters
                        0.01, 0.01,                                   // alpha, lambda
                        1, 0.5, {1, 1},                               // additional parameters
                        bisam::ComputationStrategy::SPLIT_SEQUENTIAL, // Use optimized strategy
                        8                                             // Use current thread count
                    );

                    // End timing for this individual run
                    auto end_time = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    double execution_time_ms = duration.count() / 1000.0; // Convert to milliseconds

                    // std::cout << result.indicator_means.t() << std::endl;

                    // Store timing result
                    results_storage.add_timing_result(dataset.name, RUN_NAME, run + 1, execution_time_ms,
                                                      dataset.n, dataset.t, dataset.nx);

                    // Store estimation results
                    results_storage.add_estimation_results(dataset.name, RUN_NAME, run + 1, result, dataset);

                    // Always do validation but control printing
                    ValidationResults validation = BisamValidator::validate_results(
                        result, dataset, 0.5, ENABLE_VALIDATION_OUTPUT);
                    results_storage.add_validation_summary(dataset.name, RUN_NAME, run + 1, validation);

                    // Handle output based on setting
                    if (ENABLE_VALIDATION_OUTPUT) {
                        if (run == start_run - 1) {
                            // First run of this session
                            first_run_validation = validation;
                            BisamValidator::print_validation_summary(first_run_validation, dataset.name);
                            validation_done = true;
                        }
                    } else {
                        // Just show timing - no validation output
                        if (run == start_run - 1) {
                            // First run of this session
                            std::cout << std::fixed << std::setprecision(1) << execution_time_ms << "ms";
                            first_run_validation = validation;
                            validation_done      = true;
                        } else {
                            std::cout << ", " << std::fixed << std::setprecision(1) << execution_time_ms << "ms";
                        }
                    }

                    completed_runs++;

                    try {
                        results_storage.save_results();
                        results_storage.clear_memory(); // Clear accumulated results to prevent memory issues
                    } catch (const std::exception &e) {
                        std::cerr << "ERROR saving results: " << e.what() << std::endl;
                        return 4; // Exit with error code 4 for save failures
                    }
                } catch (const std::exception &e) {
                    std::cerr << "\nERROR in run " << (run + 1) << ": " << e.what() << std::endl;
                    std::cerr << "Aborting remaining runs for " << dataset.name << std::endl;
                    return 2; // Exit with error code 2 for estimation failures
                } catch (...) {
                    std::cerr << "\nUNKNOWN ERROR in run " << (run + 1) << std::endl;
                    std::cerr << "Aborting remaining runs for " << dataset.name << std::endl;
                    return 3; // Exit with error code 3 for unknown errors
                }
            }

            // Print brief summary for this dataset
            if (!ENABLE_VALIDATION_OUTPUT && validation_done) {
                std::cout << " | Beta: " << (first_run_validation.beta_max_error < 0.5 ? "OK" : "FAIL")
                        << " | Break: " << (first_run_validation.break_detected_correctly ? "OK" : "FAIL") << std::endl;
            }
        }

        if (!found_dataset) {
            std::cerr << "ERROR: No dataset found with name in test list" << std::endl;
            return 5; // Exit with error code 5 for dataset not found
        }

        // Verify we completed the expected number of runs
        int expected_runs = active_test_datasets.size() * (RUNS_PER_DATASET - start_run + 1); // Adjust for start_run
        if (completed_runs < expected_runs) {
            std::cerr << "ERROR: Incomplete execution - completed " << completed_runs
                    << "/" << expected_runs << " runs (starting from run " << start_run << ")" << std::endl;
            return 6; // Exit with error code 6 for incomplete execution
        }

        // Save results after each thread count iteration to ensure data is persisted

        // Final save (redundant but ensures everything is saved)
        std::cout << "\nFinal save of all results..." << std::endl;
        try {
            results_storage.save_results();
        } catch (const std::exception &e) {
            std::cerr << "ERROR in final save: " << e.what() << std::endl;
            return 7; // Exit with error code 7 for final save failure
        }

        std::cout << "Completed multithreading performance test (1-30 cores)" << std::endl;
        std::cout << "Successfully completed " << completed_runs << "/" << expected_runs << " runs" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "FATAL ERROR: " << e.what() << std::endl;
        return 1; // Exit with error code 1 for fatal errors
    } catch (...) {
        std::cerr << "FATAL UNKNOWN ERROR" << std::endl;
        return 1; // Exit with error code 1 for unknown fatal errors
    }

    return 0;
}
