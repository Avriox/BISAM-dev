#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>

#include "include/biasm_model.h"
#include "r-testing-files/simulation_data/simulation_datasets.h"
#include "validation.h"

// ========================================================================
// CONFIGURATION SECTION
// ========================================================================
const bool ENABLE_VALIDATION_OUTPUT = false; // Set to false to disable all validation printing
const int RUNS_PER_DATASET          = 20;    // Number of times to run each dataset for timing
const bool SHOW_DETAILED_RESULTS    = false; // Set to false for cleaner output

// MCMC Settings
const int MCMC_ITERATIONS = 5000; // Total MCMC iterations
const int MCMC_BURNIN     = 500;  // Burn-in period

int main() {
    bisam::FunctionTimer timer;

    // Get all available datasets
    auto datasets = get_all_datasets();

    // Select which datasets to test (comment/uncomment as needed)
    std::vector<std::string> test_datasets = {
        "tiny",    // Very small for quick testing
        "small_a", // Small for detailed analysis
        "small_b", // Another small dataset
        // "med_a",     // Medium size
        // "large_a"    // Large size - takes longer
    };

    std::cout << "BISAM Performance and Validation Testing" << std::endl;
    std::cout << "Testing " << test_datasets.size() << " datasets with " << RUNS_PER_DATASET << " runs each" <<
            std::endl;
    std::cout << "MCMC Settings: " << MCMC_ITERATIONS << " iterations, " << MCMC_BURNIN << " burn-in" << std::endl;
    std::cout << "Validation output: " << (ENABLE_VALIDATION_OUTPUT ? "ENABLED" : "DISABLED") << std::endl;

    for (const auto &dataset: datasets) {
        // Skip datasets not in test list
        if (std::find(test_datasets.begin(), test_datasets.end(), dataset.name) == test_datasets.end()) {
            continue;
        }

        // Print dataset information once per dataset
        if (ENABLE_VALIDATION_OUTPUT) {
            BisamValidator::print_dataset_info(dataset);
        } else {
            std::cout << "\n" << std::string(50, '=') << std::endl;
            std::cout << "DATASET: " << dataset.name
                    << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx << ")" << std::endl;
            std::cout << std::string(50, '=') << std::endl;
        }

        // Variables to store validation results from first run
        ValidationResults first_run_validation;
        bool validation_done = false;

        // Run the dataset multiple times for timing
        for (int run = 0; run < RUNS_PER_DATASET; run++) {
            std::cout << "\nRun " << (run + 1) << "/" << RUNS_PER_DATASET << " for dataset " << dataset.name;
            if (!ENABLE_VALIDATION_OUTPUT) {
                std::cout << "...";
            }
            std::cout << std::endl;

            // Start timing
            timer.start_section("BISAM_" + dataset.name);


            // Run BISAM with optimized settings
            bisam::BisamResult result = bisam::estimate_model(
                dataset.data,
                0,
                1,
                2,
                // column indices: unit, time, y
                MCMC_ITERATIONS,
                MCMC_BURNIN,
                // use configuration values
                "g",
                // prior type
                100.0,
                // prior scale
                0.001,
                0.001,
                // tolerances
                1.0,
                1.0,
                1.0,
                // prior parameters
                true,
                // include indicators
                false,
                false,
                false,
                // no fixed effects for performance testing
                true,
                true,
                // IIS, SIS
                2,
                1,
                2,
                0,
                // method parameters
                100000,
                1,
                1,
                // MCMC parameters
                0.01,
                0.01,
                // alpha, lambda
                1,
                0.5,
                {1, 1},
                // additional parameters
                bisam::ComputationStrategy::SPLIT_SEQUENTIAL // Use optimized strategy
            );

            timer.end_section("BISAM_" + dataset.name);

            // Only do validation on the first run and only if enabled
            if (ENABLE_VALIDATION_OUTPUT && run == 0) {
                first_run_validation = BisamValidator::validate_results(result, dataset, 0.5);
                BisamValidator::print_validation_summary(first_run_validation, dataset.name);

                if (SHOW_DETAILED_RESULTS) {
                    BisamValidator::print_detailed_results(result, dataset);
                }
                validation_done = true;
            } else if (!ENABLE_VALIDATION_OUTPUT && run == 0) {
                // At least validate silently to check correctness
                first_run_validation = BisamValidator::validate_results(result, dataset, 0.5);
                validation_done      = true;

                // Print brief summary without detailed output
                std::cout << "  Beta accuracy: "
                        << (first_run_validation.beta_max_error < 0.5 ? "GOOD" : "NEEDS ATTENTION")
                        << " (max error: " << std::fixed << std::setprecision(4)
                        << first_run_validation.beta_max_error << ")" << std::endl;
                std::cout << "  Break detection: "
                        << (first_run_validation.break_detected_correctly ? "CORRECT" : "INCORRECT")
                        << " (prob: " << first_run_validation.detected_break_probability << ")" << std::endl;
            }

            if (!ENABLE_VALIDATION_OUTPUT) {
                std::cout << "  Run " << (run + 1) << " completed successfully" << std::endl;
            }
        }

        // Print summary for this dataset
        if (validation_done && !ENABLE_VALIDATION_OUTPUT) {
            std::cout << "\nDataset " << dataset.name << " summary:" << std::endl;
            std::cout << "  Completed " << RUNS_PER_DATASET << " runs successfully" << std::endl;
            std::cout << "  Beta estimation: "
                    << (first_run_validation.beta_max_error < 0.5 ? "GOOD" : "NEEDS ATTENTION") << std::endl;
            std::cout << "  Break detection: "
                    << (first_run_validation.break_detected_correctly ? "CORRECT" : "INCORRECT") << std::endl;
        }
    }

    // Print timing summary
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "PERFORMANCE SUMMARY (" << RUNS_PER_DATASET << " runs per dataset)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    timer.print_section_summary();

    // Print configuration summary
    std::cout << "\nConfiguration used:" << std::endl;
    std::cout << "  MCMC: " << MCMC_ITERATIONS << " iterations, " << MCMC_BURNIN << " burn-in" << std::endl;
    std::cout << "  Runs per dataset: " << RUNS_PER_DATASET << std::endl;
    std::cout << "  Validation output: " << (ENABLE_VALIDATION_OUTPUT ? "ENABLED" : "DISABLED") << std::endl;

    return 0;
}
