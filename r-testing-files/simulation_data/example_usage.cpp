// Example C++ usage of simulation datasets
#include "simulation_datasets.h"
#include <iostream>
#include <chrono>

int main() {
    auto datasets = get_all_datasets();
    
    for (const auto& dataset : datasets) {
        std::cout << "Running BISAM on dataset: " << dataset.name 
                  << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx << ")" << std::endl;
        
        // Access simulation truth
        std::cout << "  True betas: " << dataset.true_beta.t();
        if (dataset.has_const) {
            std::cout << "  True constant: " << dataset.true_const << std::endl;
        }
        if (dataset.true_breaks.n_rows > 0) {
            std::cout << "  True breaks: " << dataset.true_breaks.n_rows << " breaks" << std::endl;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // YOUR BISAM CODE HERE
        // arma::mat data_copy = dataset.data;  // Make non-const copy if needed
        // bisam_result = run_bisam(data_copy);
        
        // ACCURACY EVALUATION HERE
        // Compare bisam_result.estimated_betas with dataset.true_beta
        // Compare bisam_result.detected_breaks with dataset.true_breaks
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "  Runtime: " << duration.count() << "ms" << std::endl;
        std::cout << "  " << std::string(50, '-') << std::endl;
    }
    return 0;
}
