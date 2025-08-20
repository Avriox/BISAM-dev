// Example C++ usage of simulation datasets for root finding algorithm comparison
#include "simulation_datasets.h"
#include <iostream>
#include <chrono>

int main() {
    // Test different step sizes with same root finding algorithm
    auto stepsize_datasets = get_stepsize_datasets();
    std::cout << "Testing different step sizes:" << std::endl;
    for (const auto& dataset : stepsize_datasets) {
        std::cout << "Dataset: " << dataset.name << " (step_mean=" << dataset.step_mean << ")" << std::endl;
        std::cout << "  True breaks: " << dataset.num_breaks << " breaks" << std::endl;
        
        // YOUR BISAM CODE HERE with specific root finding algorithm
        // Compare results across different step sizes
    }
    
    // Test different time series lengths
    auto timelength_datasets = get_timelength_datasets();
    std::cout << "\nTesting different time series lengths:" << std::endl;
    for (const auto& dataset : timelength_datasets) {
        std::cout << "Dataset: " << dataset.name << " (t=" << dataset.t << ")" << std::endl;
        std::cout << "  True breaks: " << dataset.num_breaks << " breaks" << std::endl;
        
        // YOUR BISAM CODE HERE with specific root finding algorithm
        // Compare results across different time series lengths
    }
    
    return 0;
}
