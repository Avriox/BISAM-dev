# Example R usage for root finding algorithm comparison
library(here)

# Load dataset registry
registry <- readRDS('simulation_data/dataset_registry.rds')

# Split datasets by type
stepsize_datasets <- registry[grepl('stepsize', registry$name), ]
timelength_datasets <- registry[grepl('timelength', registry$name), ]

# Test different step sizes
cat('Testing different step sizes:\n')
stepsize_results <- list()

for (i in 1:nrow(stepsize_datasets)) {
  dataset_info <- stepsize_datasets[i, ]
  
  # Load dataset
  r_filename <- file.path('simulation_data', 
                          paste0('sim_', dataset_info$name, '_n', 
                                sprintf('%02d', dataset_info$n), '_t',
                                sprintf('%02d', dataset_info$t), '_nx',
                                sprintf('%02d', dataset_info$nx), '.rds'))
  
  dataset <- readRDS(r_filename)
  
  cat(sprintf('Dataset: %s (step_mean=%.2f, %d breaks)\n', 
              dataset_info$name, dataset_info$step_mean, dataset_info$num_breaks))
  cat(sprintf('  Break info: %s\n', 
              paste(sprintf('Unit %d Time %d', dataset$break_info$unit, dataset$break_info$time), collapse=', ')))
  
  # YOUR BISAM CODE HERE with specific root finding algorithm
  # Compare results across different step sizes
  
  stepsize_results[[dataset_info$name]] <- list(
    dataset_info = dataset_info,
    break_info = dataset$break_info
    # Add your results here
  )
}

# Test different time series lengths
cat('\nTesting different time series lengths:\n')
timelength_results <- list()

for (i in 1:nrow(timelength_datasets)) {
  dataset_info <- timelength_datasets[i, ]
  
  # Load dataset (similar to above)
  # YOUR BISAM CODE HERE with specific root finding algorithm
  # Compare results across different time series lengths
}
