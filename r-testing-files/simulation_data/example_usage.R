# Example R usage of simulation datasets
library(here)

# Load dataset registry
registry <- readRDS('simulation_data/dataset_registry.rds')

# Run BISAM on all datasets
results <- list()

for (i in 1:nrow(registry)) {
  dataset_info <- registry[i, ]
  
  # Load dataset
  r_filename <- file.path('simulation_data', 
                          paste0('sim_', dataset_info$name, '_n', 
                                sprintf('%02d', dataset_info$n), '_t',
                                sprintf('%02d', dataset_info$t), '_nx',
                                sprintf('%02d', dataset_info$nx), '.rds'))
  
  dataset <- readRDS(r_filename)
  sim <- dataset$sim_data
  
  cat(sprintf('Running BISAM on dataset: %s (n=%d, t=%d, nx=%d)\n', 
              dataset_info$name, dataset_info$n, dataset_info$t, dataset_info$nx))
  
  # Show simulation truth
  cat(sprintf('  True betas: [%s]\n', paste(sim$true.b, collapse=', ')))
  if (!is.null(sim$true.const)) {
    cat(sprintf('  True constant: %.3f\n', sim$true.const))
  }
  if (!is.null(sim$tr.idx) && nrow(sim$tr.idx) > 0) {
    cat(sprintf('  True breaks: %d breaks\n', nrow(sim$tr.idx)))
  }
  
  # Time the execution
  start_time <- Sys.time()
  
  # YOUR BISAM CODE HERE
  # bisam_result <- run_bisam(sim$data)
  
  # ACCURACY EVALUATION HERE
  # Compare bisam_result$estimated_betas with sim$true.b
  # Compare bisam_result$detected_breaks with sim$tr.idx
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  
  cat(sprintf('  Runtime: %.2f seconds\n', runtime))
  cat(paste(rep('-', 50), collapse=''), '\n')
  
  # Store results
  results[[dataset_info$name]] <- list(
    dataset_info = dataset_info,
    runtime = runtime,
    sim_truth = list(
      true_beta = sim$true.b,
      true_const = sim$true.const,
      true_breaks = sim$tr.idx
    )
    # bisam_result = bisam_result
  )
}

# Summary of results
cat('\nSummary of all runs:\n')
for (name in names(results)) {
  r <- results[[name]]
  cat(sprintf('%s: %.2fs\n', name, r$runtime))
}
