# =============================================================================
# BISAM Simulation Dataset Generator
# =============================================================================

# Load required libraries
if (!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)

# Source the simulation function (assuming it's in the same directory)
# source("simulation_functions.R")  # Uncomment if function is in separate file

# Simulation function (included directly for completeness)
contr_sim_breaks <- function(
    n,
    t,
    nx,
    iis=TRUE,
    sis=TRUE,
    pos.outl, # count in total position of panel setup
    pos.step, # count in total position of panel setup
    const = FALSE,
    ife = FALSE,
    tfe = FALSE,
    outl.mean,
    step.mean,
    error.sd){
  
  require(mvtnorm)
  
  # Create individual and time indicators
  n_<- c(t(matrix(rep(1:n,t), ncol = rep(t,n)))) # unit index
  t_<- c(matrix(rep(1:t,n), ncol = n)) # time index
  
  # Create obs. matrices and corresponding betas
  if(nx > 0){
    X_ <- rmvnorm(n*t,rep(0,nx),diag(nx))
    X  <- X_ # original version will be required later
    
    b_ <- sample(x = -10:10,size = nx,replace = T)
    b  <- b_
  }else{
    X_ <- NULL
    X  <- NULL
    b_ <- NULL
    b  <- NULL
  }
  
  if(const) { # Intercept
    X <- cbind(X, "const" = rep(1,n*t))
    b0<- sample(x = -10:10,1)
    b <- c(b,b0)
  }
  if(ife) { # Individual fixed effects
    IFE <- kronecker(diag(n), matrix(1, t))
    if(const) {IFE <- IFE[,-1]}
    X <- cbind(X, IFE)
    bife<-sample(x = -10:10,size = ncol(IFE),replace = T)
    b <- c(b,bife)
  }
  if(tfe) { # Time fixed effects
    TFE <- kronecker(matrix(1, n), diag(t))
    if(const) {TFE <- TFE[,-1]}
    X <- cbind(X, TFE)
    btfe<-sample(x = -10:10,size = ncol(TFE),replace = T)
    b <- c(b,btfe)
  }
  if(tfe & ife & !const) {
    warning("Both time and unit fixed effects used.\n Dropping first indiv. FE dummy to avoid perfect colinearity")
    X <- X[, -(ncol(X_)+1)]
    b <- b[-length(b)]
  }
  # IIS
  I <- diag(n*t)
  tr_a <- rep(0,n*t)
  tr_a[pos.outl] <- 1
  a <- tr_a*outl.mean
  
  # SIS
  Z <- kronecker(diag(n),lower.tri(matrix(1,nrow=t,ncol=t),diag = T)[,-c(1,t*iis)])
  tr_g <- rep(0,n*(t-1-iis))
  tr_g[(pos.step%/%t)*(t-1-iis)+pos.step%%t-1] <- 1
  g <- tr_g*step.mean
  
  e <- rnorm(n*t,0,error.sd)
  
  if(is.null(X)){
    X = matrix(0,nrow=n*t)
    b = 0
  }
  
  y <- X%*%b + I%*%a + Z%*%g + e
  
  tr.ind <- cbind(which(t(matrix(a,ncol=t,byrow = T))!=0,arr.ind = T),a[a!=0])
  tr.stp <- cbind(which(t(cbind(rep(0,n),matrix(g,ncol=(t-1-iis),byrow = T)))!=0,arr.ind = T),g[g!=0])
  if(length(tr.ind)!=0){
    rownames(tr.ind) <- paste('iis',tr.ind[,2],tr.ind[,1],sep = '.')
  }
  if(length(tr.stp)!=0){
    rownames(tr.stp) <- paste('sis',tr.stp[,2],tr.stp[,1],sep = '.')
  }   
  tr.ind <- cbind(as.matrix(tr.ind),'index'=(tr.ind[,2]-1)*t+tr.ind[,1])
  tr.stp <- cbind(as.matrix(tr.stp),'index'=(tr.stp[,2]-1)*t+tr.stp[,1])
  treated_<- rbind(tr.ind,tr.stp)[,c(3,4),drop=F]
  treated<- cbind(treated_,'rel_net_eff' = (treated_[,1]+e[treated_[,2]])/error.sd)
  
  errors <- c(e,e)
  names(errors)<- c(paste('iis',n_,t_,sep = '.'),paste('sis',n_,t_,sep = '.'))
  # make data pretty
  data      <- cbind(n_, t_, y, X_)
  if(nx>0){
    colnames(data)<- c("n", "t", "y", paste0('x',1:nx))
  }else{
    colnames(data)<- c("n", "t", "y")
  }
  
  
  # valuate exact treatment timing
  colnames(treated)<-c('size','index','rel_net_eff')
  
  # structure output
  sim       <- list()
  sim$data  <- data
  sim$true.b<- b_
  sim$errors<- errors
  if(const){
    sim$true.const <- b0
  }
  if(ife){
    sim$true.ife <- bife
  }
  if(tfe){
    sim$true.tfe <- btfe
  }
  sim$tr.idx<- treated
  
  return(sim)
}

# Create/clear simulation data directory
sim_dir <- "simulation_data"
if (dir.exists(sim_dir)) {
  # Clear existing files
  file.remove(list.files(sim_dir, full.names = TRUE))
  cat("Cleared existing simulation data directory\n")
} else {
  dir.create(sim_dir)
  cat("Created simulation data directory\n")
}

# =============================================================================
# DATASET CONFIGURATIONS
# =============================================================================

# Define all possible dataset sizes - comment out ones you don't want to generate
dataset_configs <- list(
  # Very small (for quick testing/validation)
  list(n=3,  t=15, nx=2, name="tiny"),
  list(n=3,  t=20, nx=3, name="small_a"),
  list(n=5,  t=15, nx=2, name="small_b"),
  
  # Small (for detailed analysis)
  list(n=5,  t=20, nx=3, name="small_c"),
  list(n=5,  t=25, nx=3, name="small_d"),
  list(n=7,  t=20, nx=4, name="small_e"),
  list(n=8,  t=25, nx=3, name="small_f"),
  
  # Medium (for main comparisons)
  list(n=10, t=20, nx=4, name="med_a"),
  list(n=8,  t=30, nx=4, name="med_b"),
  list(n=10, t=30, nx=5, name="med_c"),
  list(n=12, t=35, nx=4, name="med_d"),
  list(n=15, t=30, nx=5, name="med_e"),
  
  # Large (realistic climate research sizes)
  list(n=15, t=40, nx=5, name="large_a"),  # 15 countries, 40 years
  list(n=20, t=35, nx=6, name="large_b"),  # 20 countries, 35 years
  list(n=25, t=40, nx=6, name="large_c"),  # 25 countries, 40 years
  list(n=30, t=45, nx=7, name="large_d"),  # 30 countries, 45 years
  
  # Very large (ambitious real-world sizes)
  list(n=35, t=50, nx=8, name="xlarge_a"), # 35 countries, 50 years
  list(n=40, t=50, nx=8, name="xlarge_b"), # 40 countries, 50 years
  list(n=50, t=55, nx=9, name="xlarge_c"), # 50 countries, 55 years
  
  # Extreme (stress testing)
  list(n=60, t=60, nx=10, name="extreme_a"), # 60 countries, 60 years
  list(n=75, t=65, nx=10, name="extreme_b"), # 75 countries, 65 years
  list(n=100, t=50, nx=12, name="extreme_c") # 100 countries, 50 years - EU + others
  
  # Uncomment below for even more extreme testing
  # list(n=150, t=70, nx=15, name="massive_a"), # All UN countries
  # list(n=200, t=80, nx=20, name="massive_b")  # Theoretical maximum
)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Fixed simulation parameters
set.seed(192837612)  # For reproducibility
const <- TRUE        # Include constant
ife <- FALSE        # No individual fixed effects (for performance)
tfe <- FALSE        # No time fixed effects (for performance)
iis <- TRUE         # Include indicator saturation
sis <- TRUE         # Include stepshift saturation
error.sd <- 1.0     # Error standard deviation

# Break parameters (matching your original setup)
outl.mean <- 0      # No outliers
step.mean <- 10     # Large step size

# =============================================================================
# DATA GENERATION FUNCTIONS
# =============================================================================

generate_dataset <- function(config) {
  n <- config$n
  t <- config$t
  nx <- config$nx
  name <- config$name
  
  cat(sprintf("Generating dataset: %s (n=%d, t=%d, nx=%d)\n", name, n, t, nx))
  
  # Generate break positions (matching your original approach)
  pos.outl <- c()  # No outliers (p.outl = 0.0)
  
  # Generate one step break per dataset (like your pos.step <- c(10))
  # Place it roughly in the middle of the time series
  step_time <- round(t * 0.5)  # Middle of time series
  step_unit <- sample(1:n, 1)  # Random unit
  pos.step <- (step_unit - 1) * t + step_time
  
  # Generate simulation data
  sim_data <- contr_sim_breaks(
    n = n,
    t = t,
    nx = nx,
    iis = iis,
    sis = sis,
    pos.outl = pos.outl,
    pos.step = pos.step,
    const = const,
    ife = ife,
    tfe = tfe,
    outl.mean = outl.mean,
    step.mean = step.mean,
    error.sd = error.sd
  )
  
  return(list(
    sim_data = sim_data,
    config = config,
    n = n, t = t, nx = nx, name = name
  ))
}

save_for_r <- function(dataset, sim_dir) {
  filename <- sprintf("%s/sim_%s_n%02d_t%02d_nx%02d.rds", 
                      sim_dir, dataset$name, dataset$n, dataset$t, dataset$nx)
  saveRDS(dataset, filename)
  cat(sprintf("  Saved R data: %s\n", basename(filename)))
}

save_for_cpp <- function(dataset, sim_dir) {
  n <- dataset$n
  t <- dataset$t
  nx <- dataset$nx
  name <- dataset$name
  sim <- dataset$sim_data
  
  filename <- sprintf("%s/sim_%s_n%02d_t%02d_nx%02d.h", 
                      sim_dir, name, n, t, nx)
  
  # Create C++ header file
  cat(sprintf("// Auto-generated simulation data: %s (n=%d, t=%d, nx=%d)\n", name, n, t, nx),
      file = filename)
  cat(sprintf("#ifndef SIM_%s_H\n", toupper(name)), file = filename, append = TRUE)
  cat(sprintf("#define SIM_%s_H\n\n", toupper(name)), file = filename, append = TRUE)
  cat("#include <armadillo>\n", file = filename, append = TRUE)
  cat("#include <vector>\n", file = filename, append = TRUE)
  cat("#include <string>\n\n", file = filename, append = TRUE)
  
  # Dataset metadata
  cat(sprintf("const int DATASET_%s_N = %d;\n", toupper(name), n), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_T = %d;\n", toupper(name), t), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_NX = %d;\n", toupper(name), nx), file = filename, append = TRUE)
  cat(sprintf("const std::string DATASET_%s_NAME = \"%s\";\n\n", toupper(name), name), file = filename, append = TRUE)
  
  # Main data matrix
  cat(sprintf("arma::mat DATASET_%s_DATA = {\n", toupper(name)), file = filename, append = TRUE)
  data <- sim$data
  for (i in 1:nrow(data)) {
    row_str <- paste0("        {", paste(sprintf("%.8f", data[i,]), collapse = ", "), "}")
    if (i < nrow(data)) row_str <- paste0(row_str, ",")
    cat(row_str, "\n", file = filename, append = TRUE)
  }
  cat("    };\n\n", file = filename, append = TRUE)
  
  # True betas
  if (!is.null(sim$true.b)) {
    cat(sprintf("arma::vec DATASET_%s_TRUE_BETA = {", toupper(name)), file = filename, append = TRUE)
    cat(paste(sprintf("%.8f", sim$true.b), collapse = ", "), file = filename, append = TRUE)
    cat("};\n\n", file = filename, append = TRUE)
  }
  
  # True constant (if present)
  if (!is.null(sim$true.const)) {
    cat(sprintf("const double DATASET_%s_TRUE_CONST = %.8f;\n", toupper(name), sim$true.const), 
        file = filename, append = TRUE)
    cat(sprintf("const bool DATASET_%s_HAS_TRUE_CONST = true;\n\n", toupper(name)), 
        file = filename, append = TRUE)
  } else {
    cat(sprintf("const double DATASET_%s_TRUE_CONST = 0.0;\n", toupper(name)), 
        file = filename, append = TRUE)
    cat(sprintf("const bool DATASET_%s_HAS_TRUE_CONST = false;\n\n", toupper(name)), 
        file = filename, append = TRUE)
  }
  
  # True break indices and magnitudes
  if (!is.null(sim$tr.idx) && nrow(sim$tr.idx) > 0) {
    cat(sprintf("arma::mat DATASET_%s_TRUE_BREAKS = {\n", toupper(name)), file = filename, append = TRUE)
    tr_idx <- sim$tr.idx
    for (i in 1:nrow(tr_idx)) {
      row_str <- paste0("        {", paste(sprintf("%.8f", tr_idx[i,]), collapse = ", "), "}")
      if (i < nrow(tr_idx)) row_str <- paste0(row_str, ",")
      cat(row_str, "\n", file = filename, append = TRUE)
    }
    cat("    };\n\n", file = filename, append = TRUE)
  } else {
    cat(sprintf("arma::mat DATASET_%s_TRUE_BREAKS;\n\n", toupper(name)), file = filename, append = TRUE)
  }
  
  # Simulation parameters for reference
  cat(sprintf("const double DATASET_%s_ERROR_SD = %.8f;\n", toupper(name), error.sd), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_STEP_MEAN = %.8f;\n", toupper(name), step.mean), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_OUTL_MEAN = %.8f;\n", toupper(name), outl.mean), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_CONST = %s;\n", toupper(name), ifelse(const, "true", "false")), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_IFE = %s;\n", toupper(name), ifelse(ife, "true", "false")), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_TFE = %s;\n\n", toupper(name), ifelse(tfe, "true", "false")), file = filename, append = TRUE)
  
  cat(sprintf("#endif // SIM_%s_H\n", toupper(name)), file = filename, append = TRUE)
  
  cat(sprintf("  Saved C++ header: %s\n", basename(filename)))
}

# =============================================================================
# GENERATE ALL DATASETS
# =============================================================================

cat("Starting dataset generation...\n\n")

# Store metadata for C++ automation
cpp_datasets <- data.frame(
  name = character(0),
  n = integer(0),
  t = integer(0),
  nx = integer(0),
  filename = character(0),
  stringsAsFactors = FALSE
)

# Generate each dataset
for (config in dataset_configs) {
  dataset <- generate_dataset(config)
  
  # Save in both formats
  save_for_r(dataset, sim_dir)
  save_for_cpp(dataset, sim_dir)
  
  # Add to C++ metadata
  cpp_datasets <- rbind(cpp_datasets, data.frame(
    name = config$name,
    n = config$n,
    t = config$t,
    nx = config$nx,
    filename = sprintf("sim_%s_n%02d_t%02d_nx%02d.h", config$name, config$n, config$t, config$nx),
    stringsAsFactors = FALSE
  ))
  
  cat("\n")
}

# =============================================================================
# CREATE C++ AUTOMATION FILES
# =============================================================================

# Create master header file for C++
cpp_master_file <- file.path(sim_dir, "simulation_datasets.h")
cat("// Auto-generated master header for all simulation datasets\n", file = cpp_master_file)
cat("#ifndef SIMULATION_DATASETS_H\n", file = cpp_master_file, append = TRUE)
cat("#define SIMULATION_DATASETS_H\n\n", file = cpp_master_file, append = TRUE)
cat("#include <vector>\n", file = cpp_master_file, append = TRUE)
cat("#include <string>\n", file = cpp_master_file, append = TRUE)
cat("#include <armadillo>\n\n", file = cpp_master_file, append = TRUE)

# Include all dataset headers
for (i in 1:nrow(cpp_datasets)) {
  cat(sprintf("#include \"%s\"\n", cpp_datasets$filename[i]), file = cpp_master_file, append = TRUE)
}

cat("\n// Dataset metadata structure\n", file = cpp_master_file, append = TRUE)
cat("struct DatasetInfo {\n", file = cpp_master_file, append = TRUE)
cat("    std::string name;\n", file = cpp_master_file, append = TRUE)
cat("    int n, t, nx;\n", file = cpp_master_file, append = TRUE)
cat("    arma::mat data;\n", file = cpp_master_file, append = TRUE)
cat("    arma::vec true_beta;\n", file = cpp_master_file, append = TRUE)
cat("    arma::mat true_breaks;\n", file = cpp_master_file, append = TRUE)
cat("    double true_const;\n", file = cpp_master_file, append = TRUE)
cat("    double error_sd;\n", file = cpp_master_file, append = TRUE)
cat("    double step_mean;\n", file = cpp_master_file, append = TRUE)
cat("    bool has_const;\n", file = cpp_master_file, append = TRUE)
cat("};\n\n", file = cpp_master_file, append = TRUE)

# Create dataset registry function
cat("// Get all available datasets\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<DatasetInfo> get_all_datasets() {\n", file = cpp_master_file, append = TRUE)
cat("    return {\n", file = cpp_master_file, append = TRUE)

for (i in 1:nrow(cpp_datasets)) {
  name_upper <- toupper(cpp_datasets$name[i])
  comma <- if(i < nrow(cpp_datasets)) "," else ""
  cat(sprintf("        {\"%s\", DATASET_%s_N, DATASET_%s_T, DATASET_%s_NX, DATASET_%s_DATA,\n", 
              cpp_datasets$name[i], name_upper, name_upper, name_upper, name_upper), 
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_TRUE_BETA, DATASET_%s_TRUE_BREAKS,\n", name_upper, name_upper),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_TRUE_CONST,\n", name_upper),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_ERROR_SD, DATASET_%s_STEP_MEAN, DATASET_%s_HAS_TRUE_CONST}%s\n", 
              name_upper, name_upper, name_upper, comma),
      file = cpp_master_file, append = TRUE)
}

cat("    };\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)
cat("#endif // SIMULATION_DATASETS_H\n", file = cpp_master_file, append = TRUE)

# Create R dataset registry
r_registry_file <- file.path(sim_dir, "dataset_registry.rds")
saveRDS(cpp_datasets, r_registry_file)

# =============================================================================
# CREATE EXAMPLE USAGE FILES
# =============================================================================

# Create example C++ usage
cpp_example <- file.path(sim_dir, "example_usage.cpp")
cat("// Example C++ usage of simulation datasets\n", file = cpp_example)
cat("#include \"simulation_datasets.h\"\n", file = cpp_example, append = TRUE)
cat("#include <iostream>\n", file = cpp_example, append = TRUE)
cat("#include <chrono>\n\n", file = cpp_example, append = TRUE)
cat("int main() {\n", file = cpp_example, append = TRUE)
cat("    auto datasets = get_all_datasets();\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    for (const auto& dataset : datasets) {\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"Running BISAM on dataset: \" << dataset.name \n", file = cpp_example, append = TRUE)
cat("                  << \" (n=\" << dataset.n << \", t=\" << dataset.t << \", nx=\" << dataset.nx << \")\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        // Access simulation truth\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"  True betas: \" << dataset.true_beta.t();\n", file = cpp_example, append = TRUE)
cat("        if (dataset.has_const) {\n", file = cpp_example, append = TRUE)
cat("            std::cout << \"  True constant: \" << dataset.true_const << std::endl;\n", file = cpp_example, append = TRUE)
cat("        }\n", file = cpp_example, append = TRUE)
cat("        if (dataset.true_breaks.n_rows > 0) {\n", file = cpp_example, append = TRUE)
cat("            std::cout << \"  True breaks: \" << dataset.true_breaks.n_rows << \" breaks\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        }\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        auto start = std::chrono::high_resolution_clock::now();\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        // YOUR BISAM CODE HERE\n", file = cpp_example, append = TRUE)
cat("        // arma::mat data_copy = dataset.data;  // Make non-const copy if needed\n", file = cpp_example, append = TRUE)
cat("        // bisam_result = run_bisam(data_copy);\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        // ACCURACY EVALUATION HERE\n", file = cpp_example, append = TRUE)
cat("        // Compare bisam_result.estimated_betas with dataset.true_beta\n", file = cpp_example, append = TRUE)
cat("        // Compare bisam_result.detected_breaks with dataset.true_breaks\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        auto end = std::chrono::high_resolution_clock::now();\n", file = cpp_example, append = TRUE)
cat("        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        std::cout << \"  Runtime: \" << duration.count() << \"ms\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"  \" << std::string(50, '-') << std::endl;\n", file = cpp_example, append = TRUE)
cat("    }\n", file = cpp_example, append = TRUE)
cat("    return 0;\n", file = cpp_example, append = TRUE)
cat("}\n", file = cpp_example, append = TRUE)

# Create example R usage
r_example <- file.path(sim_dir, "example_usage.R")
cat("# Example R usage of simulation datasets\n", file = r_example)
cat("library(here)\n\n", file = r_example, append = TRUE)
cat("# Load dataset registry\n", file = r_example, append = TRUE)
cat("registry <- readRDS('simulation_data/dataset_registry.rds')\n\n", file = r_example, append = TRUE)
cat("# Run BISAM on all datasets\n", file = r_example, append = TRUE)
cat("results <- list()\n\n", file = r_example, append = TRUE)
cat("for (i in 1:nrow(registry)) {\n", file = r_example, append = TRUE)
cat("  dataset_info <- registry[i, ]\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Load dataset\n", file = r_example, append = TRUE)
cat("  r_filename <- file.path('simulation_data', \n", file = r_example, append = TRUE)
cat("                          paste0('sim_', dataset_info$name, '_n', \n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$n), '_t',\n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$t), '_nx',\n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$nx), '.rds'))\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  dataset <- readRDS(r_filename)\n", file = r_example, append = TRUE)
cat("  sim <- dataset$sim_data\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  cat(sprintf('Running BISAM on dataset: %s (n=%d, t=%d, nx=%d)\\n', \n", file = r_example, append = TRUE)
cat("              dataset_info$name, dataset_info$n, dataset_info$t, dataset_info$nx))\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Show simulation truth\n", file = r_example, append = TRUE)
cat("  cat(sprintf('  True betas: [%s]\\n', paste(sim$true.b, collapse=', ')))\n", file = r_example, append = TRUE)
cat("  if (!is.null(sim$true.const)) {\n", file = r_example, append = TRUE)
cat("    cat(sprintf('  True constant: %.3f\\n', sim$true.const))\n", file = r_example, append = TRUE)
cat("  }\n", file = r_example, append = TRUE)
cat("  if (!is.null(sim$tr.idx) && nrow(sim$tr.idx) > 0) {\n", file = r_example, append = TRUE)
cat("    cat(sprintf('  True breaks: %d breaks\\n', nrow(sim$tr.idx)))\n", file = r_example, append = TRUE)
cat("  }\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Time the execution\n", file = r_example, append = TRUE)
cat("  start_time <- Sys.time()\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # YOUR BISAM CODE HERE\n", file = r_example, append = TRUE)
cat("  # bisam_result <- run_bisam(sim$data)\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # ACCURACY EVALUATION HERE\n", file = r_example, append = TRUE)
cat("  # Compare bisam_result$estimated_betas with sim$true.b\n", file = r_example, append = TRUE)
cat("  # Compare bisam_result$detected_breaks with sim$tr.idx\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  end_time <- Sys.time()\n", file = r_example, append = TRUE)
cat("  runtime <- as.numeric(difftime(end_time, start_time, units = 'secs'))\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  cat(sprintf('  Runtime: %.2f seconds\\n', runtime))\n", file = r_example, append = TRUE)
cat("  cat(paste(rep('-', 50), collapse=''), '\\n')\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Store results\n", file = r_example, append = TRUE)
cat("  results[[dataset_info$name]] <- list(\n", file = r_example, append = TRUE)
cat("    dataset_info = dataset_info,\n", file = r_example, append = TRUE)
cat("    runtime = runtime,\n", file = r_example, append = TRUE)
cat("    sim_truth = list(\n", file = r_example, append = TRUE)
cat("      true_beta = sim$true.b,\n", file = r_example, append = TRUE)
cat("      true_const = sim$true.const,\n", file = r_example, append = TRUE)
cat("      true_breaks = sim$tr.idx\n", file = r_example, append = TRUE)
cat("    )\n", file = r_example, append = TRUE)
cat("    # bisam_result = bisam_result\n", file = r_example, append = TRUE)
cat("  )\n", file = r_example, append = TRUE)
cat("}\n\n", file = r_example, append = TRUE)
cat("# Summary of results\n", file = r_example, append = TRUE)
cat("cat('\\nSummary of all runs:\\n')\n", file = r_example, append = TRUE)
cat("for (name in names(results)) {\n", file = r_example, append = TRUE)
cat("  r <- results[[name]]\n", file = r_example, append = TRUE)
cat("  cat(sprintf('%s: %.2fs\\n', name, r$runtime))\n", file = r_example, append = TRUE)
cat("}\n", file = r_example, append = TRUE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("DATASET GENERATION COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat(sprintf("Generated %d datasets in '%s/' directory\n", length(dataset_configs), sim_dir))
cat("\nDataset sizes:\n")
for (i in 1:nrow(cpp_datasets)) {
  cat(sprintf("  %s: n=%d, t=%d, nx=%d (total obs: %d)\n", 
              cpp_datasets$name[i], cpp_datasets$n[i], cpp_datasets$t[i], 
              cpp_datasets$nx[i], cpp_datasets$n[i] * cpp_datasets$t[i]))
}

cat("\nFiles created:\n")
cat("  - R data files: sim_*.rds\n")
cat("  - C++ headers: sim_*.h\n")
cat("  - Master C++ header: simulation_datasets.h\n")
cat("  - Dataset registry: dataset_registry.rds\n")
cat("  - Example usage: example_usage.cpp, example_usage.R\n")

cat("\nTo use in C++:\n")
cat("  #include \"simulation_data/simulation_datasets.h\"\n")
cat("  auto datasets = get_all_datasets();\n")

cat("\nTo use in R:\n")
cat("  registry <- readRDS('simulation_data/dataset_registry.rds')\n")
cat("  dataset <- readRDS('simulation_data/sim_[name]_n[XX]_t[XX]_nx[XX].rds')\n")

cat("\n", paste(rep("=", 60), collapse=""), "\n")