# =============================================================================
# BISAM Simulation Dataset Generator - Fixed for Root Finding Algorithm Comparison
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
# DATASET CONFIGURATIONS - Updated for Root Finding Algorithm Comparison
# =============================================================================

dataset_configs <- list(
  # Very small (for quick testing/validation)
  list(n=3,  t=15, nx=2, step_mean=5.00, name="tiny"),
  list(n=3,  t=20, nx=3, step_mean=5.00, name="small_a"),
  list(n=5,  t=15, nx=2, step_mean=5.00, name="small_b"),
  
  # Small (for detailed analysis)
  list(n=5,  t=20, nx=3, step_mean=5.00, name="small_c"),
  list(n=5,  t=25, nx=3, step_mean=5.00, name="small_d"),
  list(n=7,  t=20, nx=4, step_mean=5.00, name="small_e"),
  list(n=8,  t=25, nx=3, step_mean=5.00, name="small_f"),
  
  # Medium (for main comparisons)
  list(n=10, t=20, nx=4, step_mean=5.00, name="med_a"),
  list(n=8,  t=30, nx=4, step_mean=5.00, name="med_b"),
  list(n=10, t=30, nx=5, step_mean=5.00, name="med_c"),
  list(n=12, t=35, nx=4, step_mean=5.00, name="med_d"),
  list(n=15, t=30, nx=5, step_mean=5.00, name="med_e"),
  
  # Large (realistic climate research sizes)
  list(n=15, t=40, nx=5, step_mean=5.00, name="large_a"),  # 15 countries, 40 years
  list(n=20, t=35, nx=6, step_mean=5.00, name="large_b"),  # 20 countries, 35 years
  list(n=25, t=40, nx=6, step_mean=5.00, name="large_c"),  # 25 countries, 40 years
  list(n=30, t=45, nx=7, step_mean=5.00, name="large_d"),  # 30 countries, 45 years
  
  # Very large (ambitious real-world sizes)
  list(n=35, t=50, nx=8, step_mean=5.00, name="xlarge_a"), # 35 countries, 50 years
  list(n=40, t=50, nx=8, step_mean=5.00, name="xlarge_b"), # 40 countries, 50 years
  list(n=50, t=55, nx=9, step_mean=5.00, name="xlarge_c"), # 50 countries, 55 years
  
  # Extreme (stress testing)
  list(n=60, t=60, nx=10, step_mean=5.00, name="extreme_a"), # 60 countries, 60 years
  list(n=75, t=65, nx=10, step_mean=5.00, name="extreme_b"), # 75 countries, 65 years
  list(n=100, t=50, nx=12, step_mean=5.00, name="extreme_c"), # 100 countries, 50 years - EU + others
  
  # =============================================================================
  # SET 1: BREAK SIZE COMPARISON (n=10, t=30, varying step sizes)
  # Fixed: n=10, t=30, nx=3
  # Variable: step_mean (break size)
  # ~50% of units have breaks in the middle third
  # =============================================================================
  
  list(n=10, t=30, nx=3, step_mean=0.50, name="rootfind_stepsize_050"),
  list(n=10, t=30, nx=3, step_mean=0.75, name="rootfind_stepsize_075"),
  list(n=10, t=30, nx=3, step_mean=1.00, name="rootfind_stepsize_100"),
  list(n=10, t=30, nx=3, step_mean=1.50, name="rootfind_stepsize_150"),
  list(n=10, t=30, nx=3, step_mean=3.00, name="rootfind_stepsize_300"),
  list(n=10, t=30, nx=3, step_mean=5.00, name="rootfind_stepsize_500"),
  
  # =============================================================================
  # SET 2: TIME SERIES LENGTH COMPARISON (n=10, fixed step_mean=5, varying t)
  # Fixed: n=10, step_mean=5, nx=3
  # Variable: t (time series length)
  # ~50% of units have breaks in the middle third
  # =============================================================================
  
  list(n=10, t=10, nx=3, step_mean=5.0, name="rootfind_timelength_t010"),
  list(n=10, t=15, nx=3, step_mean=5.0, name="rootfind_timelength_t015"),
  list(n=10, t=20, nx=3, step_mean=5.0, name="rootfind_timelength_t020"),
  list(n=10, t=25, nx=3, step_mean=5.0, name="rootfind_timelength_t025"),
  list(n=10, t=30, nx=3, step_mean=5.0, name="rootfind_timelength_t030"),
  list(n=10, t=40, nx=3, step_mean=5.0, name="rootfind_timelength_t040"),
  list(n=10, t=50, nx=3, step_mean=5.0, name="rootfind_timelength_t050"),
  list(n=10, t=60, nx=3, step_mean=5.0, name="rootfind_timelength_t060"),
  list(n=10, t=70, nx=3, step_mean=5.0, name="rootfind_timelength_t070"),
  list(n=10, t=80, nx=3, step_mean=5.0, name="rootfind_timelength_t080"),
  list(n=10, t=90, nx=3, step_mean=5.0, name="rootfind_timelength_t090"),
  list(n=10, t=100, nx=3, step_mean=5.0, name="rootfind_timelength_t100"),
  list(n=10, t=110, nx=3, step_mean=5.0, name="rootfind_timelength_t110"),
  list(n=10, t=120, nx=3, step_mean=5.0, name="rootfind_timelength_t120")
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

# Break parameters - these will be overridden per dataset
outl.mean <- 0      # No outliers (as requested)
# step.mean will be set per dataset configuration

# =============================================================================
# DATA GENERATION FUNCTIONS
# =============================================================================

generate_dataset <- function(config) {
  n <- config$n
  t <- config$t
  nx <- config$nx
  name <- config$name
  step_mean_config <- config$step_mean  # Get step_mean from config
  
  cat(sprintf("Generating dataset: %s (n=%d, t=%d, nx=%d, step_mean=%.2f)\n", 
              name, n, t, nx, step_mean_config))
  
  # Generate break positions
  pos.outl <- c()  # No outliers (as requested)
  
  # Calculate number of breaks as approximately half the units (minimum 1, maximum n-1)
  num_breaks <- max(1, min(n-1, round(n * 0.5)))
  
  cat(sprintf("  Creating %d breaks for %d units\n", num_breaks, n))
  
  # Generate breaks across different units in the middle third of time series
  set.seed(as.numeric(charToRaw(name)[1]) + n + t)  # Reproducible but different per dataset
  
  # Select random units for breaks
  break_units <- sample(1:n, num_breaks, replace = FALSE)
  
  # For each selected unit, place break in middle third of time series at different positions
  middle_start <- max(2, round(t * 0.33))  # Start of middle third, but not too early
  middle_end <- min(t-1, round(t * 0.67))  # End of middle third, but not too late
  
  # Ensure we have valid range
  if (middle_start >= middle_end) {
    middle_start <- max(2, round(t * 0.4))
    middle_end <- min(t-1, round(t * 0.6))
  }
  if (middle_start >= middle_end) {
    middle_start <- 2
    middle_end <- t-1
  }
  
  pos.step <- c()
  for (i in 1:num_breaks) {
    unit <- break_units[i]
    # Distribute breaks across the middle third
    if (num_breaks > 1) {
      spread_factor <- (i-1) / (num_breaks-1)
      break_time <- round(middle_start + spread_factor * (middle_end - middle_start))
    } else {
      break_time <- round((middle_start + middle_end) / 2)
    }
    
    # Add some randomness to avoid identical timing
    jitter_range <- max(1, round((middle_end - middle_start) * 0.1))
    break_time <- break_time + sample(-jitter_range:jitter_range, 1)
    break_time <- max(middle_start, min(middle_end, break_time))
    
    pos.step <- c(pos.step, (unit - 1) * t + break_time)
  }
  
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
    step.mean = step_mean_config,  # Use step_mean from config
    error.sd = error.sd
  )
  
  # Add break information to output for verification
  break_info <- data.frame(
    unit = break_units,
    time = (pos.step - 1) %% t + 1,
    position = pos.step,
    step_mean = rep(step_mean_config, num_breaks),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  Breaks placed at: %s\n", 
              paste(sprintf("Unit %d Time %d", break_info$unit, break_info$time), collapse=", ")))
  
  return(list(
    sim_data = sim_data,
    config = config,
    break_info = break_info,
    n = n, t = t, nx = nx, name = name,
    num_breaks = num_breaks
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
  step_mean_used <- dataset$config$step_mean
  
  filename <- sprintf("%s/sim_%s_n%02d_t%02d_nx%02d.h", 
                      sim_dir, name, n, t, nx)
  
  # Create C++ header file
  cat(sprintf("// Auto-generated simulation data: %s (n=%d, t=%d, nx=%d, step_mean=%.2f)\n", 
              name, n, t, nx, step_mean_used),
      file = filename)
  cat(sprintf("#ifndef SIM_%s_H\n", toupper(gsub("[^A-Z0-9]", "_", toupper(name)))), file = filename, append = TRUE)
  cat(sprintf("#define SIM_%s_H\n\n", toupper(gsub("[^A-Z0-9]", "_", toupper(name)))), file = filename, append = TRUE)
  cat("#include <armadillo>\n", file = filename, append = TRUE)
  cat("#include <vector>\n", file = filename, append = TRUE)
  cat("#include <string>\n\n", file = filename, append = TRUE)
  
  # Dataset metadata
  name_cpp <- toupper(gsub("[^A-Z0-9]", "_", toupper(name)))
  cat(sprintf("const int DATASET_%s_N = %d;\n", name_cpp, n), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_T = %d;\n", name_cpp, t), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_NX = %d;\n", name_cpp, nx), file = filename, append = TRUE)
  cat(sprintf("const std::string DATASET_%s_NAME = \"%s\";\n", name_cpp, name), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_STEP_MEAN_USED = %.8f;\n\n", name_cpp, step_mean_used), file = filename, append = TRUE)
  
  # Main data matrix
  cat(sprintf("arma::mat DATASET_%s_DATA = {\n", name_cpp), file = filename, append = TRUE)
  data <- sim$data
  for (i in 1:nrow(data)) {
    row_str <- paste0("        {", paste(sprintf("%.8f", data[i,]), collapse = ", "), "}")
    if (i < nrow(data)) row_str <- paste0(row_str, ",")
    cat(row_str, "\n", file = filename, append = TRUE)
  }
  cat("    };\n\n", file = filename, append = TRUE)
  
  # True betas
  if (!is.null(sim$true.b)) {
    cat(sprintf("arma::vec DATASET_%s_TRUE_BETA = {", name_cpp), file = filename, append = TRUE)
    cat(paste(sprintf("%.8f", sim$true.b), collapse = ", "), file = filename, append = TRUE)
    cat("};\n\n", file = filename, append = TRUE)
  }
  
  # True constant (if present)
  if (!is.null(sim$true.const)) {
    cat(sprintf("const double DATASET_%s_TRUE_CONST = %.8f;\n", name_cpp, sim$true.const), 
        file = filename, append = TRUE)
    cat(sprintf("const bool DATASET_%s_HAS_TRUE_CONST = true;\n\n", name_cpp), 
        file = filename, append = TRUE)
  } else {
    cat(sprintf("const double DATASET_%s_TRUE_CONST = 0.0;\n", name_cpp), 
        file = filename, append = TRUE)
    cat(sprintf("const bool DATASET_%s_HAS_TRUE_CONST = false;\n\n", name_cpp), 
        file = filename, append = TRUE)
  }
  
  # True break indices and magnitudes
  if (!is.null(sim$tr.idx) && nrow(sim$tr.idx) > 0) {
    cat(sprintf("arma::mat DATASET_%s_TRUE_BREAKS = {\n", name_cpp), file = filename, append = TRUE)
    tr_idx <- sim$tr.idx
    for (i in 1:nrow(tr_idx)) {
      row_str <- paste0("        {", paste(sprintf("%.8f", tr_idx[i,]), collapse = ", "), "}")
      if (i < nrow(tr_idx)) row_str <- paste0(row_str, ",")
      cat(row_str, "\n", file = filename, append = TRUE)
    }
    cat("    };\n\n", file = filename, append = TRUE)
  } else {
    cat(sprintf("arma::mat DATASET_%s_TRUE_BREAKS;\n\n", name_cpp), file = filename, append = TRUE)
  }
  
  # Break information for verification
  if (!is.null(dataset$break_info)) {
    cat(sprintf("// Break placement info: %d breaks\n", nrow(dataset$break_info)), file = filename, append = TRUE)
    for (i in 1:nrow(dataset$break_info)) {
      cat(sprintf("// Unit %d, Time %d, Position %d\n", 
                  dataset$break_info$unit[i], dataset$break_info$time[i], dataset$break_info$position[i]), 
          file = filename, append = TRUE)
    }
    cat("\n", file = filename, append = TRUE)
  }
  
  # Simulation parameters for reference
  cat(sprintf("const double DATASET_%s_ERROR_SD = %.8f;\n", name_cpp, error.sd), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_STEP_MEAN = %.8f;\n", name_cpp, step_mean_used), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_OUTL_MEAN = %.8f;\n", name_cpp, outl.mean), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_CONST = %s;\n", name_cpp, ifelse(const, "true", "false")), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_IFE = %s;\n", name_cpp, ifelse(ife, "true", "false")), file = filename, append = TRUE)
  cat(sprintf("const bool DATASET_%s_TFE = %s;\n", name_cpp, ifelse(tfe, "true", "false")), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_NUM_BREAKS = %d;\n\n", name_cpp, 
              if(!is.null(sim$tr.idx)) nrow(sim$tr.idx) else 0), file = filename, append = TRUE)
  
  cat(sprintf("#endif // SIM_%s_H\n", name_cpp), file = filename, append = TRUE)
  
  cat(sprintf("  Saved C++ header: %s\n", basename(filename)))
}

# =============================================================================
# GENERATE ALL DATASETS
# =============================================================================

cat("Starting dataset generation for root finding algorithm comparison...\n\n")

# Store metadata for C++ automation
cpp_datasets <- data.frame(
  name = character(0),
  n = integer(0),
  t = integer(0),
  nx = integer(0),
  step_mean = numeric(0),
  filename = character(0),
  num_breaks = integer(0),
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
    step_mean = config$step_mean,
    filename = sprintf("sim_%s_n%02d_t%02d_nx%02d.h", config$name, config$n, config$t, config$nx),
    num_breaks = dataset$num_breaks,
    stringsAsFactors = FALSE
  ))
  
  cat("\n")
}

# =============================================================================
# CREATE C++ AUTOMATION FILES
# =============================================================================

# Create master header file for C++
cpp_master_file <- file.path(sim_dir, "simulation_datasets.h")
cat("// Auto-generated master header for all simulation datasets - Root Finding Algorithm Comparison\n", file = cpp_master_file)
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
cat("    double step_mean;\n", file = cpp_master_file, append = TRUE)
cat("    arma::mat data;\n", file = cpp_master_file, append = TRUE)
cat("    arma::vec true_beta;\n", file = cpp_master_file, append = TRUE)
cat("    arma::mat true_breaks;\n", file = cpp_master_file, append = TRUE)
cat("    double true_const;\n", file = cpp_master_file, append = TRUE)
cat("    double error_sd;\n", file = cpp_master_file, append = TRUE)
cat("    bool has_const;\n", file = cpp_master_file, append = TRUE)
cat("    int num_breaks;\n", file = cpp_master_file, append = TRUE)
cat("};\n\n", file = cpp_master_file, append = TRUE)

# Create dataset registry function
cat("// Get all available datasets\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<DatasetInfo> get_all_datasets() {\n", file = cpp_master_file, append = TRUE)
cat("    return {\n", file = cpp_master_file, append = TRUE)

for (i in 1:nrow(cpp_datasets)) {
  name_cpp <- toupper(gsub("[^A-Z0-9]", "_", toupper(cpp_datasets$name[i])))
  comma <- if(i < nrow(cpp_datasets)) "," else ""
  cat(sprintf("        {\"%s\", DATASET_%s_N, DATASET_%s_T, DATASET_%s_NX,\n", 
              cpp_datasets$name[i], name_cpp, name_cpp, name_cpp), 
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_STEP_MEAN_USED, DATASET_%s_DATA,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_TRUE_BETA, DATASET_%s_TRUE_BREAKS,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_TRUE_CONST, DATASET_%s_ERROR_SD,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("         DATASET_%s_HAS_TRUE_CONST, DATASET_%s_NUM_BREAKS}%s\n", 
              name_cpp, name_cpp, comma),
      file = cpp_master_file, append = TRUE)
}

cat("    };\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)

# Add convenience functions
cat("// Get datasets by type\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<DatasetInfo> get_stepsize_datasets() {\n", file = cpp_master_file, append = TRUE)
cat("    auto all = get_all_datasets();\n", file = cpp_master_file, append = TRUE)
cat("    std::vector<DatasetInfo> result;\n", file = cpp_master_file, append = TRUE)
cat("    for (const auto& ds : all) {\n", file = cpp_master_file, append = TRUE)
cat("        if (ds.name.find(\"stepsize\") != std::string::npos) {\n", file = cpp_master_file, append = TRUE)
cat("            result.push_back(ds);\n", file = cpp_master_file, append = TRUE)
cat("        }\n", file = cpp_master_file, append = TRUE)
cat("    }\n", file = cpp_master_file, append = TRUE)
cat("    return result;\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)

cat("inline std::vector<DatasetInfo> get_timelength_datasets() {\n", file = cpp_master_file, append = TRUE)
cat("    auto all = get_all_datasets();\n", file = cpp_master_file, append = TRUE)
cat("    std::vector<DatasetInfo> result;\n", file = cpp_master_file, append = TRUE)
cat("    for (const auto& ds : all) {\n", file = cpp_master_file, append = TRUE)
cat("        if (ds.name.find(\"timelength\") != std::string::npos) {\n", file = cpp_master_file, append = TRUE)
cat("            result.push_back(ds);\n", file = cpp_master_file, append = TRUE)
cat("        }\n", file = cpp_master_file, append = TRUE)
cat("    }\n", file = cpp_master_file, append = TRUE)
cat("    return result;\n", file = cpp_master_file, append = TRUE)
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
cat("// Example C++ usage of simulation datasets for root finding algorithm comparison\n", file = cpp_example)
cat("#include \"simulation_datasets.h\"\n", file = cpp_example, append = TRUE)
cat("#include <iostream>\n", file = cpp_example, append = TRUE)
cat("#include <chrono>\n\n", file = cpp_example, append = TRUE)
cat("int main() {\n", file = cpp_example, append = TRUE)
cat("    // Test different step sizes with same root finding algorithm\n", file = cpp_example, append = TRUE)
cat("    auto stepsize_datasets = get_stepsize_datasets();\n", file = cpp_example, append = TRUE)
cat("    std::cout << \"Testing different step sizes:\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("    for (const auto& dataset : stepsize_datasets) {\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"Dataset: \" << dataset.name << \" (step_mean=\" << dataset.step_mean << \")\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"  True breaks: \" << dataset.num_breaks << \" breaks\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        // YOUR BISAM CODE HERE with specific root finding algorithm\n", file = cpp_example, append = TRUE)
cat("        // Compare results across different step sizes\n", file = cpp_example, append = TRUE)
cat("    }\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    // Test different time series lengths\n", file = cpp_example, append = TRUE)
cat("    auto timelength_datasets = get_timelength_datasets();\n", file = cpp_example, append = TRUE)
cat("    std::cout << \"\\nTesting different time series lengths:\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("    for (const auto& dataset : timelength_datasets) {\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"Dataset: \" << dataset.name << \" (t=\" << dataset.t << \")\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        std::cout << \"  True breaks: \" << dataset.num_breaks << \" breaks\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("        \n", file = cpp_example, append = TRUE)
cat("        // YOUR BISAM CODE HERE with specific root finding algorithm\n", file = cpp_example, append = TRUE)
cat("        // Compare results across different time series lengths\n", file = cpp_example, append = TRUE)
cat("    }\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    return 0;\n", file = cpp_example, append = TRUE)
cat("}\n", file = cpp_example, append = TRUE)

# Create example R usage
r_example <- file.path(sim_dir, "example_usage.R")
cat("# Example R usage for root finding algorithm comparison\n", file = r_example)
cat("library(here)\n\n", file = r_example, append = TRUE)
cat("# Load dataset registry\n", file = r_example, append = TRUE)
cat("registry <- readRDS('simulation_data/dataset_registry.rds')\n\n", file = r_example, append = TRUE)
cat("# Split datasets by type\n", file = r_example, append = TRUE)
cat("stepsize_datasets <- registry[grepl('stepsize', registry$name), ]\n", file = r_example, append = TRUE)
cat("timelength_datasets <- registry[grepl('timelength', registry$name), ]\n\n", file = r_example, append = TRUE)
cat("# Test different step sizes\n", file = r_example, append = TRUE)
cat("cat('Testing different step sizes:\\n')\n", file = r_example, append = TRUE)
cat("stepsize_results <- list()\n\n", file = r_example, append = TRUE)
cat("for (i in 1:nrow(stepsize_datasets)) {\n", file = r_example, append = TRUE)
cat("  dataset_info <- stepsize_datasets[i, ]\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Load dataset\n", file = r_example, append = TRUE)
cat("  r_filename <- file.path('simulation_data', \n", file = r_example, append = TRUE)
cat("                          paste0('sim_', dataset_info$name, '_n', \n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$n), '_t',\n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$t), '_nx',\n", file = r_example, append = TRUE)
cat("                                sprintf('%02d', dataset_info$nx), '.rds'))\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  dataset <- readRDS(r_filename)\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  cat(sprintf('Dataset: %s (step_mean=%.2f, %d breaks)\\n', \n", file = r_example, append = TRUE)
cat("              dataset_info$name, dataset_info$step_mean, dataset_info$num_breaks))\n", file = r_example, append = TRUE)
cat("  cat(sprintf('  Break info: %s\\n', \n", file = r_example, append = TRUE)
cat("              paste(sprintf('Unit %d Time %d', dataset$break_info$unit, dataset$break_info$time), collapse=', ')))\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # YOUR BISAM CODE HERE with specific root finding algorithm\n", file = r_example, append = TRUE)
cat("  # Compare results across different step sizes\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  stepsize_results[[dataset_info$name]] <- list(\n", file = r_example, append = TRUE)
cat("    dataset_info = dataset_info,\n", file = r_example, append = TRUE)
cat("    break_info = dataset$break_info\n", file = r_example, append = TRUE)
cat("    # Add your results here\n", file = r_example, append = TRUE)
cat("  )\n", file = r_example, append = TRUE)
cat("}\n\n", file = r_example, append = TRUE)
cat("# Test different time series lengths\n", file = r_example, append = TRUE)
cat("cat('\\nTesting different time series lengths:\\n')\n", file = r_example, append = TRUE)
cat("timelength_results <- list()\n\n", file = r_example, append = TRUE)
cat("for (i in 1:nrow(timelength_datasets)) {\n", file = r_example, append = TRUE)
cat("  dataset_info <- timelength_datasets[i, ]\n", file = r_example, append = TRUE)
cat("  \n", file = r_example, append = TRUE)
cat("  # Load dataset (similar to above)\n", file = r_example, append = TRUE)
cat("  # YOUR BISAM CODE HERE with specific root finding algorithm\n", file = r_example, append = TRUE)
cat("  # Compare results across different time series lengths\n", file = r_example, append = TRUE)
cat("}\n", file = r_example, append = TRUE)

# =============================================================================
# SUMMARY
# =============================================================================

# Get counts for stepsize and timelength datasets
stepsize_datasets <- cpp_datasets[grepl('stepsize', cpp_datasets$name), ]
timelength_datasets <- cpp_datasets[grepl('timelength', cpp_datasets$name), ]

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("DATASET GENERATION COMPLETE - ROOT FINDING ALGORITHM COMPARISON\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat(sprintf("Generated %d datasets in '%s/' directory\n", length(dataset_configs), sim_dir))
cat("\nDataset breakdown:\n")
cat(sprintf("  Step size comparison: %d datasets (step_mean: %.2f to %.2f)\n", 
            nrow(stepsize_datasets), min(stepsize_datasets$step_mean), max(stepsize_datasets$step_mean)))
cat(sprintf("  Time length comparison: %d datasets (t: %d to %d)\n", 
            nrow(timelength_datasets), min(timelength_datasets$t), max(timelength_datasets$t)))
cat(sprintf("  Other datasets: %d datasets for validation and stress testing\n", 
            nrow(cpp_datasets) - nrow(stepsize_datasets) - nrow(timelength_datasets)))

cat("\nKey characteristics:\n")
cat("  - Adaptive number of breaks (approximately 50% of units)\n")
cat("  - Breaks placed in middle third of time series at different positions\n")
cat("  - No outliers (outl.mean = 0)\n")
cat("  - All datasets have unique names for easy identification\n")

cat("\nFiles created:\n")
cat("  - R data files: sim_*.rds\n")
cat("  - C++ headers: sim_*.h\n")
cat("  - Master C++ header: simulation_datasets.h\n")
cat("  - Dataset registry: dataset_registry.rds\n")
cat("  - Example usage: example_usage.cpp, example_usage.R\n")

cat("\nTo use in C++:\n")
cat("  auto stepsize_datasets = get_stepsize_datasets();\n")
cat("  auto timelength_datasets = get_timelength_datasets();\n")

cat("\n", paste(rep("=", 60), collapse=""), "\n")