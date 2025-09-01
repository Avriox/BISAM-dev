# =============================================================================
# BISAM Simulation Dataset Generator - Fixed for Root Finding Algorithm Comparison
# =============================================================================

# Load required libraries
if (!require(mvtnorm, quietly = TRUE)) install.packages("mvtnorm")
if (!require(ggplot2, quietly = TRUE)) install.packages("ggplot2")
if (!require(dplyr, quietly = TRUE)) install.packages("dplyr")
if (!require(tidyr, quietly = TRUE)) install.packages("tidyr")

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# CONTROL VARIABLE FOR NUMBER OF DATASETS PER CONFIGURATION
# =============================================================================
NUM_DATASETS_PER_CONFIG <- 10  # Change this to control how many datasets per step size

# =============================================================================
# Controlled Simulation
# =============================================================================

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
  
    print( pos.outl)
    print(pos.step)

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


# =============================================================================
# Create/clear simulation data directory
# =============================================================================

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

# Fixed breaks (default) for rootfind_* datasets (you can override per dataset)
# Example fixed breaks (n=10, works best for t >= 19)
fixed_breaks_default <- data.frame(
  unit = c(8, 2, 7, 4, 6),
  time = c(11, 11, 16, 18, 19)
  # Positions will be computed as (unit - 1) * t + time per dataset
)

dataset_configs <- list(
  # =============================================================================
  # SET 1: BREAK SIZE COMPARISON (n=10, t=30, varying step sizes)
  # Fixed: n=10, t=30, nx=3; Variable: step_mean (break size)
  # =============================================================================

  # list(n=10, t=30, nx=3, step_mean=0.50, name="rootfind_stepsize_050", breaks=fixed_breaks_default),
  # list(n=10, t=30, nx=3, step_mean=0.75, name="rootfind_stepsize_075", breaks=fixed_breaks_default),
  # list(n=10, t=30, nx=3, step_mean=1.00, name="rootfind_stepsize_100", breaks=fixed_breaks_default),
  # list(n=10, t=30, nx=3, step_mean=1.50, name="rootfind_stepsize_150", breaks=fixed_breaks_default),
  list(n=8, t=30, nx=3, step_mean=3.00, name="rootfind_stepsize_300", breaks=fixed_breaks_default)
  # list(n=10, t=30, nx=3, step_mean=5.00, name="rootfind_stepsize_500", breaks=fixed_breaks_default)
)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Fixed simulation parameters
const <- FALSE        # Include constant
ife <- FALSE          # No individual fixed effects (for performance)
tfe <- FALSE          # No time fixed effects (for performance)
iis <- TRUE           # Include indicator saturation
sis <- TRUE           # Include stepshift saturation
error.sd <- 1.0       # Error standard deviation

# Break parameters - these will be overridden per dataset
outl.mean <- 0        # No outliers (as requested)

# =============================================================================
# DATA GENERATION FUNCTIONS
# =============================================================================

resolve_breaks <- function(br, n, t) {
  # Returns a 1-based integer vector of positions in 1..(n*t)
  nt <- n * t
  if (is.null(br)) return(NULL)

  # If numeric vector: assume pos_step
  if (is.numeric(br)) {
    pos <- as.integer(br)
    pos <- pos[pos >= 1 & pos <= nt]
    return(unique(pos))
  }

  # If function: call with (n, t)
  if (is.function(br)) {
    pos <- br(n, t)
    pos <- as.integer(pos)
    pos <- pos[pos >= 1 & pos <= nt]
    return(unique(pos))
  }

  # If list with pos_step or unit/time
  if (is.list(br)) {
    if (!is.null(br$pos_step)) {
      pos <- as.integer(br$pos_step)
      pos <- pos[pos >= 1 & pos <= nt]
      return(unique(pos))
    }
    if (!is.null(br$unit) && !is.null(br$time)) {
      unit <- as.integer(br$unit)
      time <- as.integer(br$time)
      stopifnot(length(unit) == length(time))
      pos <- (unit - 1L) * t + time
      pos <- pos[pos >= 1 & pos <= nt]
      return(unique(pos))
    }
  }

  # If data.frame or matrix with unit/time columns
  if (is.data.frame(br) || is.matrix(br)) {
    if (all(c("unit", "time") %in% colnames(as.data.frame(br)))) {
      df <- as.data.frame(br)
      unit <- as.integer(df$unit)
      time <- as.integer(df$time)
      pos <- (unit - 1L) * t + time
      pos <- pos[pos >= 1 & pos <= nt]
      return(unique(pos))
    }
  }

  warning("Unrecognized breaks specification; ignoring and falling back to random.")
  return(NULL)
}

generate_dataset <- function(config, dataset_index) {
  n <- config$n
  t <- config$t
  nx <- config$nx
  name <- config$name
  step_mean_config <- config$step_mean

  # Create unique name for this dataset instance
  dataset_name <- sprintf("%s_%03d", name, dataset_index)
  
  cat(sprintf("  [%03d/%03d] Generating: %s (n=%d, t=%d, nx=%d, step_mean=%.2f)\n",
              dataset_index, NUM_DATASETS_PER_CONFIG, dataset_name, n, t, nx, step_mean_config))

  # Set unique seed for this dataset (but breaks remain the same)
  set.seed(192837615 + dataset_index * 1000 + as.numeric(charToRaw(name)[1]))

  # Generate break positions (same for all datasets in this config)
  pos.outl <- integer(0) # No outliers (as requested)

  # Use fixed breaks from config
  pos.step <- NULL
  if (!is.null(config$breaks)) {
    pos.step <- resolve_breaks(config$breaks, n, t)
  }

  # Generate simulation data
  sim_data <- contr_sim_breaks(
    n = n, t = t, nx = nx,
    iis = iis, sis = sis,
    pos.outl = pos.outl,
    pos.step = pos.step,
    const = const, ife = ife, tfe = tfe,
    outl.mean = outl.mean,
    step.mean = step_mean_config,
    error.sd = error.sd
  )

  # Add break information to output for verification (keeps full data index)
  break_info <- NULL
  if (!is.null(pos.step) && length(pos.step) > 0) {
    break_info <- data.frame(
      unit = ((pos.step - 1) %/% t) + 1,
      time = ((pos.step - 1) %% t) + 1,
      position = pos.step, # full data index (1-based)
      step_mean = rep(step_mean_config, length(pos.step)),
      stringsAsFactors = FALSE
    )
  }

  return(list(
    sim_data = sim_data,
    config = config,
    break_info = break_info,
    n = n, t = t, nx = nx, 
    name = dataset_name,  # Use the unique name with index
    num_breaks = if (!is.null(break_info)) nrow(break_info) else 0
  ))
}

save_for_cpp <- function(dataset, sim_dir) {
  n <- dataset$n
  t <- dataset$t
  nx <- dataset$nx
  name <- dataset$name
  sim <- dataset$sim_data
  step_mean_used <- dataset$config$step_mean

  name_cpp <- toupper(gsub("[^A-Z0-9]", "", name))
  guard <- paste0("SIM_", name_cpp, "_H")

  filename <- sprintf("%s/sim_%s_n%02d_t%02d_nx%02d.h",
                      sim_dir, name, n, t, nx)

  # Create C++ header file
  cat(sprintf("// Auto-generated simulation data: %s (n=%d, t=%d, nx=%d, step_mean=%.2f)\n",
              name, n, t, nx, step_mean_used),
      file = filename)
  cat(sprintf("#ifndef %s\n", guard), file = filename, append = TRUE)
  cat(sprintf("#define %s\n\n", guard), file = filename, append = TRUE)
  cat("#include <armadillo>\n", file = filename, append = TRUE)
  cat("#include <vector>\n", file = filename, append = TRUE)
  cat("#include <string>\n\n", file = filename, append = TRUE)

  # Dataset metadata
  cat(sprintf("const int DATASET_%s_N = %d;\n", name_cpp, n), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_T = %d;\n", name_cpp, t), file = filename, append = TRUE)
  cat(sprintf("const int DATASET_%s_NX = %d;\n", name_cpp, nx), file = filename, append = TRUE)
  cat(sprintf("const std::string DATASET_%s_NAME = \"%s\";\n", name_cpp, name), file = filename, append = TRUE)
  cat(sprintf("const double DATASET_%s_STEP_MEAN_USED = %.8f;\n\n", name_cpp, step_mean_used), file = filename, append = TRUE)

  # Main data matrix
  cat(sprintf("arma::mat DATASET_%s_DATA = {\n", name_cpp), file = filename, append = TRUE)
  data <- sim$data
  for (i in 1:nrow(data)) {
    row_str <- paste0("  {", paste(sprintf("%.8f", as.numeric(data[i, ])), collapse = ", "), "}")
    if (i < nrow(data)) row_str <- paste0(row_str, ",")
    cat(row_str, "\n", file = filename, append = TRUE)
  }
  cat("};\n\n", file = filename, append = TRUE)

  # True betas
  if (!is.null(sim$true.b)) {
    cat(sprintf("arma::vec DATASET_%s_TRUE_BETA = {", name_cpp), file = filename, append = TRUE)
    cat(paste(sprintf("%.8f", sim$true.b), collapse = ", "), file = filename, append = TRUE)
    cat("};\n\n", file = filename, append = TRUE)
  } else {
    cat(sprintf("arma::vec DATASET_%s_TRUE_BETA;\n\n", name_cpp), file = filename, append = TRUE)
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
      row_str <- paste0("  {", paste(sprintf("%.8f", as.numeric(tr_idx[i, ])), collapse = ", "), "}")
      if (i < nrow(tr_idx)) row_str <- paste0(row_str, ",")
      cat(row_str, "\n", file = filename, append = TRUE)
    }
    cat("};\n\n", file = filename, append = TRUE)
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
              if (!is.null(sim$tr.idx)) nrow(sim$tr.idx) else 0), file = filename, append = TRUE)

  cat(sprintf("#endif // %s\n", guard), file = filename, append = TRUE)
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

plot_combined_timeseries <- function(dataset, title_suffix = "") {
  # Extract data and reshape for plotting
  sim_data <- dataset$sim_data
  data_df <- as.data.frame(sim_data$data)
  
  # Add break information if available
  break_info <- dataset$break_info
  
  # Create base title
  base_title <- sprintf("%s (step_mean=%.2f)", 
                       dataset$config$name, 
                       dataset$config$step_mean)
  if (title_suffix != "") {
    base_title <- paste(base_title, title_suffix)
  }
  
  # Determine which units have breaks
  units_with_breaks <- c()
  if (!is.null(break_info) && nrow(break_info) > 0) {
    units_with_breaks <- unique(break_info$unit)
  }
  
  # Add break status to data
  data_df$has_break <- data_df$n %in% units_with_breaks
  data_df$unit_type <- ifelse(data_df$has_break, 
                              paste("Unit", data_df$n, "(with break)"), 
                              paste("Unit", data_df$n, "(no break)"))
  data_df$line_type <- ifelse(data_df$has_break, "solid", "dotted")
  
  # Create color palette that distinguishes break vs no-break units
  all_units <- unique(data_df$n)
  n_units <- length(all_units)
  
  # Use different color families for break vs no-break
  colors_with_breaks <- rainbow(sum(all_units %in% units_with_breaks), start = 0, end = 0.6)
  colors_without_breaks <- rainbow(sum(!all_units %in% units_with_breaks), start = 0.7, end = 1)
  
  # Create color mapping
  unit_colors <- c()
  color_idx_with <- 1
  color_idx_without <- 1
  
  for (unit in sort(all_units)) {
    if (unit %in% units_with_breaks) {
      unit_colors[paste("Unit", unit, "(with break)")] <- colors_with_breaks[color_idx_with]
      color_idx_with <- color_idx_with + 1
    } else {
      unit_colors[paste("Unit", unit, "(no break)")] <- colors_without_breaks[color_idx_without]
      color_idx_without <- color_idx_without + 1
    }
  }
  
  # Combined plot: Time series by unit with break status in legend
  p1 <- ggplot(data_df, aes(x = t, y = y, group = n, color = unit_type, linetype = line_type)) +
    geom_line(alpha = 0.8, size = 0.7) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = unit_colors, name = "Unit (Break Status)") +
    scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted"),
                         guide = "none") +  # Hide linetype legend since it's in the color legend
    labs(title = paste("Combined Time Series by Unit -", base_title),
         x = "Time", y = "y") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.key.width = unit(1.5, "cm"))  # Make legend lines longer to see line types
  
  # Add break markers if available
  if (!is.null(break_info) && nrow(break_info) > 0) {
    break_points <- data_df[data_df$n %in% break_info$unit & data_df$t %in% break_info$time, ]
    if (nrow(break_points) > 0) {
      p1 <- p1 + 
        geom_vline(data = break_info, aes(xintercept = time), 
                   color = "red", linetype = "dashed", alpha = 0.5) +
        geom_point(data = break_points, aes(x = t, y = y), 
                   color = "red", size = 2, shape = 17)
    }
  }
  
  return(p1)
}

plot_concatenated_timeseries <- function(dataset, title_suffix = "") {
  # Extract data and reshape for plotting
  sim_data <- dataset$sim_data
  data_df <- as.data.frame(sim_data$data)
  
  # Add break information if available
  break_info <- dataset$break_info
  
  # Create base title
  base_title <- sprintf("%s (step_mean=%.2f)", 
                       dataset$config$name, 
                       dataset$config$step_mean)
  if (title_suffix != "") {
    base_title <- paste(base_title, title_suffix)
  }
  
  # Create concatenated time series data
  # Sort by unit then by time to ensure proper ordering
  data_sorted <- data_df[order(data_df$n, data_df$t), ]
  
  # Create new time index that concatenates all units
  n_units <- length(unique(data_df$n))
  t_periods <- length(unique(data_df$t))
  
  # Create continuous time index
  data_concat <- data_sorted
  data_concat$time_concat <- 1:nrow(data_concat)
  
  # Create unit boundary markers for visual separation
  unit_boundaries <- seq(t_periods, by = t_periods, length.out = n_units - 1)
  
  # Concatenated plot: All time series end-to-end
  p2 <- ggplot(data_concat, aes(x = time_concat, y = y)) +
    geom_line(alpha = 0.8, color = "darkblue") +
    geom_point(size = 0.3, alpha = 0.6, color = "darkblue") +
    labs(title = paste("Concatenated Time Series (All Units End-to-End) -", base_title),
         x = "Concatenated Time Index", y = "y") +
    theme_minimal()
  
  # Add unit boundary lines
  if (length(unit_boundaries) > 0) {
    p2 <- p2 + geom_vline(xintercept = unit_boundaries + 0.5, 
                          color = "gray", linetype = "dotted", alpha = 0.7)
  }
  
  # Add break markers if available (adjusted for concatenated time)
  if (!is.null(break_info) && nrow(break_info) > 0) {
    # Calculate concatenated positions for breaks
    break_concat_positions <- c()
    break_y_values <- c()
    
    for (i in 1:nrow(break_info)) {
      unit_id <- break_info$unit[i]
      break_time <- break_info$time[i]
      
      # Find concatenated position: (unit-1)*t_periods + time
      concat_pos <- (unit_id - 1) * t_periods + break_time
      break_concat_positions <- c(break_concat_positions, concat_pos)
      
      # Get y value at break point
      break_y <- data_concat$y[data_concat$time_concat == concat_pos]
      if (length(break_y) > 0) {
        break_y_values <- c(break_y_values, break_y[1])
      }
    }
    
    if (length(break_concat_positions) > 0) {
      break_data <- data.frame(
        x = break_concat_positions,
        y = if(length(break_y_values) == length(break_concat_positions)) break_y_values else rep(mean(data_concat$y), length(break_concat_positions))
      )
      
      p2 <- p2 + 
        geom_vline(data = break_data, aes(xintercept = x), 
                   color = "red", linetype = "dashed", alpha = 0.7) +
        geom_point(data = break_data, aes(x = x, y = y), 
                   color = "red", size = 2, shape = 17)
    }
  }
  
  # Add unit labels on x-axis
  unit_midpoints <- seq(t_periods/2, by = t_periods, length.out = n_units)
  unit_labels <- paste("Unit", 1:n_units)
  
  p2 <- p2 + 
    scale_x_continuous(
      breaks = unit_midpoints,
      labels = unit_labels,
      minor_breaks = NULL
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p2)
}

# =============================================================================
# GENERATE ALL DATASETS
# =============================================================================

cat("Starting dataset generation for root finding algorithm comparison...\n")
cat(sprintf("Will generate %d datasets for each configuration\n\n", NUM_DATASETS_PER_CONFIG))

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

# Store some sample datasets for plotting
sample_datasets <- list()

# Generate each dataset configuration
for (config in dataset_configs) {
  cat(sprintf("\nGenerating %d datasets for config: %s (step_mean=%.2f)\n", 
              NUM_DATASETS_PER_CONFIG, config$name, config$step_mean))
  
  for (dataset_index in 1:NUM_DATASETS_PER_CONFIG) {
    dataset <- generate_dataset(config, dataset_index)
    
    # Save in C++ format only
    save_for_cpp(dataset, sim_dir)
    
    # Store first dataset of each config for plotting
    if (dataset_index == 4) {
      sample_datasets <- append(sample_datasets, list(dataset), length(sample_datasets))
    }
    
    # Add to C++ metadata
    cpp_datasets <- rbind(cpp_datasets, data.frame(
      name = dataset$name,
      n = dataset$n,
      t = dataset$t,
      nx = dataset$nx,
      step_mean = dataset$config$step_mean,
      filename = sprintf("sim_%s_n%02d_t%02d_nx%02d.h", dataset$name, dataset$n, dataset$t, dataset$nx),
      num_breaks = dataset$num_breaks,
      stringsAsFactors = FALSE
    ))
  }
}

# =============================================================================
# CREATE VISUALIZATIONS
# =============================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n", sep = "")
cat("CREATING VISUALIZATIONS\n")
cat(paste(rep("=", 60), collapse=""), "\n", sep = "")

# Create plots for sample datasets
cat("Creating time series plots for sample datasets...\n")

# Plot individual datasets (all sample datasets)
for (i in 1:length(sample_datasets)) {
  plot(sample_datasets[[i]]$sim_data$data[,3])
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
cat("  std::string name;\n", file = cpp_master_file, append = TRUE)
cat("  int n, t, nx;\n", file = cpp_master_file, append = TRUE)
cat("  double step_mean;\n", file = cpp_master_file, append = TRUE)
cat("  arma::mat data;\n", file = cpp_master_file, append = TRUE)
cat("  arma::vec true_beta;\n", file = cpp_master_file, append = TRUE)
cat("  arma::mat true_breaks;\n", file = cpp_master_file, append = TRUE)
cat("  double true_const;\n", file = cpp_master_file, append = TRUE)
cat("  double error_sd;\n", file = cpp_master_file, append = TRUE)
cat("  bool has_const;\n", file = cpp_master_file, append = TRUE)
cat("  int num_breaks;\n", file = cpp_master_file, append = TRUE)
cat("};\n\n", file = cpp_master_file, append = TRUE)

# Create dataset registry function
cat("// Get all available datasets\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<DatasetInfo> get_all_datasets() {\n", file = cpp_master_file, append = TRUE)
cat("  return {\n", file = cpp_master_file, append = TRUE)

for (i in 1:nrow(cpp_datasets)) {
  name_cpp <- toupper(gsub("[^A-Z0-9]", "", cpp_datasets$name[i]))
  comma <- if (i < nrow(cpp_datasets)) "," else ""
  cat(sprintf("    {\"%s\", DATASET_%s_N, DATASET_%s_T, DATASET_%s_NX,\n",
              cpp_datasets$name[i], name_cpp, name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("     DATASET_%s_STEP_MEAN_USED, DATASET_%s_DATA,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("     DATASET_%s_TRUE_BETA, DATASET_%s_TRUE_BREAKS,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("     DATASET_%s_TRUE_CONST, DATASET_%s_ERROR_SD,\n", name_cpp, name_cpp),
      file = cpp_master_file, append = TRUE)
  cat(sprintf("     DATASET_%s_HAS_TRUE_CONST, DATASET_%s_NUM_BREAKS}%s\n",
              name_cpp, name_cpp, comma),
      file = cpp_master_file, append = TRUE)
}

cat("  };\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)

# Add convenience functions
cat("// Get datasets by step size\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<DatasetInfo> get_datasets_by_step_mean(double step_mean, double tolerance = 0.01) {\n", file = cpp_master_file, append = TRUE)
cat("  auto all = get_all_datasets();\n", file = cpp_master_file, append = TRUE)
cat("  std::vector<DatasetInfo> result;\n", file = cpp_master_file, append = TRUE)
cat("  for (const auto& ds : all) {\n", file = cpp_master_file, append = TRUE)
cat("    if (std::abs(ds.step_mean - step_mean) < tolerance) {\n", file = cpp_master_file, append = TRUE)
cat("      result.push_back(ds);\n", file = cpp_master_file, append = TRUE)
cat("    }\n", file = cpp_master_file, append = TRUE)
cat("  }\n", file = cpp_master_file, append = TRUE)
cat("  return result;\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)

cat("// Get unique step sizes\n", file = cpp_master_file, append = TRUE)
cat("inline std::vector<double> get_unique_step_sizes() {\n", file = cpp_master_file, append = TRUE)
cat("  std::vector<double> step_sizes = {", file = cpp_master_file, append = TRUE)
unique_steps <- unique(cpp_datasets$step_mean)
cat(paste(sprintf("%.2f", sort(unique_steps)), collapse = ", "), file = cpp_master_file, append = TRUE)
cat("};\n", file = cpp_master_file, append = TRUE)
cat("  return step_sizes;\n", file = cpp_master_file, append = TRUE)
cat("}\n\n", file = cpp_master_file, append = TRUE)

cat("#endif // SIMULATION_DATASETS_H\n", file = cpp_master_file, append = TRUE)

# =============================================================================
# CREATE EXAMPLE USAGE FILE
# =============================================================================

# Create example C++ usage
cpp_example <- file.path(sim_dir, "example_usage.cpp")
cat("// Example C++ usage of simulation datasets for root finding algorithm comparison\n", file = cpp_example)
cat("#include \"simulation_datasets.h\"\n", file = cpp_example, append = TRUE)
cat("#include <iostream>\n", file = cpp_example, append = TRUE)
cat("#include <chrono>\n\n", file = cpp_example, append = TRUE)
cat("int main() {\n", file = cpp_example, append = TRUE)
cat("  // Get all unique step sizes\n", file = cpp_example, append = TRUE)
cat("  auto step_sizes = get_unique_step_sizes();\n", file = cpp_example, append = TRUE)
cat("  \n", file = cpp_example, append = TRUE)
cat("  // Process each step size\n", file = cpp_example, append = TRUE)
cat("  for (double step_mean : step_sizes) {\n", file = cpp_example, append = TRUE)
cat("    std::cout << \"\\nProcessing step size: \" << step_mean << std::endl;\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    // Get all datasets for this step size\n", file = cpp_example, append = TRUE)
cat("    auto datasets = get_datasets_by_step_mean(step_mean);\n", file = cpp_example, append = TRUE)
cat("    std::cout << \"Found \" << datasets.size() << \" datasets\" << std::endl;\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    // Process each dataset\n", file = cpp_example, append = TRUE)
cat("    for (const auto& dataset : datasets) {\n", file = cpp_example, append = TRUE)
cat("      std::cout << \"  Processing: \" << dataset.name << std::endl;\n", file = cpp_example, append = TRUE)
cat("      \n", file = cpp_example, append = TRUE)
cat("      // YOUR ROOT FINDING ALGORITHM CODE HERE\n", file = cpp_example, append = TRUE)
cat("      // Example: Run your BISAM implementation with different root finding methods\n", file = cpp_example, append = TRUE)
cat("      // - Secant method\n", file = cpp_example, append = TRUE)
cat("      // - Newton-Raphson\n", file = cpp_example, append = TRUE)
cat("      // - Brent's method\n", file = cpp_example, append = TRUE)
cat("      // etc.\n", file = cpp_example, append = TRUE)
cat("      \n", file = cpp_example, append = TRUE)
cat("      // Access data:\n", file = cpp_example, append = TRUE)
cat("      // arma::mat data = dataset.data;\n", file = cpp_example, append = TRUE)
cat("      // arma::vec true_beta = dataset.true_beta;\n", file = cpp_example, append = TRUE)
cat("      // arma::mat true_breaks = dataset.true_breaks;\n", file = cpp_example, append = TRUE)
cat("      // int n = dataset.n;\n", file = cpp_example, append = TRUE)
cat("      // int t = dataset.t;\n", file = cpp_example, append = TRUE)
cat("      // int nx = dataset.nx;\n", file = cpp_example, append = TRUE)
cat("    }\n", file = cpp_example, append = TRUE)
cat("    \n", file = cpp_example, append = TRUE)
cat("    // Aggregate and compare results for this step size\n", file = cpp_example, append = TRUE)
cat("    // ...\n", file = cpp_example, append = TRUE)
cat("  }\n", file = cpp_example, append = TRUE)
cat("  \n", file = cpp_example, append = TRUE)
cat("  return 0;\n", file = cpp_example, append = TRUE)
cat("}\n", file = cpp_example, append = TRUE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n", sep = "")
cat("DATASET GENERATION COMPLETE - ROOT FINDING ALGORITHM COMPARISON\n")
cat(paste(rep("=", 60), collapse=""), "\n", sep = "")
cat(sprintf("Generated %d total datasets in '%s/' directory\n", nrow(cpp_datasets), sim_dir))
cat(sprintf("  - %d unique configurations\n", length(dataset_configs)))
cat(sprintf("  - %d datasets per configuration\n", NUM_DATASETS_PER_CONFIG))

cat("\nStep size breakdown:\n")
for (step in sort(unique(cpp_datasets$step_mean))) {
  count <- sum(cpp_datasets$step_mean == step)
  cat(sprintf("  Step size %.2f: %d datasets\n", step, count))
}

cat("\nAll datasets have:\n")
cat("  - Fixed break positions (same across all datasets for same config)\n")
cat("  - Unique data realizations (different random seeds)\n")
cat("  - n=10, t=30, nx=3\n")
cat("  - 5 breaks at: Unit 8 Time 11, Unit 2 Time 11, Unit 7 Time 16, Unit 4 Time 18, Unit 6 Time 19\n")

cat("\nPlots created for each sample dataset:\n")
cat("  - Combined time series: All units on same plot with break markers\n")
cat("  - Concatenated time series: All units placed end-to-end in one long line\n")

cat("\nFiles created:\n")
cat(" - C++ headers: sim_*.h\n")
cat(" - Master C++ header: simulation_datasets.h\n")
cat(" - Example usage: example_usage.cpp\n")

cat("\nTo use in C++:\n")
cat("  #include \"simulation_datasets.h\"\n")
cat("  auto datasets_050 = get_datasets_by_step_mean(0.50);\n")
cat("  auto datasets_100 = get_datasets_by_step_mean(1.00);\n")
cat("  // etc.\n")

cat("\n", paste(rep("=", 60), collapse=""), "\n", sep = "")