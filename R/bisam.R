#' Bayesian Indicator Saturation Model
#'
#' @param data A matrix with columns for unit indices, time indices, y variable, and X variables
#' @param i_index Column index for unit identifier (1-based R indexing)
#' @param t_index Column index for time identifier (1-based R indexing)
#' @param y_index Column index for dependent variable (1-based R indexing)
#' @param Ndraw Number of MCMC draws
#' @param Nburn Number of burn-in iterations
#' @param b_prior Prior type for coefficients ("g" or "f" or "hs")
#' @param lambda_b Parameter for g/f prior
#' @param c0 Parameter for inverse gamma prior on variance
#' @param C0 Parameter for inverse gamma prior on variance
#' @param va Parameter
#' @param vb Parameter
#' @param tau Parameter for model selection
#' @param priorDelta Prior object for model selection (modelbbprior or modelbinomprior)
#' @param geweke Boolean for Geweke test
#' @param use_phiinit Boolean for using initial phi values. TRUE -> Use our sampling, FALSE use rnlp sampling
#' @param const_val Boolean for including a constant
#' @param ife Boolean for individual fixed effects
#' @param tfe Boolean for time fixed effects
#' @param iis Boolean for indicator saturation
#' @param sis Boolean for stepshift saturation
#' @param computation_strategy Computation strategy (STANDARD, SPLIT_SEQUENTIAL, SPLIT_PARALLEL)
#' @param new_par_method New parameter method
#' @param new_par_hesstype New parameter hessian type
#' @param new_par_optim_method New parameter optimization method
#' @param new_par_optim_maxit New parameter optimization max iterations
#' @param new_par_B New parameter B
#' @param new_par_knownphi New parameter known phi
#' @param new_par_r New parameter r
#' @param new_par_alpha New parameter alpha
#' @param new_par_lambda New parameter lambda
#'
#' @return A list containing MCMC results and summary statistics
#' @export
estimate_model <- function(
    data,
    i_index = 1,
    t_index = 2,
    y_index = 3,
    Ndraw = 5000,
    Nburn = 500,
    b_prior = "g",
    lambda_b = 100.0,
    c0 = 0.001,
    C0 = 0.001,
    va = 1.0,
    vb = 1.0,
    tau = 1.0,
    priorDelta = modelbbprior(1, 1),
    geweke = FALSE,
    use_phiinit = TRUE,
    const_val = FALSE,
    ife = FALSE,
    tfe = FALSE,
    iis = TRUE,
    sis = TRUE,
#     new_par_include_vars = c(),
    new_par_method = 1,
    new_par_hesstype = 1,
    new_par_optim_method = 1,
    new_par_optim_maxit = 100,
    new_par_B = 10,
    new_par_knownphi = 0,
    new_par_r = 1,
    new_par_alpha = 0.05,
    new_par_lambda = 1.0,
    computation_strategy = SPLIT_SEQUENTIAL()  # Default to SPLIT_SEQUENTIAL
) {
    # Convert data to matrix if it's a data frame
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    # Adjust indices for R to C++ (0-based indexing)
    i_index_cpp <- i_index - 1
    t_index_cpp <- t_index - 1
    y_index_cpp <- y_index - 1

    a <- NA_real_
    b <- NA_real_
    p <- NA_real_

    # Extract the prior variables depending on the type
    if (inherits(priorDelta, "modelbbprior")) {
          a = priorDelta$a,
          b = priorDelta$b,
          p = NA_real_
      } else if (inherits(priorDelta, "modelbinomprior")) {
          a = NA_real_,
          b = NA_real_,
          p = priorDelta$p
      } else {
        stop("priorDelta must be either a modelbbprior or modelbinomprior object")
      }




    # Call the C++ function
    result <- rcpp_estimate_model(
        data, i_index_cpp, t_index_cpp, y_index_cpp, Ndraw, Nburn, b_prior,
        lambda_b, c0, C0, va, vb, tau, use_phiinit,
        const_val, ife, tfe, iis, sis,
#          new_par_include_vars,
        new_par_method, new_par_hesstype, new_par_optim_method, new_par_optim_maxit,
        new_par_B, new_par_knownphi, new_par_r, new_par_alpha, new_par_lambda,
        computation_strategy
    )

    return(result)
}

#' Constants for computation strategy
#'
#' @export
STANDARD <- function() {
  comp_strategy_standard()
}

#' @export
SPLIT_SEQUENTIAL <- function() {
  comp_strategy_split_sequential()
}

#' @export
SPLIT_PARALLEL <- function() {
  comp_strategy_split_parallel()
}



#' Create a Beta-Binomial prior object
#'
#' @param a First parameter of the Beta-Binomial prior (shape1)
#' @param b Second parameter of the Beta-Binomial prior (shape2)
#' @return A modelbbprior object
#' @export
modelbbprior <- function(a = 1, b = 1) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1 || length(b) != 1) {
    stop("Parameters 'a' and 'b' must be single numeric values")
  }
  if (a <= 0 || b <= 0) {
    stop("Parameters 'a' and 'b' must be positive")
  }

  structure(
    list(a = as.double(a), b = as.double(b)),
    class = "modelbbprior"
  )
}

#' Create a Binomial prior object
#'
#' @param p Probability parameter for the Binomial prior
#' @return A modelbinomprior object
#' @export
modelbinomprior <- function(p = 0.5) {
  if (!is.numeric(p) || length(p) != 1) {
    stop("Parameter 'p' must be a single numeric value")
  }
  if (p < 0 || p > 1) {
    stop("Parameter 'p' must be between 0 and 1")
  }

  structure(
    list(p = as.double(p)),
    class = "modelbinomprior"
  )
}