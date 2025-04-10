# detach(package:BISAM, unload = TRUE)
# remove.packages("BISAM")
# q(save="yes")
# library(devtools)
# install_github("Avriox/BISAM", force = T)
rm(list = ls())

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("sparseMatrixStats", force = TRUE)

#remove.packages("mombf")
#devtools::install("../mombf/mombf")


# setwd("~/uni/wu/sis_project")

source("./contr_sim_breaks_fun.R")

# source("./non_local_ism_fun.R")

set.seed(192837612)
n <- 3 # number of sim. observations
t <- 20 # number of sim. time periods
nx <- 3 # number of regressors
const <- FALSE # inclusion of a constant
ife <- F # inclusion of indiv. fixed effects
tfe <- F # inclusion of time fixed effects
iis <- T # inclusion of indicator saturation
sis <- T # inclusion of stepshift saturation
p.outl <- 0.0 # probability of outlier in a Series
p.step <- 0.0 # probability of a stepshift in a Series
outl.mean <- 0 # mean of size of outlier
outl.sd <- 0 # variance of size of outlier
step.mean <- 10 # mean of size of stepshift
step.sd <- 0.00 # variance of size of stepshift
error.sd <- 1 # variance of the error

pos.outl <- 0
pos.step <- c(10)

# Function to simulate new data or load from file
simulate_or_load_data <- function(filepath) {
  # if (file.exists(filepath)) {
  #   load(filepath)
  #   cat("Data loaded from file.\n")
  # } else {
  
  sim <- contr_sim_breaks(
    n = n, t = t, nx = nx, iis = iis, sis = sis,
    const = const, ife = ife, tfe = tfe,
    pos.outl = pos.outl, pos.step = pos.step,
    outl.mean = outl.mean, step.mean = step.mean,
    error.sd = error.sd
  )
  
  # save(sim, file=filepath)
  cat("Data simulated and saved to file.\n")
  # }
  return(sim)
}

# Define the file path where the data will be stored
data_filepath <- "sim_data2.RData"

# Simulate or load data
sim <- simulate_or_load_data(data_filepath)

data <- sim$data
const <- FALSE
tfe <- FALSE
ife <- FALSE
sis <- sis
i_index <- 1
t_index <- 2
y_index <- 3
Ndraw <- 5000L
Nburn <- 500L
b_prior <- "g"
lambda_b <- 100
c0 <- 0.001
C0 <- 0.001
va <- 1
vb <- 1
tau <- 1

geweke <- FALSE

if(geweke){
  library(extraDistr)
  lambda_b = lambda_g = 1.234
  Ndraw= 100000
  c0=3
  C0=1
  sis=TRUE
  b_prior="g"
}
#Setup
require(mvtnorm)
library(Matrix)
require(dplyr)
library(mombf)
# library(BISAM)
library(glmnet)
# require(Rcpp)
# require(RcppArmadillo)
# sourceCpp("./00_code/03_ssvs_project/01_code/working/sis_proj_JKB.cpp")


start_time <- Sys.time()

colnames(data)[-c(y_index, i_index, t_index)] <- 
  paste0("beta.", colnames(data[, -c(y_index, i_index, t_index)]))

y  <- as.matrix(data[,y_index])
X  <- X_ <- as.matrix(data[,-c(y_index,i_index,t_index)])

n <- length(unique(data[,i_index]))
t <- length(unique(data[,t_index]))
N <- n*t

#Build X matrix
if(const) {X <- cbind(X, "const" = 1)}
if(ife) { # Individual fixed effects
  IFE <- kronecker(diag(n), matrix(1, t))
  colnames(IFE) <- paste0("ife.", unique(data[, i_index]))
  if(const) {IFE <- IFE[,-1]}
  X <- cbind(X, IFE)
}
if(tfe) { # Time fixed effects
  TFE <- kronecker(matrix(1, n), diag(t))
  colnames(TFE) <- paste0("tfe.", unique(data[, t_index]))
  if(const) {TFE <- TFE[,-1]}
  X <- cbind(X, TFE)
}
if(tfe & ife & !const) {
  warning("Both time and unit fixed effects used.\nDropping first indiv. FE to avoid perfect colinearity")
  X <- X[, -(ncol(X_)+1)]
}

#Find matrix dimensions
p <- ncol(X)

z <- lower.tri(matrix(1,nrow=t,ncol=t),diag = T)[,-c(1:2,t)]
mode(z) <- "integer"
Z <- kronecker(diag(n),z)
mode(Z) <- "integer"
sis_grid <- expand.grid(unique(data[, t_index])[-c(1:2,t)], unique(data[, i_index]))
colnames(Z) <- paste0("sis.", sis_grid[, "Var2"], ".", sis_grid[, "Var1"])

r <- ncol(Z)

XX <- crossprod(X)
ZZ <- crossprod(Z)

#Prior Parameters:

if(b_prior=="g") b0 <- rep(0,p)
if(b_prior=="f") b0 <- solve(XX)%*%t(X)%*%y
# if(b_prior=="flasso"){
#   lasso_model <- glmnet(X, y, alpha = 1, intercept = FALSE)
#   best_lambda <- cv.glmnet(X, y, alpha = 1)$lambda.min
#   b0 <- (coef(lasso_model, s = best_lambda) %>% as.array())[-1]
#   b0 <- b0
# }
B0     <- solve(XX)*lambda_b
B0_inv <- XX/lambda_b

#Starting Values
b_i      <- rep(1,p)
g_i      <- rep(1,r)
g_incl_i <- g_i # always non zero for rnlp within Gibbs; last value conditional on inclusion
s2_i     <- 1

w_i    <- rep(0,r)
z_cols <- rep(1:ncol(z),n)
w_1    <- z_cols * w_i # columns where w_i is 1

#Starting Values for HS-Plug-In
lamb_b <- rep(1,p)
tau_b  <- 1
nu_b   <- rep(1,p)
xi_b   <- 1


#Prep Calculations
cN     <- c0+N/2
XX_inv <- solve(XX)
Pr_g   <- rep(NA,n*(t-2))

#Store
if(Nburn>=Ndraw){
  stop('The number of burn-in exceeds number of draws')
} else {
  Nstore  <- Ndraw-Nburn
}
#main
b_store <- matrix(NA,nrow=Nstore,ncol=p)
g_store <- matrix(NA,nrow=Nstore,ncol=r)
s2_store<- matrix(NA,nrow=Nstore,ncol=1)
colnames(b_store)<-colnames(X)
colnames(g_store)<-colnames(Z)
colnames(s2_store)<-'sigma2'
#selection
w_store <- matrix(NA,nrow=Nstore,ncol=r)
colnames(w_store)<-colnames(Z)
#Start Gibbs Sampler
pb <- txtProgressBar(min = 0, max = Ndraw, style = 3)

#$$$$$$$$$$$$$$$$$$$$$ START LOOP $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for(i in (1-Nburn):Nstore){
  #i = (1-Nburn)
  
  #================ Geweke Test =================.====
  if(geweke){
    e <- rnorm(N,0,s2_i^0.5)
    y <- X%*%b_i + Z%*%g_i + e
  }
  
  #=================== draw p(s2|a,b,g,y) ==========.====
  # cN
  CN  <- C0+1/2*crossprod(y-X%*%b_i-Z%*%g_i)
  s2_i<- 1/rgamma(1,shape=cN, rate=CN)
  
  s2_i <- sqrt(2)
  
  #=================== draw p(b|a,g,s2,y) ==========.====
  if(b_prior == "hs"){
    # draw p(xi,nu|b,s2,y) ========.====
    nu_b <- sapply(lamb_b,function(ll) 1/rgamma(1, shape = 1, rate = 1+1/ll))
    nu_b[nu_b > 1e+8] <- 1e+8
    nu_b[nu_b < 1e-8] <- 1e-8
    xi_b <- 1/rgamma(1, shape = 1, rate = 1+1/tau_b)
    if(xi_b>1e+8) xi_b<-1e+8
    if(xi_b<1e-8) xi_b<-1e-8
    # draw p(tau,lamb|xi,nu,b,s2,y).====
    lamb_b<- mapply(function(vvi,bbi) 1/rgamma(1, shape = 1, rate = 1/vvi + bbi^2/(2*tau_b*s2_i)), nu_b, b_i)
    lamb_b[lamb_b > 1e+8] <- 1e+8
    lamb_b[lamb_b < 1e-8] <- 1e-8
    tau_b <- 1/rgamma(1, shape = (p+1)/2, rate = 1/xi_b + sum(b_i^2/lamb_b)/(2*s2_i))
    if(tau_b>1e+8) tau_b<-1e+8
    if(tau_b<1e-8) tau_b<-1e-8
    
    A <- XX + 1/tau_b * diag(1/(lamb_b))
    A_inv <- solve(A)
    BN<- s2_i * A_inv
    bN<- A_inv %*% crossprod(X,y-Z%*%g_i)
  } else if(b_prior == "g" || b_prior == "f"){
    BN_inv <- (1/s2_i+1/lambda_b)*XX
    BN     <- s2_i*lambda_b/(s2_i+lambda_b)*XX_inv
    bN     <- 1/(s2_i+lambda_b)*(lambda_b*XX_inv%*%crossprod(X,y-Z%*%g_i)+s2_i*b0)
  }
  
  b_i <- try((bN + t(chol(BN)) %*% rnorm(p)), silent = TRUE)
  if(is(b_i, "try-error")) {
    b_i <- t(rmvnorm(1,bN,BN))
  }
  
  b_i <- c(1,-1,2)
  #================ draw p(w|a,b,s2,y) =========.====
  y_hat <- y-X%*%b_i
  
  nn <- 2
  if (i==1-Nburn){
    initpar = "auto"
  }else {
    initpar = g_i
  }
  
  
  # y = y_hat;
  # x = Z;
  # groups=1:ncol(x);
  # nknots=9;
  # center = FALSE;
  # scale = FALSE;
  # enumerate = FALSE;
  # includevars=rep(FALSE,ncol(x));
  # niter = nn;
  # thinning=1;
  # burnin = nn-1;
  # family = "normal";
  # priorCoef = imomprior(tau = tau);
  # priorDelta = modelbinomprior(0.50);
  # phi = s2_i;
  # deltaini = as.logical(w_i);
  # initSearch = 'none';
  # method = 'ALA';
  # hess = "asymp";
  # initpar = 'auto';#g_i; # fixed in initParameters.R function getthinit
  # #begin for debugging
  # adj.overdisp='intercept';
  # optimMethod="auto";
  # optim_maxit=10;
  # B=10^5;
  # priorVar=igprior(.01,.01);
  # priorSkew=momprior(tau=0.348);
  # #end for debugging
  # XtXprecomp = TRUE;
  # verbose = FALSE;
  
  
  w_i_mod <- modelSelection( # in modelSelection.R -> General model selection routines
    y = y_hat,
    x = Z,
    center = FALSE,
    scale = FALSE,
    enumerate = FALSE,
    niter = nn,
    burnin = nn-1,
    family = "normal",
    priorCoef = imomprior(tau = tau),
    # priorDelta = modelbinomprior(0.50),
    priorDelta = modelbbprior(1,1),
    phi = s2_i,
    deltaini = as.logical(w_i),
    initSearch = 'none',
    method = 'ALA',
    hess = "asymp",
    initpar = initpar, # fixed in initParameters.R function getthinit
    #begin for debugging
    adj.overdisp='intercept',
    optimMethod="auto",
    B=10^5,
    priorVar=igprior(.01,.01),
    #end for debugging
    XtXprecomp = TRUE,
    verbose = FALSE #,
    #   split_x = FALSE, # added by Jakob 
    #   n_observations = n, # added by Jakob 
    #   n_timeperiods = t # added by Jakob 
  )
  
  w_i <- w_i_mod$postSample
  # w_i_mod$margpp
  # w_i_mod$
  # sum(w_i)
  # w_i_mod$postSample |> is.na() |> any()
  # w_i_mod <- mombf::postModeBlockDiag()
  
  # when i-mom approx ready, consider block search:
  # w_i <- mombf::modelsearchBlockDiag(
  #   y = y_hat,
  #   x = Z,
  #   priorCoef = BISAM::imomprior(tau = tau),
  #   priorDelta = BISAM::modelbinomprior(0.50),
  #   blocksize = t-3,
  #   verbose = TRUE
  # )
  
  if(geweke){
    if(i == (1-Nburn)){
      cat("\nCAREFUL: No Variable Selection in geweke test!\n")
      pb <- txtProgressBar(min = 0, max = Ndraw, style = 3)
    }
    w_i <- rep(1,r)
  }
  
  w_1_Z <- 1:r * w_i
  #================ draw p(g|w,a,b,s2,y) =========.====
  # g_i <- rnlp(msfit=w_i_mod, niter=3, burnin=2)[-c(1,r+2)]
  # g_i <- coef(w_i_mod,niter = 3,burnin = 2, meanonly = TRUE)[-c(1,r+2)] # in modelSelection.R, fixed to only exit means
  if(sum(w_i)>0){
    g_i    <- rep(0,r)
    # set.seed(22)
    g_draw <- rnlp(
      niter = 5,
      burnin = 4,
      thinning = 1,
      y = y_hat,
      x = Z[,w_1_Z,drop=F],
      priorCoef = imomprior(tau=tau), 
      priorVar = igprior(100000, 100000), 
      outcometype = "Continuous", 
      family = "normal"
      #,
      # tau = tau
      # ,
      # a_phi = 1, # variance parameter
      # b_phi = 1, # variance parameter
      # prior = 1#, # imom
    )
      #thini = g_incl_i,
      #phiini = s2_i)
    g_draw
    
    # g_draw <- BISAM::rnlp(
    #   niter = 10,
    #   burnin = 9,
    #   thinning = 1,
    #   y = y_hat,
    #   x = Z[,w_1_Z,drop=F],
    #   priorCoef = BISAM::imomprior(tau = tau),
    #   priorVar = BISAM::igprior(1000000,1000000),
    #   outcometype = "Continuous",
    #   family = "normal"
    #   )[,-r-1]
    
    # if(any(is.nan(g_draw))){
    #   g_draw <- BISAM::rnlp_wrapper(
    #     niter = 2,
    #     burnin = 1,
    #     thinning = 1,
    #     y = y_hat,
    #     x = Z[,w_1_Z,drop=F],
    #     tau = 3*1e+4,
    #     a_phi = 1,
    #     b_phi = 1,
    #     prior = 1,
    #     thini = g_incl_i,
    #     phiini = s2_i)
    # }
    g_i[w_1_Z]      <- g_draw
    g_incl_i[w_1_Z] <- g_draw
  } else{
    g_i <- rep(0,r)
  }
  if(any(is.nan(g_i))){print("NaNs Produced")}
  # g_i <- rep(0,r)
  #=================== store and adjust tau ===========.====
  
  if(i<=0){ # adjust tau to fit FP rate
    # if(sum(g_i) < n*t*0.05){
    #   tau <- tau * 5-sum(g_i)
    # }
  } else{ # store
    b_store[i,] <-b_i
    g_store[i,] <-g_i
    w_store[i,] <-w_i
    s2_store[i,]<-s2_i
  }
  setTxtProgressBar(pb, (i + Nburn))
}
timer <- Sys.time() - start_time
close(pb)
cat("Finished after ", format(round(timer, 2)), ".\n", sep = "")