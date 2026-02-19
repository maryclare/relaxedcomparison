# load libraries

library(bestsubset)
glmnet.control(fdev=0)
library(devtools)
# install_github("maryclare/powopt")
library(powopt)
# install_github("duanning1/powwrap")
# library(powwrap)
library(L0Learn)
rm(list=ls())
dir <- "/home/nduan/Research/"

## Define three functions powCD_wrapper, coef.powCD, predict.powCD

# We need number of tuning values of q, ''nq''. The default of ''nq'' should be 1!
# We need number of tuning values of q, ''nq''. The default of ''nq'' should be 1!
powCD_wrapper <- function(x, y, nlambda, nq=1) {
  # Standardize data (based on glmnet code)
  yscale <-sqrt(sum(y^2)) 
  ystd <- y/yscale
  xscales <- apply(x, 2, function(xx) {
    sqrt(mean((xx - mean(xx))^2))
  })
  xstd <- sweep(x, 2, xscales, "/")
  sv <- svd(xstd)
  bsv <- c(sv$v%*%diag(ifelse(sv$d > 0, 1/sv$d, 0))%*%t(sv$u)%*%y)
  bm <- crossprod(xstd, ystd)/diag(crossprod(xstd))
  # q <- 0.5
  
  # Set up the q sequence if nq > 1, otherwise use the single value of q
  # q_seq <- if (nq > 1) seq(0.1, 1, length.out = nq) else q
  q_seq <- seq(0.1, 1, length.out = nq)
  
  # Initialize the coefficient matrix to hold all values of beta.pow for different q
  p <- ncol(x)
  total_cols <- nlambda * length(q_seq)
  beta.pow_all <- matrix(0, p, total_cols)
  
  # Initialize the column names for the final beta.pow_all matrix
  col_names <- character(total_cols)
  
  # Initialize column counter for beta.pow_all
  col_counter <- 1
  
  # Loop over each q in the q sequence
  for (q_val in q_seq) {
    # Calculate lambda_max for the current q
    n <- nrow(x)
    
    lammax <- max((2*(1 - q_val))^(1 - q_val)*(2-q_val)^(q_val-2)*abs(crossprod(xstd, ystd))^(2-q_val)*diag(crossprod(xstd))^(q_val - 1))
    eps <- 10^(-4)
    # eps = lam*sum(|b|^q_val)
    # eps/sum(|b|^q_val) = lam
    
    # lammin <- eps/sum(abs(bsv)^q_val)
    lammin <- eps/sum(abs(bm)^q_val)
    
    mylams <- exp(seq(log(lammax), 
                      log(lammin),
                      length.out = nlambda))
    
    betas <- matrix(nrow = ncol(xstd), ncol = length(mylams))
    start <- rep(0, ncol(xstd))
    for (k in 1:ncol(betas)) {
      myest <- powCD(y = ystd, X = xstd, q = q_val, sigma.sq = n, 
                     lambda = mylams[k]/n, tol = 10^(-7), # setting tol = 10^(-7) for convergence
                     start = start)
      start <- myest
      betas[, k] <- myest
    }
    obetas <- sweep(betas, 1, xscales, "/")*yscale
    
    # Add the results to the main beta.pow_all matrix and name the columns
    beta.pow_all[, col_counter:(col_counter + nlambda - 1)] <- obetas
    
    # Create column names for this range of lambda values
    for (i in 1:nlambda) {
      col_names[col_counter] <- paste0("lambda = ", round(mylams[i], 4), 
                                       ", q = ", round(q_val, 2))
      col_counter <- col_counter + 1
    }
  }
  
  # Set the column names of the beta.pow_all matrix
  colnames(beta.pow_all) <- col_names
  
  # Create the output object (same as before, but now with all betas for different q's)
  obj <- list("coeff" = beta.pow_all, "lambda_seq" = mylams)
  
  # Assign class
  class(obj) <- "powCD"
  
  # Return the object
  return(obj)
}


coef.powCD <- function(out) {
  # Extract the 'coeff' matrix from the 'out' object
  coeff_matrix <- out[["coeff"]]
  
  if (is.null(out$coeff)) stop("coef.powCD: Coefficients are NULL!")
  # Return the coefficient matrix
  return(coeff_matrix)
}

predict.powCD <- function(out, xnew) {
  # Perform matrix multiplication between new data and the coefficient matrix
  prediction <- xnew %*% out[["coeff"]]
  
  if (is.null(xnew)) stop("predict.powCD: xnew is NULL!")
  if (is.null(out$coeff)) stop("predict.powCD: Coefficients are NULL!")
  
  # Return the predictions
  return(prediction)
}

# Adopted from sim.hi5.R

## Low-dimensional simulation
# library(bestsubset)
# library(L0Learn)
# glmnet.control(fdev=0)

# Set some overall simulation parameters
n = 50; p = 50 # Size of training set, and number of predictors
nval = n # Size of validation set

#We want to increase nrep from 10 to 100 for low settings
nrep = 100 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
type.vec = c(1:3,5) # Simulation settings to consider
rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios 
stem = paste0("sim.n",n,".p",p)

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)
reg.funs[["Forward stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
#reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE,
#                                             time.limit=1800,
#                                             params=list(Threads=4))
reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,
                                                  nrelax=10,nlam=50)
reg.funs[["PowCD"]] = function(x,y) powCD_wrapper(x,y,nlambda=50, nq = 10)

# Also incorporate L0Learn algorithms from Hazimeh and Mazumder (2017)
# reg.funs[["L0Learn 1"]] = function(x,y) L0Learn.fit(x,y,penalty="L0",
#                                                     algorithm="CDPSI",
#                                                     nLambda=50)
# reg.funs[["L0Learn 2"]] = function(x,y) L0Learn.fit(x,y,penalty="L0L1",
#                                                     algorithm="CDPSI",
#                                                     nGamma=10,nLambda=50)

## NOTE: the loop below was not run in serial, it was in fact was split up
## and run on a Linux cluster

# Get job indicator i.e. which row of the above matrix we're in
slurm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Total combinations
n_beta = length(type.vec)
n_rho = length(rho.vec)
n_snr = length(snr.vec)
total_combinations = n_beta * n_rho * n_snr

# Adjust SLURM task ID to start from 0
slurm_id <- slurm_id - 1

# Calculate the indices for each loop based on the slurm_id
beta_index <- slurm_id %/% (n_rho * n_snr) + 1
rho_index <- (slurm_id %% (n_rho * n_snr)) %/% n_snr + 1
snr_index <- (slurm_id %% n_snr) + 1

# Extract the corresponding values for this job
beta.type <- type.vec[beta_index]
rho <- rho.vec[rho_index]
snr <- snr.vec[snr_index]

# Create file name
name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
file = paste0("Rep_rds_lo_2025_new/", name, ".snr", round(snr, 2), ".rds")

# Output the details
cat("..... NEW SIMULATION .....\n")
cat("--------------------------\n")
cat(paste0("File: ", file, "\n\n"))

# Run the simulation for this specific combination
sim.master(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=seed, s=s,
           verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)

# Add file to the list
# file.list = c(file.list, file)

cat("\n")
