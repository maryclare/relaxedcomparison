# library(bestsubset)
library(glmnet)
glmnet.control(fdev=0)
# library(devtools)
# Install the full XCode application from App Store. Command line tools are not enough.
# Do not uninstall XCode!!! It will mess up your R package installation.
# install_github("maryclare/powopt")
library(powopt)
rm(list=ls())

enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}

sim.xy = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1) {
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  
  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }
  else Sigma = diag(1,p)
  
  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  enlist(x,y,xval,yval,Sigma,beta,sigma)
}

n = 70; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set

# nrep should be 500, but we use 10 for speed
nrep = 500 # Number of repetitions

# nrep = 500 # Number of repetitions
seed = 1 # Random number generator seed
set.seed(seed)
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type

xy.obj = sim.xy(n,p,nval,rho=0,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 50 # nlam = 31

r <- 1
cat(r,"... ")
eps = rnorm(n)*sigma
y = mu + eps

x_standardized <- apply(x, 2, function(z) {(z)/(sqrt(mean(z^2) - mean(z)^2))})

# Calculate <x_j, y> for each j on standardized x
x_y_dot_standardized <- colSums(x_standardized * y)

# Find the maximum absolute value of <x_j, y> on standardized x
max_abs_x_y_dot_standardized <- max(abs(x_y_dot_standardized))

q = 0.25

# Calculate lambda_max
lambda_max_standardized_25 <- (max_abs_x_y_dot_standardized / (2-q) / n) ^(2-q) * (2-2*q)^(1-q)


lamseq <- c(lambda_max_standardized_25, rep(NA, nlam - 1))
nz <- c(30, rep(NA, nlam - 1))

beta.pow_25 = matrix(0,p, nlam)
beta.pow_25[,1] = powCD(X = x_standardized, y = y, 
                        sigma.sq = n, lambda = lamseq[1], 
                        q = q, rand.restart = 0, start = rep(0,p))


epsilon <- 0.001

for (i in 2:nlam) {  # Iterate over lambda sequence starting from the second element
  cat("i=", i, "\n")  # Print current index for debugging
  
  # Propose a new lambda slightly smaller than the previous one (log-scale decrement)
  lamnew <- exp(log(lamseq[i - 1]) - epsilon)
  # lamnew <- exp(log(lamseq[i - 1]) - epsilon + runif(1, -delta_rand, delta_rand))
  
  # Compute beta estimate at new lambda using powCD, starting from previous solution
  betacurr <- powCD(X = x_standardized, y = y, 
                    sigma.sq = n, lambda = lamnew, q = q, 
                    rand.restart = 0, start = beta.pow_25[,i-1])
  
  # Keep decreasing lambda until the sparsity pattern changes
  while (sum(betacurr == 0) == sum(beta.pow_25[,i-1] == 0)) {
    lamnew <- exp(log(lamnew) - epsilon)  # Further decrease lambda
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = n, lambda = lamnew, q = q, 
                     rand.restart = 0, start = betacurr)  # Recompute beta
  }
  
  # Check how the number of zero coefficients has changed
  if (sum(betacurr == 0) != (sum(beta.pow_25[,i-1] == 0) - 1)) {
    # If not exactly one fewer zero, we need to search the correct lambda interval
    lower <- lamnew
    upper <- lamseq[i - 1]
  } else {
    # If exactly one fewer zero, then lamnew is the desired lambda
    lower <- upper <- lamnew
  }
  
  # Binary search to find a lambda that gives exactly one fewer zero coefficient
  while (sum(betacurr == 0) != (sum(beta.pow_25[,i-1] == 0) - 1) & round(lower, 8) != round(upper, 8)) {
    mid <- (lower + upper)/2  # Midpoint lambda
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = n, lambda = mid, q = q, 
                     rand.restart = 0, start = beta.pow_25[,i-1])  # Recompute beta
    
    if ((sum(betacurr == 0) < (sum(beta.pow_25[,i-1] == 0) - 1))) {
      # If too few zeros, lambda is too small → increase lower bound
      lower <- mid
    } else {
      # Otherwise decrease upper bound
      upper <- mid
    }
  }
  
  lamseq[i] <- lamnew  # Save final lambda
  beta.pow_25[,i] <- betacurr  # Save corresponding beta
  nz[i] <- sum(betacurr == 0)  # Record number of zero coefficients
  cat("nz=", nz[i], "\n")  # Print for debug
  
  if (nz[i] == 0) {
    break  # Stop if no more zero coefficients (fully dense)
  }
}

par(mfrow = c(1, 2))
plot(log(lamseq), nz)

abline(h = 0:30, lty = 3)
plot(log(lamseq))
nz
unique(nz)
length(na.omit(unique(nz)))

plot(nz,log(lamseq))


##########################################################################################
# Set the target number of zeros.
# For example, if p = 30 and we want 28 zeros, that means only 2 nonzeros.
target_nz <- 28         # Target number of zeros
max_iter <- 500         # Maximum iterations to prevent infinite loops
epsilon <- 0.001        # Step size parameter (same as before)
found <- FALSE          # Flag to indicate if we have found a solution with target nz

# Start from the initial solution at lambda_max (which corresponds to the sparsest solution)
curr_lambda <- lamseq[1]
curr_beta <- beta.pow_25[,1]
curr_nz <- sum(curr_beta == 0)
iter <- 1

cat("Starting search for target nz =", target_nz, "\n")
while (iter <= max_iter && curr_nz != target_nz) {
  # Update lambda using a log-scale decrement
  new_lambda <- exp(log(curr_lambda) - epsilon)
  # Compute the candidate beta using the updated lambda, starting from the previous beta
  new_beta <- powCD(X = x_standardized, y = y, 
                    sigma.sq = n, lambda = new_lambda, q = q, 
                    rand.restart = 0, start = curr_beta)
  new_nz <- sum(new_beta == 0)
  
  cat("Iter", iter, ": lambda =", new_lambda, "nz =", new_nz, "\n")
  
  # Check if the target nz is between the current and the new solution's nz values
  if ((curr_nz > target_nz && new_nz < target_nz) || 
      (curr_nz < target_nz && new_nz > target_nz)) {
    cat("Target nz is between", curr_nz, "and", new_nz, ", starting binary search.\n")
    # Set the lower and upper bounds for lambda in the binary search
    lower <- new_lambda
    upper <- curr_lambda
    # Binary search for lambda that yields exactly the target nz
    while (round(lower, 12) != round(upper, 12)) {
      # Compute the midpoint and add a small random perturbation (optional)
      mid <- (lower + upper) / 2 + runif(1, -epsilon*0.1, epsilon*0.1)
      mid_beta <- powCD(X = x_standardized, y = y, 
                        sigma.sq = n, lambda = mid, q = q, 
                        rand.restart = 0, start = curr_beta)
      mid_nz <- sum(mid_beta == 0)
      
      cat("  Binary search: lambda =", mid, "nz =", mid_nz, "\n")
      
      # If the midpoint gives the target nz, save and break out of the loop
      if (mid_nz == target_nz) {
        new_beta <- mid_beta
        new_lambda <- mid
        found <- TRUE
        break
      } else if ((mid_nz < target_nz && curr_nz > target_nz) || 
                 (mid_nz > target_nz && curr_nz < target_nz)) {
        # If the midpoint nz is on the same side as the new solution, update lower bound
        lower <- mid
      } else {
        # Otherwise, update the upper bound
        upper <- mid
      }
    }
    if (found || sum(new_beta == 0) == target_nz) {
      cat("Found solution with target nz =", target_nz, "at lambda =", new_lambda, "\n")
      curr_beta <- new_beta
      curr_lambda <- new_lambda
      curr_nz <- sum(curr_beta == 0)
      break
    }
  } else {
    # If target nz is not bracketed, continue decreasing lambda using the original step
    curr_lambda <- new_lambda
    curr_beta <- new_beta
    curr_nz <- new_nz
  }
  
  iter <- iter + 1
}

if (curr_nz == target_nz) {
  cat("Final solution: lambda =", curr_lambda, "with nz =", curr_nz, "\n")
} else {
  cat("Did not converge to target nz =", target_nz, "after", max_iter, "iterations.\n")
}
