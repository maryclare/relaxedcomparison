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
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type

set.seed(seed)
xy.obj = sim.xy(n,p,nval,rho=0,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 300

r <- 1
cat(r,"... ")
eps = rnorm(n)*sigma
y = mu + eps



############################################################################################################
X <- matrix(rnorm(n*p, 0, 0.2), n, p)
X <- scale(X)
u <- svd(X)$u
X <- u
Y <- rnorm(n, 0, 1)
############################################################################################################

x_standardized <- apply(x, 2, function(z) {(z)/(sqrt(mean(z^2) - mean(z)^2))})

# Calculate <x_j, y> for each j on standardized x
x_y_dot_standardized <- colSums(x_standardized * y)

# Find the maximum absolute value of <x_j, y> on standardized x
max_abs_x_y_dot_standardized <- max(abs(x_y_dot_standardized))

q = 0.25

# Calculate lambda_max
lambda_max_standardized_25 <- (max_abs_x_y_dot_standardized / (2-q) / n) ^(2-q) * (2-2*q)^(1-q)


# Extract lambda_max from the glmnet model
# lambda_max_glmnet <- fit.lasso$lambda[1]

epsilon <- 0.00001

# Calculate lambda_min
lambda_min_25 <- epsilon * lambda_max_standardized_25

# Generate the sequence of lambda values on a log scale
#a <- 0.5
#lambda_seq_manual_25 <- seq((lambda_max_standardized_25)^(a), 
#                            (lambda_min_25)^(a), length.out = nlam)^(1/a)

# lambda_seq_manual_25 <- exp(seq(log(lambda_max_standardized_25), 
#                                log(lambda_min_25), length.out = nlam))

sorted_abs_x_y_dot_standardized <- sort((abs(x_y_dot_standardized) / (2-q) / n)^(2-q) * (2 - 2*q)^(1-q), 
                                        decreasing = TRUE)
refine_lambda_seq <- function(seq, num_points = 30) {
  refined_seq <- c()
  for (i in seq_along(seq)[-length(seq)]) {
    refined_seq <- c(refined_seq, seq(from = seq[i], to = seq[i + 1], length.out = num_points + 2)[-1])
  }
  return(refined_seq)
}
num_lambda <- 30  # Number of points between each pair
lambda_seq_manual_25 <- refine_lambda_seq(sorted_abs_x_y_dot_standardized, num_lambda)

beta.pow_25 = matrix(0,p, length(lambda_seq_manual_25))
beta.pow_25[,1] = powCD(X = x_standardized, y = y, 
                        sigma.sq = n, lambda = lambda_seq_manual_25[1], 
                        q = q, rand.restart = 0, start = rep(0,p))
for (i in 2:length(lambda_seq_manual_25)) {
  beta.pow_25[,i] = powCD(X = x_standardized, y = y, 
                          sigma.sq = n, lambda = lambda_seq_manual_25[i], q = q, 
                          max.iter = 10000, tol = 1e-7,
                          rand.restart = 0, start = beta.pow_25[,i-1])
}

nzs.pow_25 = colSums(beta.pow_25 != 0)
plot(log(lambda_seq_manual_25), nzs.pow_25)
unique(nzs.pow_25)
length(unique(nzs.pow_25))
