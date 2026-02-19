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
nlam = 31

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

for (i in 2:nlam) {
  cat("i=", i, "\n")
  lamnew <- exp(log(lamseq[i - 1]) - epsilon)
  betacurr <- powCD(X = x_standardized, y = y, 
                    sigma.sq = n, lambda = lamnew, q = q, 
                    rand.restart = 0, start = beta.pow_25[,i-1])
  
  while (sum(betacurr == 0) == sum(beta.pow_25[,i-1] == 0)) {
    lamnew <- exp(log(lamnew) - epsilon)
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = n, lambda = lamnew, q = q, 
                     rand.restart = 0, start = betacurr)
  }
  if (sum(betacurr == 0) != (sum(beta.pow_25[,i-1] == 0) - 1)) {
    lower <- lamnew
    upper <- lamseq[i - 1]
  } else {
    lower <- upper <- lamnew
  }
  while (sum(betacurr == 0) != (sum(beta.pow_25[,i-1] == 0) - 1) & round(lower, 8) != round(upper, 8)) {
    # cat("lower=", lower, "\n")
    # cat("upper=", upper, "\n")
    mid <- (lower + upper)/2
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = n, lambda = mid, q = q, 
                     rand.restart = 0, start = beta.pow_25[,i-1])
    # cat("nz=", sum(betacurr == 0), "\n")
    if ((sum(betacurr == 0) < (sum(beta.pow_25[,i-1] == 0) - 1))) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  lamseq[i] <- lamnew
  beta.pow_25[,i] <- betacurr
  nz[i] <- sum(betacurr == 0)
  cat("nz=", nz[i], "\n")
  # plot(lamseq[1:i], nz[1:i])
  if (nz[i] == 0) {
    break
  }
}

par(mfrow = c(1, 2))
plot(lamseq, nz)
abline(h = 0:30, lty = 3)
plot(log(lamseq))
