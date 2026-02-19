
library(powopt)
library(tidyverse)

n <- 200
p <- 50
prop <- 0.2
beta <- numeric(p)
beta[1:(prop * p)] <- 1

rho <- 0.4
Sigma <- toeplitz(x = rho^(0:(p-1)))
R <- chol(Sigma)

X <- matrix(rnorm(n*p, 0, 0.2), n, p) %*% R
colnorms <- apply(X, 2, function(t) sqrt(sum(t^2)))
x_standardized <- t(t(X) / colnorms)


y <- X %*% beta + rnorm(n, 0, 1)

x_y_dot_standardized <- t(x_standardized) %*% y

# Find the maximum absolute value of <x_j, y> on standardized x
max_abs_x_y_dot_standardized <- max(abs(x_y_dot_standardized))

q = 0.5
a <- (2-q)*q^(1/(q-2)) * (2*(1-q))^((q-1)/(2-q))

# Calculate lambda_max
omega <- max_abs_x_y_dot_standardized / a
lambda_max_standardized <- omega^(2-q) / q

nlam = p + 1
lamseq <- c(lambda_max_standardized, rep(NA, nlam - 1))
nz <- c(p, rep(NA, nlam - 1))

beta.pow = matrix(0,p, nlam)
beta.pow[,1] = powCD(X = x_standardized, y = y, 
                        sigma.sq = n, lambda = lamseq[1], 
                        q = q, rand.restart = 0, start = rep(0,p))


epsilon <- 0.001
# lambda min
beta_hat <- lm(y ~ x_standardized - 1)$coef
lambda.min <- epsilon / sum(abs(beta_hat)^q)
# smallest increment in log-lambda if we can have at most 10^5 lambda values 
delta <- (log(lambda_max_standardized) - log(lambda.min)) / 10^5

b <- 1
for (i in 2:nlam) {
  cat("i=", i, "\n")
  # calculate the next value where we expect to see a new non-zero
  res <- y - x_standardized %*% beta.pow[,i-1]
  x_res_dot_standardized <- t(x_standardized) %*% res
  max_abs_x_res_dot_standardized <- max(abs(x_res_dot_standardized)[which(beta.pow[,i-1] == 0)])
  omega <- max_abs_x_res_dot_standardized / a
  est_lambda <- omega^(2-q) / q
  
  # next lambda
  increment <- (log(lamseq[i-1]) - log(est_lambda)) 
  if(increment < delta) increment <- delta 
    
  lamnew <- exp(log(lamseq[i - 1]) - increment)
  betacurr <- powCD(X = x_standardized, y = y, 
                    sigma.sq = 1, lambda = lamnew, q = q, # I dont understand what sigma.sq mean, is this the variance of Y?
                    rand.restart = 0, start = beta.pow[,i-1])
  b <- b + 1
  
  while (sum(betacurr == 0) == sum(beta.pow[,i-1] == 0)) {
    # calculate the next value where we expect to see a new non-zero
    lambda_curr <- lamnew
    res <- y - x_standardized %*% betacurr
    x_res_dot_standardized <- t(x_standardized) %*% res
    max_abs_x_res_dot_standardized <- max(abs(x_res_dot_standardized)[which(beta.pow[,i-1] == 0)])
    omega <- max_abs_x_res_dot_standardized / a
    est_lambda <- omega^(2-q) / q
    
    # next lambda
    increment <- (log(lambda_curr) - log(est_lambda)) 
    # cat(increment, "'\n")
    if(increment < delta) increment <- delta 
    
    lamnew <- exp(log(lambda_curr) - increment)
    
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = 1, lambda = lamnew, q = q, 
                     rand.restart = 0, start = beta.pow[,i-1])
    b <- b+1
  }
  
  if (sum(betacurr == 0) != (sum(beta.pow[,i-1] == 0) - 1)) {
    lower <- lamnew
    upper <- lambda_curr
  } else {
    lower <- upper <- lamnew
  }
  while (sum(betacurr == 0) != (sum(beta.pow[,i-1] == 0) - 1) & round(lower, 8) != round(upper, 8)) {
    # cat("lower=", lower, "\n")
    # cat("upper=", upper, "\n")
    mid <- (lower + upper)/2
    betacurr = powCD(X = x_standardized, y = y, 
                     sigma.sq = 1, lambda = mid, q = q, 
                     rand.restart = 0, start = beta.pow[,i-1])
    b <- b+1
    # cat("nz=", sum(betacurr == 0), "\n")
    if ((sum(betacurr == 0) < (sum(beta.pow[,i-1] == 0) - 1))) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  lamseq[i] <- lamnew
  beta.pow[,i] <- betacurr
  nz[i] <- sum(betacurr == 0)
  cat("nz=", nz[i], "\n")
  # plot(lamseq[1:i], nz[1:i])
  if (nz[i] == 0) {
    break
  }
}

cat("Num of evaluations = ", b , "\n")

par(mfrow = c(1, 2))
plot(log(lamseq), colSums(beta.pow!=0))
abline(h = 0:p, lty = 3)
plot(log(lamseq))

