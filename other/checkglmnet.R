# Simulate some simple regression data with the first 5 coefficients
# being nonzero
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/sim.R")
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/common.R")

library(glmnet)
library(powopt)

set.seed(0)
n = 100
p = 20
xy.obj = sim.xy(n,p,nval=0,s=5,beta.type=2,snr=1)
x = xy.obj$x
y = xy.obj$y

mean(x)
mean(y)

# I believe that the glmnet options used in the simulations
# in the paper
intercept <- FALSE 
standardize <- TRUE
lambda.min.ratio <- 0.0001 
alpha <- 1
nlambda <- 50
dfmax <- 20
lambda <- NULL

obj = glmnet(x, y, alpha=alpha, nlambda=nlambda, dfmax=dfmax,
             lambda.min.ratio=lambda.min.ratio, lambda=lambda,
             intercept=intercept, standardize=standardize)

# Standardize data (based on glmnet code)
yscale <-sqrt(sum(y^2)) 
ystd <- y/yscale
xscales <- apply(x, 2, function(xx) {
  sqrt(mean((xx - mean(xx))^2))
})
xstd <- sweep(x, 2, xscales, "/")

# Fit to standardized data without standardize option
objstd = glmnet(xstd, ystd, alpha=alpha, nlambda=nlambda, dfmax=dfmax,
             lambda.min.ratio=lambda.min.ratio, lambda=lambda,
             intercept=intercept, standardize=FALSE)

# Both solutions give same deviance ratios suggesting 
# they are equivalent
plot(obj$dev.ratio)
points(objstd$dev.ratio, col = "blue")

# We can relate the lambda sequences
plot(obj$lambda, objstd$lambda*sqrt(sum(y^2)))
abline(a = 0, b = 1)

# Recovered maximum lambda value
max(objstd$lambda) 
max(abs(crossprod(xstd, ystd))/n) 
# We need to divide by n to get there's because they include the 1/n

# Verify that manually constructing the lambda sequence 
# reproduces their results
mylams <- exp(seq(log(max(abs(crossprod(xstd, ystd)))), 
                  log(max(abs(crossprod(xstd, ystd)))*0.0001),
                  length.out = 50))

betas <- matrix(nrow = ncol(xstd), ncol = length(mylams))
start <- rep(0, ncol(xstd))
for (k in 1:ncol(betas)) {
myest <- powCD(y = ystd, X = xstd, q = 1, sigma.sq = n, 
               lambda = mylams[k]/n, tol = 10^(-14),
               start = start)
start <- myest
betas[, k] <- myest
}
plot(c(betas[, 1:length(objstd$lambda)]), c(as.numeric(objstd$beta)))
abline(a = 0, b = 1, col = "blue")
# glmnet automatically stops at some point, not sure exactly when,
# hard to decipher from the code. But that's ok. 
# We don't need to understand exactly when they stop, just
# what sequence they use (which I think we do understand now)

# We can recover the solutions on the original scale
plot(c(as.numeric(obj$beta)),
     c(as.numeric(sweep(betas[, 1:length(objstd$lambda)], 1, xscales, "/")*yscale)))
abline(a = 0, b = 1, col = "blue")
# This means that that our lq functions should standardize 
# the data first, compute lq solutions, then back transform
