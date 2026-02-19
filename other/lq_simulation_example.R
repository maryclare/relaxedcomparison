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
lambda.min.ratio <- 0.0001 
nlambda <- 50

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
q <- 0.5
lammax <- max((2*(1 - q))^(1 - q)*(2-q)^(q-2)*abs(crossprod(xstd, ystd))^(2-q)*diag(crossprod(xstd))^(q - 1))
eps <- 10^(-7)
# eps = lam*sum(|b|^q)
# eps/sum(|b|^q) = lam

# lammin <- eps/sum(abs(bsv)^q)
lammin <- eps/sum(abs(bm)^q)

mylams <- exp(seq(log(lammax), 
                  log(lammin),
                  length.out = 50))

betas <- matrix(nrow = ncol(xstd), ncol = length(mylams))
start <- rep(0, ncol(xstd))
for (k in 1:ncol(betas)) {
  myest <- powCD(y = ystd, X = xstd, q = q, sigma.sq = n, 
                 lambda = mylams[k]/n, tol = 10^(-14),
                 start = start)
  start <- myest
  betas[, k] <- myest
}
obetas <- sweep(betas, 1, xscales, "/")*yscale

plot(log(mylams), betas[1, ], ylim = range(betas), type = "n")
for (k in 1:nrow(betas)) {
  lines(log(mylams), betas[k, ], col = k)
}
plot(apply(betas == 0, 2, sum)/p, ylim = c(0, 1))
range(apply(betas == 0, 2, sum)/p)
# solving
# -2 x'y beta + beta^2 x'x + lam*|beta|^q
# equiv to
# -2 (x'y/ x'x) beta + beta^2 + (lam/(x'x))*|beta|^q

# can find oo s.t. we get a zero
# -2 (x'y/ x'x) beta + beta^2 + (oo^(2-q)/q)*|beta|^q
# then solve for (oo^(2-q)/q) = lam/(x'x)
# oo = (2*(1 - q))^((1 - q)/(2 - q))*q^(1/(2 - q))*(|x'y|/ x'x)/(2-q)
# (oo^(2-q)/q) = (2*(1 - q))^(1 - q)*((|x'y|/ x'x)/(2-q))^(2-q)
#              = (2*(1 - q))^(1 - q)*(2-q)^(q-2)*|x'y|^(2-q)*(x'x)^(q - 2)
# (oo^(2-q)/q)*(x'x) = (2*(1 - q))^(1 - q)*(2-q)^(q-2)*|x'y|^(2-q)*(x'x)^(q - 1)
# (2*(1 - q))^(1 - q)*(2-q)^(q-2)*|x'y|^(2-q)*(x'x)^(q - 1)