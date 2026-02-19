## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
# library(devtools)
# install_github(repo="ryantibs/best-subset", subdir="bestsubset")
library(bestsubset)
glmnet.control(fdev=0)
# library(devtools)
# Install the full XCode application from App Store. Command line tools are not enough.
# Do not uninstall XCode!!! It will mess up your R package installation.
# install_github("maryclare/powopt")

library(BeSS)
# library(powopt)
library(powoptRcpp)
rm(list=ls())


## -----------------------------------------------------------------------------
n = 70; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set

# nrep should be 500, but we use 10 for speed
nrep = 5 # Number of repetitions

# nrep = 500 # Number of repetitions
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type


## -----------------------------------------------------------------------------
set.seed(seed)
xy.obj = sim.xy(n,p,nval,rho=0.35,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
# y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 300
# nrel = 9

ip.las = matrix(0,nrep,p+1)
ip.fs = matrix(0,nrep,p+1)
ip.bs = matrix(0,nrep,p+1)
ip.bs2 = matrix(0,nrep,p+1)
ip.pow_25 = matrix(0,nrep,p+1)
ip.pow_50 = matrix(0,nrep,p+1)
ip.pow_75 = matrix(0,nrep,p+1)


for (r in 1:nrep) {
  cat(r,"... ")
  eps = rnorm(n)*sigma
  y = mu + eps
  
  fit.lasso = lasso(x,y,intercept=FALSE,nlam=nlam,nrel=1) # we don't use nrel = nrel (relaxed lasso) as in fig.df.R, instead we use nrel = 1 which means the regular lasso
  beta.las = as.matrix(coef(fit.lasso))
  nzs.las = colSums(beta.las != 0)
  j = nlam - rev(match(p:0, rev(nzs.las), NA))
  ind = j + 1
  yhat.las = (x %*% beta.las)[,ind]
  ip.las[r,] = colSums(yhat.las * eps)
  
  # Create lambdas for the power penalty
  # n <- 70  # example number of observations
  # alpha <- 1
  
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

  
  # q = 0.25
  q = 0.25
  
  lammax_25 <- max((2*(1 - q))^(1 - q)*(2-q)^(q-2)*abs(crossprod(xstd, ystd))^(2-q)*diag(crossprod(xstd))^(q - 1))
# eps <- 10^(-4)
# eps = lam*sum(|b|^q)
# eps/sum(|b|^q) = lam

# lammin <- eps/sum(abs(bsv)^q)
lammin_25 <- 10^(-4)/sum(abs(bm)^q)

mylams_25 <- exp(seq(log(lammax_25), 
                  log(lammin_25),
                  length.out = nlam))

  betas_25 <- matrix(nrow = ncol(xstd), ncol = length(mylams_25))
start_25 <- rep(0, ncol(xstd))
for (k in 1:ncol(betas_25)) {
  myest_25 <- powCD(y = ystd, X = xstd, q = q, sigma.sq = n, 
                 lambda = mylams_25[k]/n, tol = 10^(-14),
                 start = start_25)
  start_25 <- myest_25
  betas_25[, k] <- myest_25
}
obetas_25 <- sweep(betas_25, 1, xscales, "/")*yscale

  nzs.pow_25 = colSums(obetas_25 != 0)
  k_25 = nlam - rev(match(p:0, rev(nzs.pow_25), NA))
  ind.pow_25 = k_25 + 1
  yhat.pow_25 = (x %*% obetas_25)[,ind.pow_25]
  
  
  ip.pow_25[r,] = colSums(yhat.pow_25 * eps)
  
  
  # q = 0.50
  q = 0.5
  
  lammax_50 <- max((2*(1 - q))^(1 - q)*(2-q)^(q-2)*abs(crossprod(xstd, ystd))^(2-q)*diag(crossprod(xstd))^(q - 1))
# eps <- 10^(-4)
# eps = lam*sum(|b|^q)
# eps/sum(|b|^q) = lam

# lammin <- eps/sum(abs(bsv)^q)
lammin_50 <- 10^(-4)/sum(abs(bm)^q)

mylams_50 <- exp(seq(log(lammax_50), 
                  log(lammin_50),
                  length.out = nlam))

  betas_50 <- matrix(nrow = ncol(xstd), ncol = length(mylams_50))
start_50 <- rep(0, ncol(xstd))
for (k in 1:ncol(betas_50)) {
  myest_50 <- powCD(y = ystd, X = xstd, q = q, sigma.sq = n, 
                 lambda = mylams_50[k]/n, tol = 10^(-14),
                 start = start_50)
  start_50 <- myest_50
  betas_50[, k] <- myest_50
}
obetas_50 <- sweep(betas_50, 1, xscales, "/")*yscale

  nzs.pow_50 = colSums(obetas_50 != 0)
  k_50 = nlam - rev(match(p:0, rev(nzs.pow_50), NA))
  ind.pow_50 = k_50 + 1
  yhat.pow_50 = (x %*% obetas_50)[,ind.pow_50]
  
  
  ip.pow_50[r,] = colSums(yhat.pow_50 * eps)
  
  
  
  # q = 0.75
  q = 0.75
  
  lammax_75 <- max((2*(1 - q))^(1 - q)*(2-q)^(q-2)*abs(crossprod(xstd, ystd))^(2-q)*diag(crossprod(xstd))^(q - 1))
# eps <- 10^(-4)
# eps = lam*sum(|b|^q)
# eps/sum(|b|^q) = lam

# lammin <- eps/sum(abs(bsv)^q)
lammin_75 <- 10^(-4)/sum(abs(bm)^q)

mylams_75 <- exp(seq(log(lammax_75), 
                  log(lammin_75),
                  length.out = nlam))

  betas_75 <- matrix(nrow = ncol(xstd), ncol = length(mylams_75))
start_75 <- rep(0, ncol(xstd))
for (k in 1:ncol(betas_75)) {
  myest_75 <- powCD(y = ystd, X = xstd, q = q, sigma.sq = n, 
                 lambda = mylams_75[k]/n, tol = 10^(-14),
                 start = start_75)
  start_75 <- myest_75
  betas_75[, k] <- myest_75
}
obetas_75 <- sweep(betas_75, 1, xscales, "/")*yscale

  nzs.pow_75 = colSums(obetas_75 != 0)
  k_75 = nlam - rev(match(p:0, rev(nzs.pow_75), NA))
  ind.pow_75 = k_75 + 1
  yhat.pow_75 = (x %*% obetas_75)[,ind.pow_75]
  
  
  ip.pow_75[r,] = colSums(yhat.pow_75 * eps)
  
  
  
  
  
  yhat.fs = predict(fs(x,y,intercept=FALSE))[, 1:31] # At some point we should figure out why it's giving us 32
  yhat.bs = bestsubset::predict.bs(bestsubset::bs(x,y,intercept=FALSE))
  fit_bess = bess(x,y, method = "sequential", s.list= 1:min(n,p))
  temp = predict(fit_bess,newdata = x)
  temp2 = cbind(0,t(temp))
  yhat.bs2 = temp2
  
# intercept?
  ip.fs[r,] = colSums(yhat.fs * eps)
  ip.bs[r,] = colSums(yhat.bs * eps)
  ip.bs2[r,] = colSums(yhat.bs2 * eps)
}

df.pow_25 = colMeans(ip.pow_25, na.rm=TRUE) / sigma^2
df.pow_50 = colMeans(ip.pow_50, na.rm=TRUE) / sigma^2
df.pow_75 = colMeans(ip.pow_75, na.rm=TRUE) / sigma^2

df.las = colMeans(ip.las, na.rm=TRUE) / sigma^2
# df.las = matrix(df.las, p+1, nrel, byrow=TRUE)
df.fs = colMeans(ip.fs, na.rm=TRUE) / sigma^2
df.bs = colMeans(ip.bs, na.rm=TRUE) / sigma^2
df.bs2 = colMeans(ip.bs2, na.rm=TRUE) / sigma^2

save(list=ls(),file="bess_sim.df_new_epsilon_June.rda")


## -----------------------------------------------------------------------------
# Run the code below to reproduce the df figure without rerunning the sims
library(bestsubset)
rm(list = ls())
load("./bess_sim.df_new_epsilon_June.rda")
# Plot the results
dat = data.frame(x=rep(0:p,7),
                 y=c(df.bs, df.bs2, df.fs,df.las,df.pow_25,df.pow_50,df.pow_75),
                 Method=factor(rep(c("Forward stepwise","Lasso","Power Penalty q = 0.25","Power Penalty q = 0.5", "Power Penalty q = 0.75"),
                                    rep(p+1,5))))
p <- ggplot(dat, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(pch=19) +
    theme_bw() + 
  theme(
    legend.position = c(0.95, 0.05),              # inside panel (x=right, y=bottom)
    legend.justification = c("right", "bottom"), # anchor point
    legend.background = element_rect(fill = "white", color = NA)
  )
print(p)
ggsave("bess_df4_June.pdf", height=6, width=8, device="pdf")


## -----------------------------------------------------------------------------
# Create the ggplot
dat_las <- data.frame(lambda = fit.lasso$lambda, nzs = nzs.las)
ggplot(dat_las, aes(x = log(lambda), y = nzs)) +
  geom_point() +  # Add points
  geom_line() +   # Add lines
  labs(
    title = "Lambda against the # of nonzero coefficients for the lasso",
    x = "log(lambda)",
    y = "Number of nonzero coefficients"
  ) +
  theme_minimal()
ggsave("lambda_nzs_lasso.pdf", height=6, width=8, device="pdf")


## -----------------------------------------------------------------------------
# Create the ggplot
dat25 <- data.frame(lambda = lambda_seq_manual_25, nzs = nzs.pow_25)
ggplot(dat25, aes(x = log(lambda), y = nzs)) +
  geom_point() +  # Add points
  geom_line() +   # Add lines
  labs(
    title = "Lambda vs. the # of nonzero coefficients for power penalty regression, q = 0.25",
    x = "log(lambda)",
    y = "Number of nonzero coefficients"
  ) +
  theme_minimal()
ggsave("lambda_nzs_25_epsilon.pdf", height=6, width=8, device="pdf")


## -----------------------------------------------------------------------------
# Create the ggplot
dat50 <- data.frame(lambda = lambda_seq_manual_50, nzs = nzs.pow_50)
ggplot(dat50, aes(x = log(lambda), y = nzs)) +
  geom_point() +  # Add points
  geom_line() +   # Add lines
  labs(
    title = "Lambda vs. the # of nonzero coefficients for power penalty regression, q = 0.5",
    x = "log(lambda)",
    y = "Number of nonzero coefficients"
  ) +
  theme_minimal()
ggsave("lambda_nzs_50_epsilon.pdf", height=6, width=8, device="pdf")


## -----------------------------------------------------------------------------
# Create the ggplot
dat75 <- data.frame(lambda = lambda_seq_manual_75, nzs = nzs.pow_75)
ggplot(dat75, aes(x = log(lambda), y = nzs)) +
  geom_point() +  # Add points
  geom_line() +   # Add lines
  labs(
    title = "Lambda vs. the # of nonzero coefficients for power penalty regression, q = 0.75",
    x = "log(lambda)",
    y = "Number of nonzero coefficients"
  ) +
  theme_minimal()
ggsave("lambda_nzs_75.pdf", height=6, width=8, device="pdf")


## -----------------------------------------------------------------------------
rm(list=ls())
load("./sim.df_new_epsilon.rda")
dat_las <- data.frame(lambda = fit.lasso$lambda, nzs = nzs.las)
dat25 <- data.frame(lambda = lambda_seq_manual_25, nzs = nzs.pow_25)
dat50 <- data.frame(lambda = lambda_seq_manual_50, nzs = nzs.pow_50)
dat75 <- data.frame(lambda = lambda_seq_manual_75, nzs = nzs.pow_75)
# Add a group identifier to each data frame
dat_las$group <- "Lasso"
dat25$group <- "q = 0.25"
dat50$group <- "q = 0.50"
dat75$group <- "q = 0.75"

# Combine all data frames into one
combined_df <- rbind(dat_las, dat25, dat50, dat75)

# Create the ggplot
ggplot(combined_df, aes(x = log(lambda), y = nzs, color = group)) +
  geom_line() +  # Add lines
  geom_point() + # Add points for emphasis
  labs(
    title = "Comparison of lambda vs the number of nonzero coefficients",
    x = "log(lambda)",
    y = "Number of nonzero coefficients",
    color = "Method"
  ) +
  theme_minimal()
# ggsave("lambda_nzs_all_epsilon.pdf", height=6, width=8, device="pdf")

