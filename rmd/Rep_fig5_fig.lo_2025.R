## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
# library(devtools)
# install_github(repo="ryantibs/best-subset", subdir="bestsubset")
library(bestsubset)
glmnet.control(fdev=0)
# library(devtools)
# install_github("maryclare/powopt")
library(powopt)
library(L0Learn)
#install.packages("slam")
#install.packages("/Library/gurobi1103/macos_universal2/R/gurobi_11.0-3_R_4.4.0.tgz", 
#                 repos = NULL, 
#                 type = "binary")
rm(list=ls())


## -----------------------------------------------------------------------------
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
                 lambda = mylams[k]/n, tol = 10^(-14),
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


## -----------------------------------------------------------------------------
# Set some overall simulation parameters
n = 50; p = 50 # Size of training set, and number of predictors
# note: if n=50 and p=50, change the tolerance to 1e-7 in powcd.
nval = n # Size of validation set
nrep = 10 # Number of repetitions for a given setting
# nrep=100 for n=10 or 50, p = 50
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
#reg.funs[["L0Learn 1"]] = function(x,y)L0Learn.fit(x,y,penalty="L0",algorithm="CDPSI",nLambda=50,intercept=FALSE)
#reg.funs[["L0Learn 2"]] = function(x,y) L0Learn.fit(x,y,penalty="L0L1",algorithm="CDPSI",nGamma=10,nLambda=50, intercept=FALSE)


## -----------------------------------------------------------------------------
file.list = c() # Vector of files for the saved rds files
for (beta.type in type.vec) {
  for (rho in rho.vec) {
    name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
    for (snr in snr.vec) {
      file = paste0("Rep_rds_new/", name, ".snr", round(snr,2), ".rds")
      cat("..... NEW SIMULATION .....\n")
      cat("--------------------------\n")
      cat(paste0("File: ", file, "\n\n"))

      sim.master(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=seed, s=s,
                 verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)

      file.list = c(file.list, file)
      cat("\n")
    }
  }
}


## -----------------------------------------------------------------------------
# From sim.lo.R
# Run the code below to reproduce the figures without rerunning the sims
library(bestsubset)
n = 100; p = 10
# file.list = system(paste0("ls Rep_rds/sim.n",n,".p",p,".*.rds"),intern=TRUE)
file.list = system(paste0("ls Rep_rds_new/sim.n", n, ".p", p, ".beta2.rho0.35.*.rds"), intern = TRUE)
method.nums = c(2,1,3,4)
method.names = c("Forward stepwise","Lasso","Relaxed lasso","PowCD")

# Validation tuning

# From fig.lo.R
plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="~/Downloads/Research_Spring_2025",
               file.name="lo.err_new", w = 8, h = 6)

plot.from.file(file.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="~/Downloads/Research_Spring_2025",
               file.name="lo.prop_new", w = 8, h = 6)

plot.from.file(file.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="~/Downloads/Research_Spring_2025",
               file.name="lo.nzs_new", w = 8, h = 6)

plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               legend.pos="bottom", make.pdf=TRUE, fig.dir="~/Downloads/Research_Spring_2025",
               file.name="lo.F_new", w = 8, h = 6)



## -----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
n = 10; p = 50; s = 5
file.list = system(paste0("ls Rep_rds_new/sim.n", n, ".p", p, ".beta2.rho0.35.*.rds"), intern = TRUE)
method.nums = c(2,1,3,4)
method.names = c("Forward stepwise","Lasso","Relaxed lasso","PowCD")
# Create individual plots for each metric
plot_error <- plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none")
plot_prop <- plot.from.file(file.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none")
plot_nzs <- plot.from.file(file.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none")
plot_F <- plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none")
# Arrange the plots in a 2x2 grid
grid_arranged <- (plot_error + plot_prop + plot_nzs + plot_F) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")
# Save the combined plot as a PDF
beta <- 2; rho <- 0.35
outfile <- sprintf("~/Downloads/Research_Spring_2025/fig5_supp_plots/lo.sim.n%d.p%d.beta%d.rho%.2f.pdf",
                   n, p, beta, rho)
ggsave(filename = outfile, plot = grid_arranged, width = 10, height = 8)


## -----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
n = 10; p = 50; s = 5
file.list = system(paste0("ls Rep_rds_lo_2025_new/sim_n10_p50/sim.n", n, ".p", p, ".beta2.rho0.35.*.rds"), intern = TRUE)
method.nums = c(2,1,3,4)
method.names = c("Forward stepwise","Lasso","Relaxed lasso","PowCD")
# Create individual plots for each metric
plot_error <- plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank())
plot_F <- plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE, w = 8, h = 6) + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank())
# Arrange the plots in a 2x2 grid
grid_arranged <- (plot_error + plot_F) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")
# Save the combined plot as a PDF
beta <- 2; rho <- 0.35
outfile <- sprintf("~/Downloads/Research_Spring_2025/fig5_supp_plots/lo.poster.pdf",
                   n, p, beta, rho)
ggsave(filename = outfile, plot = grid_arranged, width = 10, height = 4)


## -----------------------------------------------------------------------------
# From sim.lo.R
# Run the code below to reproduce the figures without rerunning the sims

library(bestsubset)
n = 50; p = 50; s = 5
file.list = system(paste0("ls Rep_rds_lo_2025_new/sim_n50_p50/sim.n",n,".p",p,".*.rds"),intern=TRUE)
# method.nums = c(3,2,1,4,5,6,7)
# method.names = c("Best subset","Forward stepwise","Lasso","PowCD", "Relaxed lasso",
#                  "L0Learn 1", "L0Learn 2")

method.nums = c(2,1,3,4)
method.names = c("Forward stepwise","Lasso","Relaxed lasso","PowCD")

# Validation tuning
plot.from.file(file.list, what="risk", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.risk.rel"))

plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.err.rel"))

plot.from.file(file.list, what="prop", tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.prop"))

plot.from.file(file.list, what="nonzero", tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.nzs"))

plot.from.file(file.list, what="F", tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.F"))

plot.from.file(file.list, what=c("error", "risk", "prop", "F", "nonzero"), tuning="val",
               method.nums=method.nums, method.names=method.names,
               main=paste0("n=",n,", p=",p,", s=",s), make.pdf=TRUE,
               fig.dir="~/Downloads/Research_Spring_2025/fig5_supp_plots",
               file.name=paste0("sim.n",n,".p",p,".val.all"))




## -----------------------------------------------------------------------------
library(gridExtra)
library(ggplot2)
n = 50; p = 50; s = 5
file.list = system(paste0("ls Rep_rds_lo_2025_new/sim_n50_p50/sim.n",n,".p",p,".*.rds"),intern=TRUE)
method.nums = c(2,1,3,4)
method.names = c("Forward stepwise","Lasso","Relaxed lasso","PowCD")
# Create individual plots for each metric
plot_error <- plot.from.file(file.list, what="error", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE) + theme(legend.position = "none")
plot_prop <- plot.from.file(file.list, what="prop", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE) + theme(legend.position = "none")
plot_nzs <- plot.from.file(file.list, what="nonzero", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE) + theme(legend.position = "none")
plot_F <- plot.from.file(file.list, what="F", rel.to=NULL, tuning="val",
               method.nums=method.nums, method.names=method.names, make.pdf=FALSE) + theme(legend.position = "none")
# Arrange the plots in a 2x2 grid
grid_arranged <- (plot_error + plot_prop + plot_nzs + plot_F) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")
# Save the combined plot as a PDF
outfile <- sprintf("~/Downloads/Research_Spring_2025/fig5_supp_plots/lo.sim.n%d.p%d.pdf",
                   n, p, beta, rho)
ggsave(filename = outfile, plot = grid_arranged, width = 10, height = 8)

