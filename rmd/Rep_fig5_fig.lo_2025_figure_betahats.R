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
n = 50; p = 100 # Size of training set, and number of predictors
# high-5 setting
# n = 50; p = 1000
nval = n # Size of validation set
nrep = 1 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
type.vec = c(1:3,5) # Simulation settings to consider
rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios

stem = paste0("sim.n",n,".p",p)

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
#reg.funs[["Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,nlam=50)
#reg.funs[["Forward stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
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
sim.master.betahat <- function (n, p, nval, reg.funs, nrep = 50, seed = NULL, verbose = FALSE, 
          file = NULL, rho = 0, s = 5, beta.type = 1, 
          snr = 1) 
{
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  N = length(reg.funs)
  reg.names = names(reg.funs); if (is.null(reg.names)) reg.names = paste("Method", 1:N)

  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = opt = runtime = vector("list", N)
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) =
    names(F1) = names(opt) = names(runtime) = reg.names

  ## container for final chosen betas (matches tun.val: row-wise min of err.val)
  # Each entry will be a matrix [nrep x ncoef]; columns include intercept name if present.
  beta.hat.final = vector("list", N)
  names(beta.hat.final) = reg.names
  ##

  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA, nrep, 1)
  }
  filled = rep(FALSE, N)

  err.null = risk.null = sigma = rep(NA, nrep)

  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n", i, nrep))
      cat("  Generating data ...\n")
    }
    xy.obj = sim.xy(n, p, nval, rho, s, beta.type, snr)
    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma

    for (j in 1:N) {
      if (verbose) cat(sprintf("  Applying regression method %i (of %i) ...\n", j, N))
      tryCatch({
        runtime[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xy.obj$x, xy.obj$y)
        })[1]

        coef.raw = as.matrix(coef(reg.obj))   # [ncoef x m], may include "(Intercept)"
        m        = ncol(coef.raw)
        ncoef    = nrow(coef.raw)
        rnames   = rownames(coef.raw)

        muhat.train = as.matrix(predict(reg.obj, xy.obj$x))
        muhat.val   = as.matrix(predict(reg.obj, xy.obj$xval))

        if (!filled[j]) {
          # expand metric matrices to width m once
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
            nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA, nrep, m)
          ## allocate beta matrix for this method
          beta.hat.final[[j]] = matrix(NA_real_, nrep, ncoef,
                                       dimnames = list(NULL, rnames))
          ##
          filled[j] = TRUE
        }

        err.train[[j]][i, ] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i, ]   = colMeans((muhat.val   - xy.obj$yval)^2)

        # For risk/selection counts, drop intercept if present
        betahat = coef.raw
        if (nrow(betahat) == p + 1) {
          intercept = TRUE
          betahat0  = betahat[1, ]
          betahat   = betahat[-1, , drop = FALSE]
        } else intercept = FALSE
        
        # View(betahat)

        delta = betahat - xy.obj$beta
        risk[[j]][i, ] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i, ] = risk[[j]][i, ] + betahat0^2

        err.test[[j]][i, ] = risk[[j]][i, ] + xy.obj$sigma^2
        prop[[j]][i, ]     = 1 - err.test[[j]][i, ]/err.null[i]

        nzs[[j]][i, ]  = colSums(betahat != 0)
        tpos           = colSums((betahat != 0) * (xy.obj$beta != 0))
        fpos[[j]][i, ] = nzs[[j]][i, ] - tpos
        fneg[[j]][i, ] = colSums((betahat == 0) * (xy.obj$beta != 0))
        F1[[j]][i, ]   = 2 * tpos/(2 * tpos + fpos[[j]][i, ] + fneg[[j]][i, ])
        opt[[j]][i, ]  = (err.test[[j]][i, ] - err.train[[j]][i, ]) / err.train[[j]][i, ]

        
        #cat(coef.raw, "\n")
        ## save the chosen column's coefficients (same rule as tun.val)
        # choose per-replicate index by row-wise min of validation error
        chosen.idx = which.min(err.val[[j]][i, ])
        if (!is.na(chosen.idx) && chosen.idx <= ncol(coef.raw)) {
          beta.hat.final[[j]][i, ] = coef.raw[, chosen.idx]  # includes intercept if present
          cat(beta.hat.final[[j]], "\n")
        } else {
          # leave NA row if no valid choice
          beta.hat.final[[j]][i, ] = NA_real_
        }
        ##
        #cat(beta.hat.final,"\n")

      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n", err$message))
          cat("    *** End error message ***\n")
        }
        ## on error, ensure a row of NAs exists so dimensions stay consistent
        if (!is.null(beta.hat.final[[j]])) beta.hat.final[[j]][i, ] = NA_real_
        ##
      })
    }
    

    # # periodic checkpoint also writes beta.hat.final
    # if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
    #   saveRDS(bestsubset:::enlist(err.train, err.val, err.test, err.null,
    #                  prop, risk, risk.null, nzs, fpos, fneg, F1,
    #                  opt, sigma, runtime, beta.hat.final),
    #           file = file)
    # }
  } # enlist is a hidden function in bestsubset package

  saveRDS(beta.hat.final,file = file)
  
  # out = bestsubset:::enlist(err.train, err.val, err.test, err.null, prop,
  #              risk, risk.null, nzs, fpos, fneg, F1, opt, sigma, runtime, beta.hat.final)
  # # enlist is a hidden function in bestsubset package
  # 
  # if (!is.null(file)) saveRDS(out, file)
  # out = bestsubset:::choose.tuning.params(out)  # still returns tun.val & tun.ora
  # # choose.tuning.params is a hidden function in bestsubset package
  # out = c(out, list(rho = rho, s = s, beta.type = beta.type, snr = snr, call = this.call))
  # class(out) = "sim"
  # if (!is.null(file)) { saveRDS(out, file); invisible(out) } else return(out)
}


## -----------------------------------------------------------------------------
#type.vec = c(1:3,5) # Simulation settings to consider
#rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
#snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios
type.vec = 3
rho.vec = 0
snr.vec = 1.22

file.list = c() # Vector of files for the saved rds files
for (beta.type in type.vec) {
  for (rho in rho.vec) {
    name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
    for (snr in snr.vec) {
      file = paste0("rds_betahat/", name, ".snr", round(snr,2), ".rds")
      cat("..... NEW SIMULATION .....\n")
      cat("--------------------------\n")
      cat(paste0("File: ", file, "\n\n"))
      # insert a timer here
      start_time <- Sys.time()
      cat("Start time: ", start_time, "\n")
      sim.master.betahat(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=1, s=s,
                 verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)
      end_time <- Sys.time()
      cat("End time: ", end_time, "\n")
      cat("Time taken: ", end_time - start_time, "\n")
      file.list = c(file.list, file)
      cat("\n")
    }
  }
}


## -----------------------------------------------------------------------------
sim.obj = readRDS(file.list[1])
# View(sim.obj$beta.hat.final$PowCD)
b.relax <- sim.obj$`Relaxed lasso`[1, ]
b.pow <- sim.obj$`PowCD`[1,]
cat(b.relax, "\n")
cat(b.pow, "\n")
# plot by index and by value on the same graph using ggplot2
library(ggplot2)
ggplot() +
  geom_point(aes(x = 1:length(b.relax), y = b.relax, color = "Relaxed lasso"), size = 1, shape = 3) +
  geom_point(aes(x = 1:length(b.pow), y = b.pow, color = "PowCD"), size = 1, shape = 4) +
  labs(title = sprintf(
    "Comparison of Coefficients from Relaxed Lasso and PowCD\nn=%d, p=%d, s=%d\nbeta.type=%d, rho=%.2f, snr=%.2f",
    n, p, s, beta.type, rho, snr
  ),
       x = "Coefficient Index",
       y = "Coefficient Value") +
  scale_color_manual(values = c("Relaxed lasso" = "blue", "PowCD" = "black")) +
  theme_minimal() +
  theme(legend.title = element_blank())
ggsave(sprintf("fig5_supp_plots/sim.n%d.p%d.beta%d.rho%.2f.snr%.2f.betahat_comparison.pdf", n, p, beta.type, rho, snr),
       width = 8, height = 6)

