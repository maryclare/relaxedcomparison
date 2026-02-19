# library(devtools)
# install_github(repo="ryantibs/best-subset", subdir="bestsubset")
# library(bestsubset)
library(colorspace)
library(glmnet)
glmnet.control(fdev=0)
# library(devtools)
# install_github("maryclare/powopt")
library(powopt)
library(L0Learn)
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/check.R")
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/sim.R")
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/common.R")
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/fs.R")
source("https://raw.githubusercontent.com/ryantibs/best-subset/refs/heads/master/bestsubset/R/lasso.R")
#install.packages("slam")
#install.packages("/Library/gurobi1103/macos_universal2/R/gurobi_11.0-3_R_4.4.0.tgz", 
#                 repos = NULL, 
#                 type = "binary")

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
    
    eps <- 10^(-4)
    if (q_val <= 1) {
      lammax <- max((2*(1 - q_val))^(1 - q_val)*(2-q_val)^(q_val-2)*abs(crossprod(xstd, ystd))^(2-q_val)*diag(crossprod(xstd))^(q_val - 1))
    } else {
      lammax <- 1/eps - min(sv$d)
    }
    
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
                     lambda = mylams[k]/n, tol = 10^(-7),
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
  
  # saveRDS(beta.hat.final,file = file)
  return(list("beta.hat.final" = beta.hat.final, "err.val" = err.val, "err.test" = err.test,
              "risk" = risk,
              "fitted.final" = lapply(beta.hat.final, function(b) {xy.obj$xval%*%t(b)}),
              "y" = xy.obj$y, "x" = xy.obj$x,
              "beta" = xy.obj$beta))
  
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

#type.vec = c(1:3,5) # Simulation settings to consider
#rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
#snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios
type.vec = 3
rho.vec = 0.7
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
      res <- sim.master.betahat(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=1, s=s,
                         verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr)
      end_time <- Sys.time()
      cat("End time: ", end_time, "\n")
      cat("Time taken: ", end_time - start_time, "\n")
      file.list = c(file.list, file)
      cat("\n")
    }
  }
}

layout(cbind(c(1, 1), c(2, 3), c(4, 5)))
par(mar = c(2, 2, 0, 0))
par(oma = c(0, 0, 1, 1))
plot(c(res$beta.hat.final$`Relaxed lasso`),
     col = rgb(1, 0, 0, 0.5), type = "h",
     ylim = range(unlist(res$beta.hat.final)))
points(which(c(res$beta.hat.final$`Relaxed lasso`) != 0),
       c(res$beta.hat.final$`Relaxed lasso`)[c(res$beta.hat.final$`Relaxed lasso`) != 0], pch = 16,
       col = rgb(1, 0, 0, 0.5))
points(c(res$beta.hat.final$`PowCD`),
     pch = 16, col = rgb(0, 0, 1, 0.5), type = "h")
points(which(c(res$beta.hat.final$`PowCD`) != 0),
       c(res$beta.hat.final$`PowCD`)[c(res$beta.hat.final$`PowCD`) != 0], pch = 16,
       col = rgb(0, 0, 1, 0.5))

# For relaxed lasso, gammas are changing faster than lambdas, gammas go from 1 to 0,
# gamma of 1 is lasso, gamma of 0 is least squares
# For this simulation it chooses the solution with first three predictors in,
# least squares estimates (no shrinkage)
forplotrl <- matrix(res$err.val$`Relaxed lasso`[1, ]/min(res$err.val$`Relaxed lasso`[1, ]), 10, 50)[10:1, ]
riskrl <- matrix(res$risk$`Relaxed lasso`[1, ]/min(res$risk$`Relaxed lasso`[1, ]), 10, 50)[10:1, ]
forplotpo <- t(matrix(res$err.val$`PowCD`[1, ]/min(res$err.val$`PowCD`[1, ]), 50, 10))
riskpo <- t(matrix(res$risk$`PowCD`[1, ]/min(res$risk$`PowCD`[1, ]), 50, 10))

brsfprl <- exp(seq(0, log(max(c(forplotrl))),
               length.out = 20))
brsfppo <- exp(seq(0, log(max(c(forplotpo))),
                   length.out = 20))
brsrirl <- exp(seq(0, log(max(c(riskrl))),
                 length.out = 20))
brsripo <- exp(seq(0, log(max(c(riskpo))),
                 length.out = 20))
# cols <- rev(sequential_hcl(length(brs) - 1, "Purple"))
cols <- rev(sequential_hcl(length(brsfprl) - 1, "Purple"))


image(1:10, 1:50, forplotrl,
      br = brsfprl, 
      col = cols, 
      axes = FALSE)
mtext(expression(gamma), 1, line = 2)
mtext(expression(lambda), 2)
wmin <- which(forplotrl  == min(forplotrl), arr.ind = TRUE)
for (i in 1:nrow(wmin)) {
  points(wmin[i, 1], wmin[i, 2], pch = 16)
}
box() 
image(1:10, 1:50, riskrl,
      br = brsrirl, col = cols, axes = FALSE)
wmin <- which(riskrl  == min(riskrl), arr.ind = TRUE)
for (i in 1:nrow(wmin)) {
  points(wmin[i, 1], wmin[i, 2], pch = 16)
}

image(1:10, 1:50, forplotpo,
      br = brsfppo, col = cols, axes = FALSE)
box()
axis(1, at = 1:10, lab = round(seq(0, 1, length.out = 10), 1), las = 2)
mtext(expression(q), 1, line = 2)
mtext(expression(paste(lambda, " Index", sep = "")), 2)

wmin <- which(forplotpo  == min(forplotpo), arr.ind = TRUE)
for (i in 1:nrow(wmin)) {
  points(wmin[i, 1], wmin[i, 2], pch = 16)
}

image(1:10, 1:50, riskpo, br = brsripo, col = cols, axes = FALSE)
wmin <- which(riskpo  == min(riskpo), arr.ind = TRUE)
for (i in 1:nrow(wmin)) {
  points(wmin[i, 1], wmin[i, 2], pch = 16)
}

res$beta.hat.final$`Relaxed lasso`
# Looks like the solution is definitely a ls solution (confirmed)

min(res$err.val$`Relaxed lasso`)
min(res$err.val$`PowCD`)

min(res$risk$`Relaxed lasso`)
min(res$risk$`PowCD`)

min(res$err.test$`Relaxed lasso`)
min(res$err.test$`PowCD`)
