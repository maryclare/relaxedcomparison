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
scale_text = paste0("n = ",n, ", p = ",p)
nval = n # Size of validation set

# nrep should be 500, but we use 10 for speed
nrep = 500 # Number of repetitions

# nrep = 500 # Number of repetitions
seed = 1 # Random number generator seed
set.seed(seed)
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type

nlam = 50 # nlam = 31


#############################################################################
# Simulation: sim.xy
rho = 0
snr = 0.7
label = NULL
xy.obj = sim.xy(n,p,nval,rho=0,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma

r <- 1
cat(r,"... ")
eps = rnorm(n)*sigma
y = mu + eps
################################################################################

#############################################################################
# Simulation 2: sim.xy
rho = 0.2
snr = 0.8
label = NULL
xy.obj = sim.xy(n,p,nval,rho=0.2,s=s,beta.type=beta.type,snr=0.8)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma

r <- 1
cat(r,"... ")
eps = rnorm(n)*sigma
y = mu + eps
################################################################################


################################################################################
# Simulation: Orthogonal X
rho = NULL
snr = NULL
sim_text = "Orthogonal X"
label = sim_text
# Generate n x p matrix with small noise, then standardize
X <- matrix(rnorm(n*p, 0, 0.2), n, p)
X <- scale(X)  # Column-wise standardization
u <- svd(X)$u  # Take left singular vectors from SVD
X <- u         # Replace X with orthonormal matrix U (n x n)

# Generate response vector Y ~ N(0,1)
Y <- rnorm(n, 0, 1)
x <- X
y <- Y
##################################################################################


x_standardized <- apply(x, 2, function(z) {(z)/(sqrt(mean(z^2) - mean(z)^2))})

# Calculate <x_j, y> for each j on standardized x
x_y_dot_standardized <- colSums(x_standardized * y)

# Find the maximum absolute value of <x_j, y> on standardized x
max_abs_x_y_dot_standardized <- max(abs(x_y_dot_standardized))

q = 0.25

# Calculate lambda_max
lambda_max_standardized_25 <- (max_abs_x_y_dot_standardized / (2-q) / n) ^(2-q) * (2-2*q)^(1-q)


lamseq <- c(lambda_max_standardized_25, rep(NA, nlam - 1))
z <- c(30, rep(NA, nlam - 1)) # a list to store number of zeros. Changed the name to z to avoid confusion

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
  z[i] <- sum(betacurr == 0)  # Record number of zero coefficients
  cat("z=", z[i], "\n")  # Print for debug
  
  if (z[i] == 0) {
    break  # Stop if no more zero coefficients (fully dense)
  }
}



nz = 30-z 


model = lm(nz ~ log(lamseq))
coefs = round(model$coefficients, 3)

# Save result as a list
result <- list(
  n = n,
  p = p,
  beta.type = beta.type,
  rho = rho,
  snr = snr,
  label = label,
  q = q,
  nz = nz,
  lamseq = lamseq,
  intercept = coefs[1],
  slope = coefs[2],
  unique_nz = length(na.omit(unique(nz)))
)

# Convert q to string with decimal point replaced
q_str <- gsub("\\.", "p", as.character(q))

if (!is.null(rho) && !is.null(snr)) {
  # rho and snr are given separately — use them directly
  rho_str <- gsub("\\.", "p", as.character(rho))
  snr_str <- gsub("\\.", "p", as.character(snr))
  
  filename <- paste0(
    "sim.n", n,
    ".p", p,
    ".beta", beta.type,
    ".rho", rho_str,
    ".snr", snr_str,
    ".q", q_str
  )
  
} else if (!is.null(sim_text)) {
  # Only sim_text is provided — use it as a label
  safe_sim <- sim_text |>
    gsub("\\.", "p", x = _) |>        # decimals → p
    gsub("[^a-zA-Z0-9]", "", x = _)   # remove all non-alphanum
  
  filename <- paste0(
    "sim.n", n,
    ".p", p,
    ".beta", beta.type,
    ".", safe_sim,
    ".q", q_str
  )
}
# Save result
save(result, file = paste0("~/Downloads/Research_Spring_2025/", filename, ".RData"))



########################################################################################
## Load the result
folder <- "~/Downloads/Research_Spring_2025"
# Load and organize all .RData files
all_files <- list.files(path = "~/Downloads/Research_Spring_2025", 
                    pattern = "^sim.*\\.RData$", full.names = TRUE)

# Filter filenames that contain "n100" and "p10"
files <- all_files[grepl("n100", all_files) & grepl("p10", all_files)]

lambda_summary <- lapply(files, function(f) {
  load(f)  # loads `result`
  
  lambda_ratio <- max(result$lamseq, na.rm = TRUE) / min(result$lamseq, na.rm = TRUE)
  
  setting_label <- if (!is.null(result$label)) {
    result$label
  } else {
    paste0("rho = ", result$rho, ", snr = ", result$snr)
  }
  
  data.frame(
    setting = setting_label,
    q = result$q,
    lambda_ratio = lambda_ratio
  )
})

lambda_df <- do.call(rbind, lambda_summary)
# Plot three graphs: one for each setting
library(ggplot2)

# free scales
ggplot(lambda_df, aes(x = q, y = lambda_ratio)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ setting, scales = "free", nrow = 1) +
  labs(
    title = bquote("Lambda Ratio (" ~ lambda[max] / lambda[min] ~ 
                     ") over q, for n = " * .(n) * ", p = " * .(p)),
    x = "q",
    y = expression(lambda[max] / lambda[min])
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = paste0(folder, "/lambda_ratio_plot_n", n, "_p", p, "_free.png"),
  width = 10,
  height = 5
)

# fixed scales
ggplot(lambda_df, aes(x = q, y = lambda_ratio)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ setting, scales = "fixed", nrow = 1) +
  labs(
    title = bquote("Lambda Ratio (" ~ lambda[max] / lambda[min] ~ 
                     ") over q, for n = " * .(n) * ", p = " * .(p)),
    x = "q",
    y = expression(lambda[max] / lambda[min])
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = paste0(folder, "/lambda_ratio_plot_n", n, "_p", p, "_fixed.png"),
  width = 10,
  height = 5
)

# log of lambda max/lambda min
ggplot(lambda_df, aes(x = q, y = lambda_ratio)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ setting, scales = "fixed") +
  scale_y_log10() +
  labs(
    title = bquote("Log Lambda Ratio (" ~ lambda[max] / lambda[min] ~ 
                     ") over q, for n = " * .(n) * ", p = " * .(p)),
    x = "q",
    y = expression(log[10](lambda[max] / lambda[min]))
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = paste0(folder, "/lambda_ratio_plot_n", n, "_p", p, "_log_fixed.png"),
  width = 10,
  height = 5
)


###############################################################################
## Load the result
folder <- "~/Downloads/Research_Spring_2025"
# Load and organize all .RData files
all_files <- list.files(path = "~/Downloads/Research_Spring_2025", 
                        pattern = "^sim.*\\.RData$", full.names = TRUE)
# Filter filenames that contain "n70", "p30", and "snr0p7"
files <- all_files[grepl("n70", all_files) & grepl("p30", all_files) & grepl("snr0p7", all_files)]

library(ggplot2)
library(gridExtra)

# List to store plots
plot_list <- list()

# Loop through up to 4 matching files
for (i in seq_len(min(4, length(files)))) {
  load(files[i])  # loads 'result' list
  
  df <- data.frame(
    log_lambda = log(result$lamseq),
    nz = result$nz
  )
  
  # Extract intercept and slope
  intercept <- result$intercept
  slope <- result$slope
  
  # Label
  label <- paste0("n = ", result$n,
                  ", p = ", result$p,
                  ", q = ", result$q,
                  "\nintercept = ", round(intercept, 3),
                  ", slope = ", round(slope, 3),
                  "\nunique nz = ", result$unique_nz)
  
  # Plot
  p <- ggplot(df, aes(x = log_lambda, y = nz)) +
    geom_point() +
    geom_abline(intercept = intercept, slope = slope, color = "blue") +
    geom_hline(yintercept = 0:30, linetype = "dotted", color = "gray") +
    labs(title = label, x = "log(lambda)", y = "Number of nonzeros") +
    theme_minimal(base_size = 14)
  
  plot_list[[i]] <- p
}

# Show in 2x2 layout
p_grid = grid.arrange(grobs = plot_list, ncol = 2)
ggsave("~/Downloads/Research_Spring_2025/simulations_April_29_sim1.png", p_grid, width = 10, height = 8)


# print only one plot
file <- all_files[17]

load(file)  # loads 'result'

df <- data.frame(
  log_lambda = log(result$lamseq),
  nz = result$nz
)

intercept <- result$intercept
slope <- result$slope

label <- paste0("n = ", result$n,
                ", p = ", result$p,
                ", q = ", result$q,
                "\nintercept = ", round(intercept, 3),
                ", slope = ", round(slope, 3),
                "\nunique nz = ", result$unique_nz)

p <- ggplot(df, aes(x = log_lambda, y = nz)) +
  geom_point() +
  geom_abline(intercept = intercept, slope = slope, color = "blue") +
  geom_hline(yintercept = 0:30, linetype = "dotted", color = "gray") +
  labs(title = label, x = "log(lambda)", y = "Number of nonzeros") +
  theme_minimal(base_size = 14)

print(p)
ggsave(p, filename = "~/Downloads/Research_Spring_2025/simulations_April_29_sim1_rho0p35_snr0p7_s5.png", width = 6, height = 5)

# par(mfrow = c(2, 2))
# plot(log(lamseq), nz)
# # lm$coefficients
# # plot the fitted line
# abline(model, col = "blue")
# abline(h = 0:30, lty = 3)
# # plot(log(lamseq))
# # nz
# # unique(nz)
# # add figure title
# subtext <- paste0(scale_text, ", ",
#                   sim_text,
#                   ", q = ", q,
#                   "\nintercept = ", coefs[1],
#                   ", slope = ", coefs[2],
#                   "\nunique nz = ", length(na.omit(unique(nz))))
# 
# title(main = subtext)
# #length(na.omit(unique(nz)))

