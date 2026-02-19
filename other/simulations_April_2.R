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
rho = 0.1
snr = 0.75
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
lambda_max_standardized <- (max_abs_x_y_dot_standardized / (2-q) / n) ^(2-q) * (2-2*q)^(1-q)


# Extract lambda_max from the glmnet model
# lambda_max_glmnet <- fit.lasso$lambda[1]

# epsilon <- 0.00001

# Calculate lambda_min
# lambda_min <- epsilon * lambda_max_standardized

# Generate the sequence of lambda values on a log scale
#a <- 0.5
#lambda_seq_manual <- seq((lambda_max_standardized)^(a), 
#                            (lambda_min)^(a), length.out = nlam)^(1/a)

# lambda_seq_manual <- exp(seq(log(lambda_max_standardized), 
#                                log(lambda_min), length.out = nlam))

sorted_abs_x_y_dot_standardized <- sort((abs(x_y_dot_standardized) / (2-q) / n)^(2-q) * (2 - 2*q)^(1-q), 
                                        decreasing = TRUE)
refine_lambda_seq <- function(seq, num_points = 30) {
  refined_seq <- c()
  for (i in seq_along(seq)[-length(seq)]) {
    refined_seq <- c(refined_seq, seq(from = seq[i], to = seq[i + 1], length.out = num_points + 2)[-1])
  }
  return(refined_seq)
}
num_lambda <- 50  # Number of points between each pair
lambda_seq_manual <- refine_lambda_seq(sorted_abs_x_y_dot_standardized, num_lambda)

beta.pow = matrix(0,p, length(lambda_seq_manual))
beta.pow[,1] = powCD(X = x_standardized, y = y, 
                        sigma.sq = n, lambda = lambda_seq_manual[1], 
                        q = q, rand.restart = 0, start = rep(0,p))
for (i in 2:length(lambda_seq_manual)) {
  beta.pow[,i] = powCD(X = x_standardized, y = y, 
                          sigma.sq = n, lambda = lambda_seq_manual[i], q = q, 
                          max.iter = 10000, tol = 1e-7,
                          rand.restart = 0, start = beta.pow[,i-1])
}

nzs.pow = colSums(beta.pow != 0)
plot(log(lambda_seq_manual), nzs.pow)
unique(nzs.pow)
length(unique(nzs.pow))

nz = nzs.pow
lamseq = lambda_seq_manual

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
    "sortedlambda.n", n,
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
    "sortedlambda.n", n,
    ".p", p,
    ".beta", beta.type,
    ".", safe_sim,
    ".q", q_str
  )
}
# Save result
save(result, file = paste0("~/Downloads/Research_Spring_2025/", filename, ".RData"))



########################################################################################
## Later, load the result
folder <- "~/Downloads/Research_Spring_2025"
# Load and organize all .RData files
all_files <- list.files(path = "~/Downloads/Research_Spring_2025",
                        pattern = "^sortedlambda.*\\.RData$", full.names = TRUE)

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
  filename = paste0(folder, "/sortedlambda_lambda_ratio_plot_n", n, "_p", p, "_free.png"),
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
  filename = paste0(folder, "/sortedlambda_lambda_ratio_plot_n", n, "_p", p, "_fixed.png"),
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
  filename = paste0(folder, "/sortedlambda_lambda_ratio_plot_n", n, "_p", p, "_log_fixed.png"),
  width = 10,
  height = 5
)


###############################################################################
# Filter filenames that contain "n70", "p30", and "snr0p7"
files <- all_files[grepl("sortedlambda", all_files) &
                     grepl("n70", all_files) & grepl("p30", all_files) & grepl("snr0p8", all_files)]

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
                  #", Orthogonal",
                  ", rho = ", result$rho,
                  ", snr = ", result$snr,
                  "\nq = ", result$q,
                  ", intercept = ", round(intercept, 3),
                  "\nslope = ", round(slope, 3),
                  ", unique nz = ", result$unique_nz)
  
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
ggsave("~/Downloads/Research_Spring_2025/sortedlambda_simulations_April_29_rho0p2_snr0p8.png", p_grid, width = 10, height = 8)



# # par(mfrow = c(2, 2))
# plot(log(lamseq), nz)
# # lm$coefficients
# # plot the fitted line
# abline(model, col = "blue")
# abline(h = 0:30, lty = 3)
# # plot(log(lamseq))
# # nz
# # unique(nz)
# # add figure title
# subtext <- paste0(scale_text,
#                   ", rho = ", rho,
#                   ", snr = ", snr,
#                   ", q = ", q,
#                   "\nintercept = ", coefs[1],
#                   ", slope = ", coefs[2],
#                   "\nunique nz = ", length(na.omit(unique(nz))))
# 
# title(main = subtext)
# #length(na.omit(unique(nz)))
