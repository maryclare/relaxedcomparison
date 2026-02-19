library(powopt)  # Load the powopt package, which includes powCD()

# Set sample size n and number of predictors p
n <- 70
p <- 30

# Generate n x p matrix with small noise, then standardize
X <- matrix(rnorm(n*p, 0, 0.2), n, p)
X <- scale(X)  # Column-wise standardization
u <- svd(X)$u  # Take left singular vectors from SVD
X <- u         # Replace X with orthonormal matrix U (n x n)

# Generate response vector Y ~ N(0,1)
Y <- rnorm(n, 0, 1)

q = 0.25  # Set Lq penalty index

# Compute regularization weights (ω_j), based on inner products with response
r <- abs(t(X) %*% Y)  # Absolute correlation between predictors and response
a <- (2 - q) * q^(1 / (q - 2)) * (2 * (1 - q))^((q - 1) / (2 - q))  # Scaling constant in theory
omega <- r / a  # Compute ω_j

# Compute λ_j values for each variable based on theoretical formula (sorted decreasingly)
lambda_val <- sort(omega^(2 - q) / q, decreasing = TRUE)

# Run powCD at the largest lambda to get the most sparse solution
powCD(X = X, y = Y, 
      sigma.sq = 1, lambda = lambda_val[1], q = q, 
      rand.restart = 0, start = rep(0, p))  # Start from all zeros

# Construct a dense lambda sequence by interpolating between adjacent lambda_vals
lambda_seq <- c(lambda_val[1])
for(i in 1:(length(lambda_val) - 1)){
  lambda_seq <- c(lambda_seq, seq(lambda_val[i], lambda_val[i+1], length.out = 30))
}

# Initialize solution matrix: each column will store solution for a different lambda
solution <- matrix(NA, nrow = p, ncol = length(lambda_seq))

# Solve for the first lambda (starting from zero)
solution[,1] <- powCD(X = X, y = Y, 
                      sigma.sq = 1, lambda = lambda_seq[1], q = q, 
                      rand.restart = 0, start = rep(0, p))

# Solve sequentially across lambda values, warm-starting from previous solution
for(i in 2:length(lambda_seq)){
  solution[,i] <- powCD(X = X, y = Y, 
                        sigma.sq = 1, lambda = lambda_seq[i], q = q, 
                        rand.restart = 0, start = solution[,i - 1])
}

# Show how many distinct sparsity levels (i.e., number of nonzero coefficients) appear
unique(colSums(solution != 0))
plot(lambda_seq, colSums(solution != 0), type = "l", xlab = "Lambda", ylab = "Number of non-zero coefficients")
