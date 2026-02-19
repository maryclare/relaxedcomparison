library(powopt)


n <- 70
p <- 30

X <- matrix(rnorm(n*p, 0, 0.2), n, p)
X <- scale(X)
u <- svd(X)$u
X <- u
Y <- rnorm(n, 0, 1)

q = 0.25

r <- abs(t(X) %*% Y) 
a <- (2-q)*q^(1/(q-2)) * (2*(1-q))^((q-1)/(2-q))
omega <- r / a
lambda_val <- sort(omega^(2-q) / q, decreasing = T)
#  lambda_val <-sort( (abs(colSums(X * Y) / n) / (2-q)) ^(2-q) * (2-2*q)^(1-q), decreasing = T)

powCD(X = X, y = Y, 
      sigma.sq = n, lambda = lambda_val[1], q = q, 
      rand.restart = 0, start = rep(0, p))

lambda_seq <- c(lambda_val[1])
for(i in 1:(length(lambda_val)-1)){
  lambda_seq <- c(lambda_seq, seq(lambda_val[i], lambda_val[i+1], length.out = 30))
}

solution <- matrix(NA, nrow = p, ncol = length(lambda_seq))
solution[,1] <- powCD(X = X, y = Y, 
                      sigma.sq = 1, lambda = lambda_seq[1], q = q, 
                      rand.restart = 0, start = rep(0, p))
for(i in 2:length(lambda_seq)){
  solution[,i] <- powCD(X = X, y = Y, 
                        sigma.sq = 1, lambda = lambda_seq[i], q = q, 
                        rand.restart = 0, start = solution[,i-1])
}
unique(colSums(solution != 0))

