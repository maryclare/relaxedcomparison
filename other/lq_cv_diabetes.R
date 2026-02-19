
library(powopt)
library(lars)

mypowCD <- function(ytr, Xtr, yme = 0, 
                    ys = 1, 
                    Xme = rep(0, ncol(Xtr)), 
                    Xs = rep(1, ncol(Xtr)), yte, Xte, nlam = 50,
                    qs = seq(0.1, 1, by = 0.1), epsmin = 10^(-7), 
                    tol = 10^(-14)) {
  
  n <- nrow(Xtr)
  p <- ncol(Xtr)
  ystd <- (ytr - yme)/ys
  Xstd <- sweep(sweep(Xtr, 2, Xme, "-"), 2, Xs, "/")
  Xstdte <- sweep(sweep(Xte, 2, Xme, "-"), 2, Xs, "/")
  
  bm <- crossprod(Xstd, ystd)/diag(crossprod(Xstd))
  bsp <- array(NA, c(length(mylams), length(qs)))
  temse <- array(NA, c(length(mylams), length(qs)))
  
  for (q in qs) {
  lammax <- max((2*(1 - q))^(1 - q)*(2-q)^(q-2)*abs(crossprod(Xstd, ystd))^(2-q)*diag(crossprod(Xstd))^(q - 1))
  lammin <- epsmin/sum(abs(bm)^q)
  
  mylams <- exp(seq(log(lammax), 
                    log(lammin),
                    length.out = 49))
  mylams <- c((mylams[1]/mylams[2])*mylams[1], mylams)
  
  start <- rep(0, ncol(Xstd))
  for (k in 1:length(mylams)) {
    myest <- powCD(y = ystd, X = Xstd, q = q, sigma.sq = n, 
                   lambda = mylams[k]/n, tol = tol,
                   start = start)
    start <- myest
    bsp[k, which(q == qs)] <- sum(myest == 0)
    temse[k, which(q == qs)] <- mean((yte - (ys*(Xstdte%*%myest) + yme))^2)
  }
  }
  return(list("bsp" = bsp, "temse" = temse))
}

data(diabetes)

# Either parametrization produces similar results
y <- (diabetes$y - mean(diabetes$y))/sd(diabetes$y)  
X <- diabetes$x2 

n <- nrow(X)
p <- ncol(X)
x <- X%*%diag(rep(sqrt(n), p))
X <- x

ngroups <- 10
bsp <- cve <- array(NA, dim = c(50, 10, ngroups))
cv.groups <- rep(1:ngroups, each = ceiling(length(y)/ngroups))[1:length(y)]
table(cv.groups)

for (g in 1:ngroups) {
  cat("g=", g, "\n")
ytr <- y[cv.groups != g]
Xtr <- X[cv.groups != g, ]
yte <- y[cv.groups == g]
Xte <- X[cv.groups == g, ]

yme = mean(ytr) 
ys = sqrt(mean((ytr - mean(ytr))^2)) + 1
Xme = apply(Xtr, 2, mean)
Xs = apply(Xtr, 2, function(x) {sqrt(mean((x - mean(x))^2))}) + 1

gf <- mypowCD(ytr = ytr, Xtr = Xtr, yte = yte, Xte = Xte,
              yme = yme, ys = ys, Xme = Xme, Xs = Xs, tol = 10^(-5))
cve[, , g] <- gf$temse
bsp[, , g] <- gf$bsp
}

ecv <- apply(cve, c(1, 2), mean)
which(ecv == min(ecv), arr.ind = TRUE)
bsp[9, 5, ]
which(ecv[, 10] == min(ecv[, 10]), arr.ind = TRUE)
bsp[8, 10, ]

boxplot(t(cve[, 8, ]))
boxplot(t(cve[, 9, ]))
boxplot(t(cve[, 10, ]))

boxplot(t(bsp[, 1, ]))
boxplot(t(bsp[, 5, ]))
boxplot(t(bsp[, 10, ]))


