# examples of using the functions powCD_wrapper, coef.powCD, and predict.powCD

library(devtools)
# install_github("maryclare/powopt")
library(powopt)
# install_github("duanning1/powwrap")
library(powwrap)

x = matrix(rnorm(100), ncol = 10)
y = rnorm(10)
out = powCD_wrapper(x, y, nlambda = 10, nq = 5)
coefficients = coef.powCD(out)
newx = matrix(rnorm(100), ncol = 10)
predictions = predict.powCD(out, newx)
