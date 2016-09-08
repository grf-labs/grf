p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)

D = data.frame(X=X, Y=Y, W=W)

write.table(D, file="smoke_test/causal.dat", row.names=FALSE, quote=FALSE)
