p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
Y = (2 * rbinom(n, 1, 0.5) - 1) * (1 + 9 * (X[,1] > 0))

D = data.frame(X=X, Y=Y)

write.table(D, file="smoke_test/quantile.dat", row.names=FALSE)
