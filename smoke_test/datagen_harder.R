p = 10
n = 1000
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
X[1:10,1] = round(seq(-1, 1, length=10), 2)
Y = round(rnorm(n) * (1 + 9 * (X[,1] > 0)), 3)

D = data.frame(X=X, Y=Y)

write.table(D, file="smoke_test/harder.dat", row.names=FALSE)
