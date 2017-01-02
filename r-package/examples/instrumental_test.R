p = 20
n = 100

X = matrix(2 * runif(n * p) - 1, n, p)
I = rbinom(n, 1, 1/(1 + exp(X[,2])))
A = rbinom(n, 1, 1/(1 + exp(X[,3])))
W = I * A
Y = (W - 1/2) * tau(X[,1]) + 3/2 * (2*A - 1) * (X[,1] > 1/3) + 2 * rnorm(n)

D = data.frame(X=X, Y=Y, W=W, I=I)

M = as.matrix(sapply(D, as.numeric))

forest <- instrumental_train(M, 20, 21, 22, matrix(NA, 0, 0), character(0), 5, 4, FALSE, 2, 5, TRUE, FALSE, 0.63)

instrumental_predict(forest, M, matrix(NA, 0, 0), character(0), 3)