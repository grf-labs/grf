p = 20
n = 100

tau = function(x) (2 * as.numeric(x > -1/3))

X = matrix(2 * runif(n * p) - 1, n, p)
I = rbinom(n, 1, 1/(1 + exp(X[,2])))
A = rbinom(n, 1, 1/(1 + exp(X[,3])))
W = I * A
Y = (W - 1/2) * tau(X[,1]) + 3/2 * (2*A - 1) * (X[,1] > 1/3) + 2 * rnorm(n)

D = data.frame(X=X, Y=Y)

M = as.matrix(sapply(D, as.numeric))

forest <- regression_train(M, 20, matrix(NA, 0, 0), character(0), 5, 4, FALSE, 1, 5, TRUE, FALSE, 0.63, numeric(0), 42)

regression_predict(forest, M, matrix(NA, 0, 0), character(0), 3)