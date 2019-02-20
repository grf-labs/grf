library(grf)

n <- 2000
p <- 2
X <- matrix(rnorm(n * p), n, p)
mu = X[,1]
Y <- mu + 0.1 * rnorm(n)

e = 1/(1+exp(-1*X[,1]))
w = runif(n) <= e
sample.weights <- (1-e[w])/e[w] 

forest = regression_forest(X[w,], Y[w])
forest.weighted = regression_forest(X[w,], Y[w], sample.weights)

mu.forest = predict(forest, X[!w,], estimate.variance=TRUE)
mu.forest.weighted = predict(forest.weighted, X[!w,], estimate.variance=TRUE)
z.forest = abs(mu.forest$predictions - mu[!w])/ sqrt(mu.forest$variance.estimates)
z.forest.weighted = abs(mu.forest.weighted$predictions - mu[!w])/ sqrt(mu.forest.weighted$variance.estimates)

mean((mu.forest$predictions - mu[!w]))^2
mean((mu.forest.weighted$predictions - mu[!w]))^2

mean(abs(z.forest <= 1.96))
mean(abs(z.forest.weighted <= 1.96))
