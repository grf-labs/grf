library(grf)
###
# Here we consider a problem in which we want to impute outcomes that are missing at random.
# We observe (W_i Y_i, W_i) where Y is real and W is binary, and we will use a regression 
# forest to estimate mu(x) = E[Y_i | X_i=x, W_i = 1] using the units with W_i=1.
# Doing this without weights gives mu minimizing an empirical version of E[(Y_i - mu(X_i))^2 | W_i=1],
# squared error over the population of units with nonmissing outcomes,
# while a more natural measure of imputation quality is E[(Y_i - mu(X_i))^2 | W_i=0],
# squared error over the population of units for which we will be imputing outcomes.
# With the inverse propensity weights w(x)=(1-e(x))/e(x), we instead minimize 
# unbiased estimator of this more natural measure, as
# E[w(x) (Y_i - mu(X_i))^2 | W_i=1] = E[ (Y_i - mu(X_i))^2 | W_i=0].
#
# We simulate data from this missing data model, then compare error on the imputed observations
# and check coverage.
###

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
