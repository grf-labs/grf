library(grf)

save_and_load = function(forest) {
  save('forest', file='forest.Rdata')
  load('forest.Rdata')
  forest
}

n <- 200
p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- X[,1] + 0.1 * rnorm(n)
X.new <- matrix(rnorm(n * p), n, p)

## Time making a forest with and without serialization

system.time(replicate(50, regression_forest(X, Y, serialize=TRUE)))
system.time(replicate(50, regression_forest(X, Y, serialize=FALSE)))

## save and load forest so that it is serialized w/no pointer 
forest.a = regression_forest(X, Y, serialize=TRUE)
forest.b = save_and_load(forest.a)

## Time predictions with and without deserialization
# because forest is deserialized on the first call,
# this should be slow, fast, equally fast
system.time(predict(forest, X.new))
system.time(predict(forest, X.new))
system.time(predict(nonserialized.forest, X.new))


p = 10
n = 500
i = 5
X = matrix(2 * runif(n * p) - 1, n, p)
Y = rnorm(n)  + 100 * (X[,i] > 0)
X.new = matrix(2 * runif(n * p) - 1, n, p)

forest.a = quantile_forest(X, Y)
forest.b = save_and_load(forest.a)
forest.c = save_and_load(forest.a)
predict(forest.a, X[1:4,])
predict(forest.b, X[1:4,])
predict(forest.c, X[1:4,])

all(predict(forest.b,X.new) == predict(forest.c,X.new))
