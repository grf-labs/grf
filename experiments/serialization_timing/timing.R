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

## save and load forest so forest.b is serialized w/no xptr
forest.a = lapply(1:50, function(i) { regression_forest(X, Y, serialize=TRUE) })
forest.b = lapply(forest.a, save_and_load)

## Time predictions with and without deserialization
# because forest.b is deserialized on the first call,
# this should be fast, slow, fast
system.time(lapply(forest.a, function(forest) { predict(forest, X.new) }))
system.time(lapply(forest.b,  function(forest) { predict(forest, X.new) }))
system.time(lapply(forest.b,  function(forest) { predict(forest, X.new) }))

