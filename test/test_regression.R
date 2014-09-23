
library(randomForest)
library(ranger)

num.replications <- 100
n <- 500

ntree <- 500
mtry <- 3
nodesize <- 10

p <- 6
beta <- c(0,0,0,0,1,2)

dat <- lapply(1:num.replications, function(x) {
  x <- matrix(rbinom(p*n, 1, 0.5), nrow = n)
  y <- x %*% beta + rnorm(n, 0, 0.5)
  dat <- data.frame(y,x)
})

## Test predictions error 
result <- sapply(dat, function(x) {
  rf <- randomForest(y~., x, ntree = ntree, mtry = mtry, nodesize = nodesize)
  rg <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize)
  
  c(rf$mse[ntree], rg$prediction.error)
})

plot(result[1, ], result[2, ], pch = 20, xlab = "RF", ylab = "Ranger")
abline(0,1)