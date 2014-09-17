
library(randomForest)
library(ranger)

num.replications <- 100
n <- 500

ntree <- 500
mtry <- 3
nodesize <- 3

p <- 6
beta <- c(0,0,0,0,1,2)

result <- replicate(num.replications, {
  x <- matrix(rbinom(p*n, 1, 0.5), nrow = n)
  y <- x %*% beta + rnorm(n, 0, 0.5)
  dat <- data.frame(y,x)
  
  rf <- randomForest(y~., dat, ntree = ntree, mtry = mtry, nodesize = nodesize)
  rg <- ranger(y~., dat, num.trees = ntree, mtry = mtry, min.node.size = nodesize)
  
  c(rf$mse[n], rg$prediction.error)
})

plot(result[1, ], result[2, ], pch = 20, xlab = "RF", ylab = "Ranger")
abline(0,1)
