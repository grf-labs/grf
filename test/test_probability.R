
library(ranger)

num.replications <- 100
n <- 500

ntree <- 500
mtry <- 3
nodesize <- 50

p <- 6
beta <- c(-5,-7,2,4,1,2)

dat <- lapply(1:num.replications, function(x) {
  x <- matrix(rbinom(p*n, 1, 0.5), nrow = n)
  y.linear <- x %*% beta + rnorm(n, 0, 0.5)
  y <- rbinom(n, size=1, prob = plogis(y.linear))
  dat <- data.frame(y,x)
})

## Error
result <- sapply(dat, function(x) {
  
  rg.reg <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize)
  x$y <- as.factor(x$y)
  rg.prb <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, probability = TRUE)
  
  c(rg.reg$prediction.error, rg.prb$prediction.error)
})

plot(result[1, ], result[2, ], pch = 20, xlab = "Regression", ylab = "Probability prediction")
abline(0,1)

## Predictions
x <- dat[[5]]
x[1,1] <- 1
rg.reg <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize)

x$y <- factor(x$y)
rg.prb <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, probability = TRUE)

plot(rg.reg$predictions, rg.prb$predictions[, 2], pch = 20, xlab = "Regression", ylab = "Probability prediction")
abline(0,1)

## Gini Importance
result.vimp <- lapply(dat, function(x) {
  
  rg.reg <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, importance = "impurity")
  x$y <- as.factor(x$y)
  rg.prb <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, importance = "impurity", probability = TRUE)
  
  cbind(importance(rg.reg), importance(rg.prb))
})
vimp <- simplify2array(result.vimp)
vimp.reg <- vimp[, 1, ]
vimp.prb <- vimp[, 2, ]

par(mfrow=c(1,2))
boxplot(t(vimp.reg))
boxplot(t(vimp.prb))

## Permutation Importance
result.vimp <- lapply(dat, function(x) {
  
  rg.reg <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, importance = "permutation")
  x$y <- as.factor(x$y)
  rg.prb <- ranger(y~., x, num.trees = ntree, mtry = mtry, min.node.size = nodesize, importance = "permutation", probability = TRUE)
  
  cbind(importance(rg.reg), importance(rg.prb))
})
vimp <- simplify2array(result.vimp)
vimp.reg <- vimp[, 1, ]
vimp.prb <- vimp[, 2, ]

par(mfrow=c(1,2))
boxplot(t(vimp.reg))
boxplot(t(vimp.prb))


