
library(survival)
library(ranger)
library(randomForestSRC)
library(microbenchmark)

ntree <- 500
mtry <- 3
nodesize <- 10
formula <- Surv(time, status) ~.

n <- 500
p <- 6
beta <- c(1,0.5,0.1,0,1,2)

## Binary covariates
x <- matrix(rbinom(p*n, 1, 0.5), nrow = n)
time <- x %*% beta + rnorm(n, 0, 0.5)
time <- time - min(time)
status <- rbinom(n, 1, 0.7)
dat <- data.frame(time, status, x)
microbenchmark(
  Ranger = ranger(formula, dat, num.trees = ntree, mtry = mtry, min.node.size = nodesize), 
  SRC = rfsrc(formula, dat, ntree = ntree, mtry = mtry, nodesize = nodesize), 
  times = 1)

## Continuous covariates
x <- matrix(rnorm(p*n, 1, 5), nrow = n)
time <- x %*% beta + rnorm(n, 0, 0.5)
time <- time - min(time)
status <- rbinom(n, 1, 0.7)
dat <- data.frame(time, status, x)
microbenchmark(
  Ranger = ranger(formula, dat, num.trees = ntree, mtry = mtry, min.node.size = nodesize), 
  SRC = rfsrc(formula, dat, ntree = ntree, mtry = mtry, nodesize = nodesize), 
  times = 1)

