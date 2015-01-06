library(ranger)
library(randomForest)
library(epade)

## ======================================================================
## Simulation parameters
## ======================================================================
n <- 500
num_replicates <- 100
#prop_case <- 0.5
##prevalence <- 0.1
OR <- c(1, 1, 1, 1, 2, 4, 0.25, 0.5)
freq <- 0.5

## ======================================================================
## RF parameters
## ======================================================================
num_trees <- 500
mtry <- 4
nodesize <- 5

## ======================================================================
## Simulate data function
## ======================================================================
simulate_data <- function(n, OR, freq) {
  beta <- log(OR)
  x <- replicate(length(OR), rbinom(n, 1, freq))
  z <- x %*% beta
  prob_y <- 1/(1+exp(-z))
  y <- factor(rbinom(n, 1, prob_y))
  dat <- data.frame(y, x)
  return(dat)
}

## ======================================================================
## Prediction error
## ======================================================================
prediction_error <- replicate(num_replicates, {
  dat <- simulate_data(n, OR, freq)
  
  rf <- randomForest(y ~ ., dat, ntree = num_trees, mtry = mtry, nodesize = nodesize)
  rg <- ranger(y ~ ., dat, num.trees = num_trees, mtry = mtry, min.node.size = nodesize)
  
  c(rf$err.rate[num_trees, 1], rg$prediction.error)
})

## Scatter
plot(prediction_error[1, ], prediction_error[2, ], pch = 19)
abline(0, 1)

## Bland Altman
bland.altman.ade(prediction_error[1, ], prediction_error[2, ], fitline = 2)

## ======================================================================
## Unscaled VIMP
## ======================================================================
vimp <- replicate(num_replicates, {
  dat <- simulate_data(n, OR, freq)
  
  rf <- randomForest(y ~ ., dat, ntree = num_trees, mtry = mtry, nodesize = nodesize,
                     importance = TRUE)
  rg <- ranger(y ~ ., dat, num.trees = num_trees, mtry = mtry, min.node.size = nodesize, 
               importance = "permutation", scale.permutation.importance = FALSE)
  
  cbind(importance(rf, type = 1, scale = FALSE), importance(rg))
})

## Scatter
plot(c(vimp[ , 1, ]), c(vimp[ , 2, ]), pch = 19)
abline(0, 1)

## Bland Altman
bland.altman.ade(c(vimp[ , 1, ]), c(vimp[ , 2, ]), fitline = 2)

## ======================================================================
## Scaled VIMP
## ======================================================================
vimp <- replicate(num_replicates, {
  dat <- simulate_data(n, OR, freq)
  
  rf <- randomForest(y ~ ., dat, ntree = num_trees, mtry = mtry, nodesize = nodesize,
                     importance = TRUE)
  rg <- ranger(y ~ ., dat, num.trees = num_trees, mtry = mtry, min.node.size = nodesize, 
               importance = "permutation", scale.permutation.importance = TRUE)
  
  cbind(importance(rf, type = 1, scale = TRUE), importance(rg))
})

## Scatter
plot(c(vimp[ , 1, ]), c(vimp[ , 2, ]), pch = 19)
abline(0, 1)

## Bland Altman
bland.altman.ade(c(vimp[ , 1, ]), c(vimp[ , 2, ]), fitline = 2)
