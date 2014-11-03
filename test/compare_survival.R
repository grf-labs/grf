
library(survival)
library(ranger)
library(randomForestSRC)
library(pec)
library(prodlim)

predictSurvProb.ranger <- function(object, newdata, times, ...) {
  pred <- predict(object, newdata)
  pos <- sapply(times, function(x) {
    max(c(0, which(x >= pred$unique.death.times)))
  })    
  p <- cbind(1, predictions(pred))[, pos+1, drop = FALSE]
  return(p)
}

num_replicates <- 20
ntree <- 500
mtry <- 3
nodesize <- 10
formula <- Surv(time, status) ~.

n <- 500
p <- 6
beta <- c(1,0.5,0.1,0,1,2)

res <- replicate(num_replicates, {
  ## Binary covariates
  x <- matrix(rbinom(p*n, 1, 0.5), nrow = n)
  time <- x %*% beta + rnorm(n, 0, 0.5)
  time <- time - min(time)
  status <- rbinom(n, 1, 0.7)
  dat <- data.frame(time, status, x)
  idx <- sample(n, 2/3*n)
  dat_train <- dat[idx, ]
  dat_test <- dat[-idx, ]
  time <- sort(unique(dat_train$time))
    
  rg = ranger(formula, dat_train, num.trees = ntree, mtry = mtry, min.node.size = nodesize, write.forest = TRUE)
  src = rfsrc(formula, dat_train, ntree = ntree, mtry = mtry, nodesize = nodesize)
  
  fit.pec <- pec(list(Ranger = rg, rfsrc = src), formula = formula, data = dat_test, times = time, cens.model = "marginal")
  crps(fit.pec)[, 1]
})

boxplot(t(res[2:3,]))
