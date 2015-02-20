library(survival)
library(ranger)

## Test classification
ranger("Species ~ Petal.Length + Sepal.Length", data = iris, importance = "permutation")
ranger("Species ~ .", data = iris)

## Test regression
ranger("Sepal.Width ~ Petal.Length + Sepal.Length + Species", data = iris, 
       importance = "impurity")
ranger("Sepal.Width ~ .", data = iris)

## Test survival
data(veteran, package = "randomSurvivalForest")
xx <- ranger("Surv(time, status) ~ karno + age", data = veteran)
plot(xx$unique.death.times, xx$survival[1,])
yy <- ranger("Surv(time, status) ~ .", data = veteran)

## Test split select weights
weights <- c(0,0.2,0.5,1)
ww <- ranger("Species ~ .", data = iris, split.select.weights = weights, 
             mtry = 3, importance = "impurity", verbose = TRUE, num.trees = 1)

## Test verbose output
temp <- ranger("Sepal.Width ~ .", data = iris, verbose = TRUE)
temp <- ranger("Sepal.Width ~ .", data = iris, verbose = FALSE)

## Test write forest
temp <- ranger("Sepal.Width ~ .", data = iris, verbose = FALSE, write.forest = TRUE)
temp <- ranger("Surv(time, status) ~ .", data = veteran, verbose = FALSE, write.forest = TRUE)

## Test prediction
temp <- ranger("Species ~ .", data = iris, verbose = FALSE, write.forest = TRUE)
pred <- predict(temp$forest, data = iris)
pred.nodep <- predict(temp$forest, data = iris[, 1:4])

temp <- ranger("Sepal.Width ~ .", data = iris, verbose = FALSE, write.forest = TRUE)
pred <- predict(temp$forest, data = iris)
pred.nodep <- predict(temp$forest, data = iris[, c(1, 2:5)])

data(veteran, package = "randomSurvivalForest")
temp <- ranger("Surv(time, status) ~ karno + age", data = veteran, verbose = FALSE, write.forest = TRUE)
pred <- predict(temp$forest, data = veteran)
pred.nodep <- predict(temp$forest, data = veteran[, c(1:2, 5:8)])
pred.nodep <- predict(temp$forest, data = veteran[, c(1:3, 5:8)])
pred.nodep <- predict(temp$forest, data = veteran[, c(1:2, 4:8)])

temp <- ranger("Species ~ .", data = iris, verbose = FALSE, write.forest = TRUE, probability = TRUE)
pred <- predict(temp$forest, data = iris)
pred.nodep <- predict(temp$forest, data = iris[, 1:4])

## Test GWA mode
library(GenABEL)
##convert.snp.ped("../../gwa_data/chr1.allChunks.ped", "../../gwa_data/chr1.allChunks.map", "../../gwa_data/chr1.allChunks.raw")
chr21 <- load.gwaa.data("../../gaw_data/pheno.GenABEL", "../../gaw_data/chr21.allChunks.raw")
phdata(chr21)$hyper <- factor(phdata(chr21)$hyper)
temp <- ranger("hyper ~ .", data = chr21, write.forest = TRUE)
pred <- predict(temp, data = chr21[1:10,])

##chr1 <- load.gwaa.data("../../gwa_data/pheno.GenABEL", "../../gwa_data/chr1.allChunks.raw")
##phdata(chr1)$hyper <- factor(phdata(chr1)$hyper)
##ranger("hyper ~ .", data = chr1)

## Test memory modes
temp <- ranger("Species ~ .", data = iris, verbose = FALSE, write.forest = TRUE, memory = "double")
temp$prediction.error
pred.double <- predict(temp$forest, data = iris)
temp <- ranger("Species ~ .", data = iris, verbose = FALSE, write.forest = TRUE, memory = "float")
temp$prediction.error
pred.float <- predict(temp$forest, data = iris)
temp <- ranger("Species ~ .", data = iris, verbose = FALSE, write.forest = TRUE, memory = "char")
temp$prediction.error
pred.char <- predict(temp$forest, data = iris)

## Test non-formula interface
ranger(Species ~., data = iris)
ranger(data = iris, dependent.variable.name = "Species")

ranger(Sepal.Length ~., data = iris)
ranger(data = iris, dependent.variable.name = "Sepal.Length")

data(veteran, package = "randomSurvivalForest")
ranger(Surv(time, status) ~., data = veteran)
ranger(data = veteran, dependent.variable.name = "time", status.variable.name = "status")


## Error
ranger(data = iris, dependent.variable.name = "foo")

library(ranger)
library(randomForest)
library(randomSurvivalForest)
library(randomForestSRC)
library(survival)

## Compare results: Classification
res.rg <- ranger(Species ~ ., data = iris, num.trees = 10000, importance = "impurity")
res.rf <- randomForest(Species ~., data = iris, ntree = 10000, importance = TRUE)
res.rg
res.rf
sum(res.rg$predictions != res.rf$predicted)
res.rg$variable.importance
res.rf$importance

dat <- airquality[complete.cases(airquality), ]
dat$Month <- factor(dat$Month)
ranger(Month ~., data = dat)
randomForest(Month ~., data = dat)

## Permutation importance
res.rg <- ranger(Species ~ ., data = iris, importance = "permutation")
res.rf <- randomForest(Species ~., data = iris, importance = TRUE)
res.rg$variable.importance
importance(res.rf, type=1, scale = FALSE)

## Compare results: Regression
res.rg <- ranger(Sepal.Width ~ ., data = iris, importance = "impurity")
res.rf <- randomForest(Sepal.Width ~., data = iris, importance = TRUE)
res.rg
res.rf
res.rg$variable.importance
res.rf$importance

dat <- airquality[complete.cases(airquality), ]
ranger(Temp ~., data = dat)
randomForest(Temp ~., data = dat)

## Compare results: Survival
res.rg <- ranger(Surv(time, status) ~ ., data = veteran, num.trees = 10000, importance = "permutation")
res.src <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 10000, importance = "permute")
res.rsf <- rsf(Surv(time, status) ~ ., data = veteran, ntree = 10000, importance = "permute")
res.rg
res.src
res.rsf
res.rg$variable.importance
res.src$importance
res.rsf$importance

## Permutation importance
n <- 50
mult.rg <- replicate(n, ranger(Surv(time, status) ~ ., data = veteran, num.trees = 10000, importance = "permutation", mtry = 3)$variable.importance)
mult.src <- replicate(n, rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 10000, importance = "permute", mtry = 3)$importance)
mult.rsf <- replicate(n, rsf(Surv(time, status) ~ ., data = veteran, ntree = 10000, importance = "permute", mtry = 3)$importance)

boxplot(cbind(t(mult.rg), t(mult.src), t(mult.rsf)))

## Test with many rows
rows <- 10
n <- 100000
dat.long <-  data.frame(replicate(rows, rbinom(n, 1, 0.5)))
colnames(dat.long) <- paste("x", 1:rows, sep = "")
ranger(x1~., dat.long)

##----------------------------------------------------------------------
## Test unordered variables
##----------------------------------------------------------------------
library(ranger)
library(survival)
library(GenABEL)
rmultinomfactor <- function(n, size, prob) {
  dummy <- rmultinom(n, size, prob)
  nofac <- apply(dummy, 2, function(x) which(x == 1))
  fac <- factor(nofac)
  levels(fac) <- sample.int(length(prob))
  fac
}
n <- 1000

## Covariates
x1 <- rmultinomfactor(n, 1, rep(1, 2))
x2 <- rmultinomfactor(n, 1, rep(1, 4))
x3 <- rmultinomfactor(n, 1, rep(1, 6))
x4 <- rmultinomfactor(n, 1, rep(1, 8))
x5 <- rnorm(n)
X <- cbind(x1, x2, x3, x4, x5)
beta <- c(3,1,2,1,4)
beta0 <- -15

## Endpoint
y <- as.matrix(X) %*% matrix(beta, ncol = 1) + beta0
yfac <- as.factor(rbinom(n, size=1, prob = plogis(y)))
time <- y + rnorm(n)
time <- time - min(time)
status <- rbinom(n, 1, 0.8)

## Data
dat_reg <- data.frame(y, x1, x2, x3, x4, x5)
dat_class <- data.frame(yfac, x1, x2, x3, x4, x5)
dat_surv <- data.frame(time, status, x1, x2, x3, x4, x5)

## Classification
ranger(yfac ~ ., data = dat_class, respect.unordered.factors = FALSE)
ranger(yfac ~ ., data = dat_class, respect.unordered.factors = TRUE)
ranger(dependent.variable.name = "yfac", data = dat_class, respect.unordered.factors = FALSE)
ranger(dependent.variable.name = "yfac", data = dat_class, respect.unordered.factors = TRUE)

## Regression
ranger(y ~ ., data = dat_reg, respect.unordered.factors = FALSE)
ranger(y~ ., data = dat_reg, respect.unordered.factors = TRUE)
ranger(dependent.variable.name = "y", data = dat_reg, respect.unordered.factors = FALSE)
ranger(dependent.variable.name = "y", data = dat_reg, respect.unordered.factors = TRUE)

## Probability prediction
ranger(yfac ~ ., data = dat_class, probability = TRUE, respect.unordered.factors = FALSE)
ranger(yfac ~ ., data = dat_class, probability = TRUE, respect.unordered.factors = TRUE)
ranger(dependent.variable.name = "yfac", data = dat_class, probability = TRUE, respect.unordered.factors = FALSE)
ranger(dependent.variable.name = "yfac", data = dat_class, probability = TRUE, respect.unordered.factors = TRUE)

## Survival
ranger(Surv(time, status) ~ ., data = dat_surv, respect.unordered.factors = FALSE)
ranger(Surv(time, status) ~ ., data = dat_surv, respect.unordered.factors = TRUE)
ranger(dependent.variable.name = "time", status.variable.name = "status", data = dat_surv, respect.unordered.factors = FALSE)
ranger(dependent.variable.name = "time", status.variable.name = "status", data = dat_surv, respect.unordered.factors = TRUE)

## GWAS data
chr21 <- load.gwaa.data("~/sandbox/attic/medianImp/test.pheno", "~/sandbox/attic/medianImp/test.raw")
phdata(chr21)$trait <- factor(phdata(chr21)$trait)
phdata(chr21)$sex <- factor(phdata(chr21)$sex)
ranger("trait ~ .", data = chr21, write.forest = TRUE, respect.unordered.factors = FALSE)
ranger("trait ~ .", data = chr21, write.forest = TRUE, respect.unordered.factors = TRUE)

## Test maximum number of factor levels
n <- 30
x <- factor(1:n)
y <- rbinom(n, 1, 0.5)
dat <- data.frame(y, x)
ranger(y ~ x, data = dat, respect.unordered.factors = TRUE)
