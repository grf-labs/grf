library(ranger)
library(grf)
library(microbenchmark)

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

y <- as.matrix(X) %*% matrix(beta, ncol = 1) + beta0
yfac <- as.factor(rbinom(n, size=1, prob = plogis(y)))

## Data
dat_reg <- data.frame(y, x1, x2, x3, x4, x5)
dat_class <- data.frame(yfac, x1, x2, x3, x4, x5)

ranger.reg=ranger(y ~ ., data = dat_reg, num.trees = 5000)
grf.reg=regression.forest(X, y, num.trees = 5000)

## Run benchmark
bm.result = microbenchmark(
  ranger.reg=predict(ranger.reg, dat_reg),
  grf.reg=predict(grf.reg, X),
  times = 3, unit = "s")
print(bm.result)