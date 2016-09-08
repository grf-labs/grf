setwd("~/code/personal/ranger/ranger-r-package/ranger/")
devtools::load_all()
source("R/ranger.R")
source("R/predict.R")

p = 10
n = 100
X = matrix(round(2 * runif(n * p) - 1, 2), n, p)
W = rbinom(n, 1, 1/(1 + exp(X[,3])))
Y = round(100 * (2*W - 1) * (X[,1] > 0) + X[,2] + rnorm(n), 2)
D = data.frame(X=X, Y=Y, W=W)

#rf.causal <- ranger(Y ~ ., D, causal = TRUE, status.variable.name = "W", num.trees = 100, write.forest = TRUE)
rf.instrumental <- ranger(Y ~ ., D, instrumental = TRUE, status.variable.name = "W", instrument.variable.name = "W", 
                          num.trees = 100, write.forest = TRUE)

#causal.predictions <- predict(rf.causal, D)
instrumental.predictions <- predict(rf.instrumental, D)

#plot(causal.predictions$predictions, instrumental.predictions$predictions)
