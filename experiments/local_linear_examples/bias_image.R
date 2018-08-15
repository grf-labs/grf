set.seed(1234)

rm(list = ls())

setwd("~/git/grf/experiments/local_linear_examples")

library(grf)
library(glmnet)

mu = function(x){ log(1 + exp(6 * x)) }

# set parameters
n = 600; p = 20; sigma = sqrt(20)
X.test =  matrix(runif(n*p, -1, 1), nrow = n)

# show a grid of predictions along the first coordinate 
ticks = seq(-1, 1, length=n) 
X.test[,1] = ticks
truth = mu(ticks)

# train forest 
forest = local_linear_forest(X, Y, tune.parameters = TRUE)

# lasso to select local linear correction variables 
lasso.mod = cv.glmnet(X, Y, alpha=1)
coefs = coef(lasso.mod)
selected = which(coefs != 0)
# remove intercept and adjust indexes correspondingly 
selected = selected[2:length(selected)] - 1 

# make predictions
preds.llf = predict(forest, X.test, linear.correction.variables = selected, tune.lambda = TRUE)$predictions
preds.rf = predict(forest, X.test)$predictions
df = data.frame(cbind(ticks, truth, preds.llf, preds.rf))

pdf('rf_bias.pdf')
plot(df$ticks, df$truth, 'l', col = "black", main = "Random Forest Predictions",
     xlab = "x", ylab = "y", ylim = range(c(truth, df$preds.rf)))
points(df$ticks, df$preds.rf, col = "red")
dev.off()

pdf('llf_bias.pdf')
plot(df$ticks, df$truth, 'l', col = "black", main = "Local Linear Forest Predictions",
     xlab = "x", ylab = "y", ylim = range(c(truth, df$preds.llf)))
points(df$ticks, df$preds.llf, col = "red")
dev.off()