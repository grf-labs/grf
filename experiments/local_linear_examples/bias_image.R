library(grf)
library(glmnet)

mu = function(x){ log(1 + exp(6 * x)) }

# set parameters
n = 600; p = 20; sigma = sqrt(20)
X = matrix(runif(n*p, -1, 1), nrow = n)
Y = mu(X[,1]) + sigma*rnorm(n)

# random forest predictions
forest = regression_forest(X, Y, tune.parameters = "all")
preds.rf = predict(forest)$predictions

# local linear forest predictions
ll.forest = local_linear_forest(X, Y)

# lasso to select local linear correction variables
lasso.mod = cv.glmnet(X, Y, alpha=1)
selected = as.numeric(predict(lasso.mod, type = "nonzero"))
preds.llf = predict(ll.forest, linear.correction.variables = selected, tune.lambda = TRUE)$predictions

ticks = seq(-1, 1, length = 2000)

pdf('rf_bias.pdf')
plot(X[,1], preds.rf, col = "red", main = "Random Forest Predictions",
     xlab = "x", ylab = "y", ylim = range(c(mu(ticks), preds.rf)))
points(ticks, mu(ticks), 'l', col = "black")
dev.off()

pdf('llf_bias.pdf')
plot(X[,1], preds.llf, col = "red", main = "Local Linear Forest Predictions",
     xlab = "x", ylab = "y", ylim = range(c(mu(ticks), preds.llf)))
points(ticks, mu(ticks), 'l', col = "black")
dev.off()
