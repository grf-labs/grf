rm(list = ls())
set.seed(42)

library(grf)
library(maq) # For Qini curves.
library(ggplot2)

# Read in data and specify outcome Y, treatment W, and (numeric) matrix of covariates X.
# data = read.csv("https://raw.githubusercontent.com/grf-labs/grf/master/experiments/IJMPR/synthetic_data.csv")
data = read.csv("synthetic_data.csv") # TODO change with above line later.

Y = data$outcome
W = data$treatment
X = data[, -c(1, 2)]


# This script assumes the covariates X have named columns.
# If not provided, we make up some default names.
if (is.null(colnames(X))) colnames(X) = make.names(1:ncol(X))


# *** Estimating an average treatment effect (ATE) ***

# A simple difference in means estimate ignoring non-random assignment.
summary(lm(Y ~ W))

# A doubly robust ATE estimate (forest-based AIPW).
cf.full = causal_forest(X, Y, W)
average_treatment_effect(cf.full)

# A histogram of the estimated propensity scores.
# Overlap requires that these don't get too close to either 0 or 1.
hist(cf.full$W.hat, xlab = "Estimated propensity scores", main = "")


# *** Estimating CATEs ***

# Split data into a train and test sample.
train = sample(nrow(X), 0.6 * nrow(X))
test = -train

# Fit a CATE function on training data.
cate.forest = causal_forest(X[train, ], Y[train], W[train])

# Predict CATEs on test set.
X.test = X[test, ]
tau.hat.test = predict(cate.forest, X.test)$predictions

# A histogram of CATE estimates.
hist(tau.hat.test, xlab = "Estimated CATEs", main = "")

# On their own, the CATE point estimates are noisy.
# What we often care about is whether they capture meaningful heterogeneity.

# Here we construct groups according to which quartile of the predicted CATEs the unit belongs.
# Then, we calculate ATEs in each of these groups and see if they differ.
num.groups = 4 # 4 for quartiles, 5 for quintiles, etc.
quartile = cut(tau.hat.test,
               quantile(tau.hat.test, seq(0, 1, by = 1 / num.groups)),
               labels = 1:num.groups,
               include.lowest = TRUE)
# Create a list of test set samples by CATE quartile.
samples.by.quartile = split(seq_along(quartile), quartile)

# Look at ATEs in each of these quartiles. To calculate these we fit a separate evaluation forest.
eval.forest = causal_forest(X.test, Y[test], W[test])

# Calculate doubly robust ATEs for each group.
ate.by.quartile = lapply(samples.by.quartile, function(samples) {
  average_treatment_effect(eval.forest, subset = samples)
})

# Plot group ATEs along with 95% confidence bars.
df.plot.ate = data.frame(
  matrix(unlist(ate.by.quartile), num.groups, byrow = TRUE, dimnames = list(NULL, c("estimate","std.err"))),
  group = 1:num.groups
)

ggplot(df.plot.ate, aes(x = group, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.err, ymax = estimate + 1.96 * std.err, width = 0.2)) +
  geom_hline(yintercept = 0, lty = 3) +
  xlab("Estimated CATE quantile") +
  ylab("Average treatment effect")


# *** Evaluate heterogeneity via TOC/AUTOC ***

# In the previous section we split the data into quartiles.
# The TOC is essentially a continuous version of this exercise.

# Use eval.forest to form a doubly robust estimate of TOC/AUTOC.
rate.cate = rank_average_treatment_effect(
  eval.forest,
  tau.hat.test,
  q = seq(0.05, 1, length.out = 100)
)
# Plot the TOC.
plot(rate.cate)

# An estimate and standard error of AUTOC.
print(rate.cate)

# Get a 2-sided p-value Pr(>|t|) for RATE = 0 using a t-value.
2 * pnorm(-abs(rate.cate$estimate / rate.cate$std.err))

# [For alternatives to perform these tests without performing train/test splits,
# including relying on one-sided tests, see the following vignette for more details
# https://grf-labs.github.io/grf/articles/rate_cv.html]


# *** Evaluate CATE models via the AUTOC ***

# Causal forest is a two-step algorithm that first accounts for confounding and baseline effects
# via the propensity score e(x) and a conditional mean model m(x), then in the second step estimates
# treatment effect heterogeneity. In some settings we may want to try using different covariates (or
# possibly models) for e(x), m(x), and CATE predictions.

# Estimate m(x) = E[Y | X = x] using a regression forest.
Y.forest = regression_forest(X[train, ], Y[train], num.trees = 500)
Y.hat = predict(Y.forest)$predictions

# Estimate e(x) = E[W | X = x] using a regression forest.
W.forest = regression_forest(X[train, ], W[train], num.trees = 500)
W.hat = predict(W.forest)$predictions

# Select the covariates X with, for example, m(x) variable importance in the top 25%.
varimp.Y = variable_importance(Y.forest)
selected.vars = which(varimp.Y >= quantile(varimp.Y, 0.75))
print(colnames(X)[selected.vars])

if (length(selected.vars) <= 1) stop("You should really try and use more than just one predictor variable with forests.")

# Try and fit a CATE model using this smaller set of potential heterogeneity predictors.
X.subset = X[, selected.vars]
cate.forest.restricted = causal_forest(X.subset[train, ], Y[train], W[train],
                                       Y.hat = Y.hat, W.hat = W.hat)
# Predict CATEs on test set.
tau.hat.test.restricted = predict(cate.forest.restricted, X.test[, selected.vars])$predictions

# Compare CATE models with AUTOC.
rate.cate.compare = rank_average_treatment_effect(
  eval.forest,
  cbind(tau.hat.test, tau.hat.test.restricted)
)
# Get an estimate of the AUTOCs, as well as difference in AUTOC.
print(rate.cate.compare)

# Get a p-value for the AUTOCs and difference in AUTOCs.
data.frame(
  p.value = 2 * pnorm(-abs(rate.cate.compare$estimate / rate.cate.compare$std.err)),
  target = rate.cate.compare$target
)

# Or equivalently, we could construct a 2-sided confidence interval.
rate.cate.compare$estimate + data.frame(lower = -1.96 * rate.cate.compare$std.err,
                                        upper = 1.96 * rate.cate.compare$std.err,
                                        row.names = rate.cate.compare$target)

# [In this example the restricted CATE model does not do much better.]


# *** Policy evaluation with Qini curves ****

# We can use the `maq` package for this exercise. This package is more general
# and accepts CATE estimates from multiple treatment arms along with costs that
# denominate what we spend by assigning a unit a treatment. In this application
# we can simply treat the number of units we are considering deploying as the cost.

# Form a doubly robust estimate of a CATE-based Qini curve (using eval.forest).
num.units = nrow(X)
qini = maq(tau.hat.test,
           num.units,
           get_scores(eval.forest) * num.units,
           R = 200)

# Form a baseline Qini curve that assigns treatment uniformly.
qini.baseline = maq(tau.hat.test,
                    num.units,
                    get_scores(eval.forest) * num.units,
                    R = 200,
                    target.with.covariates = FALSE)

# Plot the Qini curve along with 95% confidence lines.
plot(qini, ylab = "PTSD cases prevented", xlab = "Units held back from deployment", xlim = c(0, num.units))
plot(qini.baseline, add = TRUE, ci.args = NULL)

# Get estimates from the curve, at for example 500 deployed units.
average_gain(qini, 500)

# Compare the benefit of targeting the 500 units predicted to benefit the most with the baseline.
difference_gain(qini, qini.baseline, 500)


# [The paper shows Qini curves embellished with ggplot. We could have retrieved
# the data underlying the curves and customized our plots further.
# For more details we refer to https://github.com/grf-labs/maq]


# *** Describing the fit CATE function ****

# Our `cate.forest` has given us some estimated function \tau(x).
# Let's have a closer look at how this function stratifies our sample in terms of "covariate" profiles.
# One way to do so is to look at histograms of our covariates by for example low / high CATE predictions.

# First, we'll use a simple heuristic to narrow down the number of predictors to look closer at.
# Here we use the variable importance metric of the fit CATE function to select 4 predictors to look closer at.
varimp.cate = variable_importance(cate.forest)
ranked.variables = order(varimp.cate, decreasing = TRUE)
top.varnames = colnames(X)[ranked.variables[1:4]]
print(top.varnames)

# Select the test set samples predicted to have low/high CATEs.
# [We could also have used the full sample for this exercise.]
low = samples.by.quartile[[1]]
high = samples.by.quartile[[num.groups]]

# Make some long format data frames for ggplot.
df.lo = data.frame(
  covariate.value = unlist(as.vector(X.test[low, top.varnames])),
  covariate.name = rep(top.varnames, each = length(low)),
  cate.estimates = "Low"
)
df.hi = data.frame(
  covariate.value = unlist(as.vector(X.test[high, top.varnames])),
  covariate.name = rep(top.varnames, each = length(high)),
  cate.estimates = "High"
)
df.plot.hist = rbind(df.lo, df.hi)

# Plot overlaid histograms of the selected covariates by low/high classification.
ggplot(df.plot.hist, aes(x = covariate.value, fill = cate.estimates)) +
  geom_histogram(alpha = 0.7, position = "identity") +
  facet_wrap(~ covariate.name, scales = "free", ncol = 2)
  labs(fill = "CATE estimates")


# *** Best linear projections (BLP) ****

# Select some potential effect modifier(s) we are interested in.
blp.vars = c("X1", "X2", "X3")

# Estimate the best linear projection on our variables.
best_linear_projection(cf.full, X[, blp.vars])


# *** Risk vs CATE-based targeting ***

# In our application it was reasonable to hypothesize that soldiers with high
# "risk" of developing PTSD also has a high treatment effect (i.e. low resilience).

# Train a risk model on the training set. First, select units with high combat stress.
train.hi = train[W[train] == 0]

# Use a regression forest to estimate P[develops PTSD | X, high combat stress]
# (We recorded Y = 1 if healthy and so 1 - Y is 1 if the outcome is PTSD)
rf.risk = regression_forest(X[train.hi, ], 1 - Y[train.hi])
risk.hat.test = predict(rf.risk, X.test)$predictions

# Compare risk vs CATE-based targeting using the AUTOC.
rate.risk = rank_average_treatment_effect(
  eval.forest,
  cbind(tau.hat.test, risk.hat.test)
)
plot(rate.risk)
print(rate.risk)

# Construct a 95% confidence interval for the AUTOCs as well as for AUTOC(cate.hat) - AUTOC(risk.hat).
rate.risk$estimate + data.frame(lower = -1.96 * rate.risk$std.err,
                                upper = 1.96 * rate.risk$std.err,
                                row.names = rate.risk$target)
