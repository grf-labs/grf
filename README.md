# generalized random forests <a href='https://grf-labs.github.io/grf/'><img src='https://raw.githubusercontent.com/grf-labs/grf/master/images/logo/grf_logo_wbg_cropped.png' align="right" height="120" /></a>

[![CRANstatus](https://www.r-pkg.org/badges/version/grf)](https://cran.r-project.org/package=grf)
[![](https://cranlogs.r-pkg.org/badges/grand-total/grf)](https://cran.r-project.org/package=grf)
[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.grf?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=2&branchName=master)

A package for forest-based statistical estimation and inference. GRF provides non-parametric methods for heterogeneous treatment effects estimation (optionally using right-censored outcomes, multiple treatment arms or outcomes, or instrumental variables), as well as least-squares regression, quantile regression, and survival regression, all with support for missing covariates.

In addition, GRF supports 'honest' estimation (where one subset of the data is used for choosing splits, and another for populating the leaves of the tree), and confidence intervals for least-squares regression and treatment effect estimation.

Some helpful links for getting started:

- The [R package documentation](https://grf-labs.github.io/grf/) contains usage examples and method reference.
- The [GRF reference](https://grf-labs.github.io/grf/REFERENCE.html) gives a detailed description of the GRF algorithm and includes troubleshooting suggestions.
- For community questions and answers around usage, see [Github issues labelled 'question'](https://github.com/grf-labs/grf/issues?q=label%3Aquestion).

The repository first started as a fork of the [ranger](https://github.com/imbs-hl/ranger) repository -- we owe a great deal of thanks to the ranger authors for their useful and free package.

### Installation

The latest release of the package can be installed through CRAN:

```R
install.packages("grf")
```

`conda` users can install from the [conda-forge](https://anaconda.org/conda-forge/r-grf) channel:

```
conda install -c conda-forge r-grf
```

The current development version can be installed from source using devtools.

```R
devtools::install_github("grf-labs/grf", subdir = "r-package/grf")
```

Note that to install from source, a compiler that implements C++11 or later is required. If installing on Windows, the RTools toolchain is also required.

### Usage Examples

The following script demonstrates how to use GRF for heterogeneous treatment effect estimation. For examples
of how to use other types of forests, please consult the R [documentation](https://grf-labs.github.io/grf/reference/index.html) on the relevant methods.

```R
library(grf)

# Generate data.
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)

# Train a causal forest.
W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
tau.forest <- causal_forest(X, Y, W)

# Estimate treatment effects for the training data using out-of-bag prediction.
tau.hat.oob <- predict(tau.forest)
hist(tau.hat.oob$predictions)

# Estimate treatment effects for the test sample.
tau.hat <- predict(tau.forest, X.test)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 2)

# Estimate the conditional average treatment effect on the full sample (CATE).
average_treatment_effect(tau.forest, target.sample = "all")

# Estimate the conditional average treatment effect on the treated sample (CATT).
average_treatment_effect(tau.forest, target.sample = "treated")

# Add confidence intervals for heterogeneous treatment effects; growing more trees is now recommended.
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 1)

# In some examples, pre-fitting models for Y and W separately may
# be helpful (e.g., if different models use different covariates).
# In some applications, one may even want to get Y.hat and W.hat
# using a completely different method (e.g., boosting).

# Generate new data.
n <- 4000
p <- 20
X <- matrix(rnorm(n * p), n, p)
TAU <- 1 / (1 + exp(-X[, 3]))
W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
Y <- pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU + rnorm(n)

forest.W <- regression_forest(X, W, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions

forest.Y <- regression_forest(X, Y, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

forest.Y.varimp <- variable_importance(forest.Y)

# Note: Forests may have a hard time when trained on very few variables
# (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
# in selection.
selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)

tau.forest <- causal_forest(X[, selected.vars], Y, W,
                            W.hat = W.hat, Y.hat = Y.hat,
                            tune.parameters = "all")

# See if a causal forest succeeded in capturing heterogeneity by plotting
# the TOC and calculating a 95% CI for the AUTOC.
train <- sample(1:n, n / 2)
train.forest <- causal_forest(X[train, ], Y[train], W[train])
eval.forest <- causal_forest(X[-train, ], Y[-train], W[-train])
rate <- rank_average_treatment_effect(eval.forest,
                                      predict(train.forest, X[-train, ])$predictions)
plot(rate)
paste("AUTOC:", round(rate$estimate, 2), "+/", round(1.96 * rate$std.err, 2))
```

### Developing

In addition to providing out-of-the-box forests for quantile regression and causal effect estimation, GRF provides a framework for creating forests tailored to new statistical tasks. If you'd like to develop using GRF, please consult the [algorithm reference](https://grf-labs.github.io/grf/REFERENCE.html) and [development guide](https://grf-labs.github.io/grf/DEVELOPING.html).

### Funding

Development of GRF is supported by the National Institutes of Health, the National Science Foundation, the Sloan Foundation, the Office of Naval Research (Grant N00014-17-1-2131) and Schmidt Futures.

### References

Susan Athey and Stefan Wager.
<b>Estimating Treatment Effects with Causal Forests: An Application.</b>
<i>Observational Studies</i>, 5, 2019.
[<a href="https://doi.org/10.1353/obs.2019.0001">paper</a>,
<a href="https://arxiv.org/abs/1902.07409">arxiv</a>]

Susan Athey, Julie Tibshirani and Stefan Wager.
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2), 2019.
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>,
<a href="https://arxiv.org/abs/1610.01271">arxiv</a>]

Yifan Cui, Michael R. Kosorok, Erik Sverdrup, Stefan Wager, and Ruoqing Zhu.
<b>Estimating Heterogeneous Treatment Effects with Right-Censored Data via Causal Survival Forests.</b>
<i>Journal of the Royal Statistical Society: Series B</i>, 85(2), 2023.
[<a href="https://doi.org/10.1093/jrsssb/qkac001">paper</a>,
<a href="https://arxiv.org/abs/2001.09887">arxiv</a>]

Rina Friedberg, Julie Tibshirani, Susan Athey, and Stefan Wager.
<b>Local Linear Forests.</b> <i>Journal of Computational and Graphical Statistics</i>, 30(2), 2020.
[<a href="https://www.tandfonline.com/doi/abs/10.1080/10618600.2020.1831930">paper</a>,
<a href="https://arxiv.org/abs/1807.11408">arxiv</a>]

Imke Mayer, Erik Sverdrup, Tobias Gauss, Jean-Denis Moyer, Stefan Wager and Julie Josse.
<b>Doubly Robust Treatment Effect Estimation with Missing Attributes.</b>
<i>Annals of Applied Statistics</i>, 14(3), 2020.
[<a href="https://projecteuclid.org/euclid.aoas/1600454872">paper</a>,
<a href="https://arxiv.org/pdf/1910.10624.pdf">arxiv</a>]

Erik Sverdrup, Maria Petukhova, and Stefan Wager
<b>Estimating Treatment Effect Heterogeneity in Psychiatry: A Review and Tutorial with Causal Forests.</b> 2024.
[<a href="TODO">arxiv</a>]

Stefan Wager and Susan Athey.
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523), 2018.
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>,
<a href="https://arxiv.org/abs/1510.04342">arxiv</a>]

Steve Yadlowsky, Scott Fleming, Nigam Shah, Emma Brunskill, and Stefan Wager.
<b>Evaluating Treatment Prioritization Rules via Rank-Weighted Average Treatment Effects.</b> 2021.
[<a href="https://arxiv.org/abs/2111.07966">arxiv</a>]
