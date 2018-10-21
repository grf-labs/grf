# The GRF Algorithm

The following guide gives an introduction to the generalized random forests algorithm as implemented in the `grf` package. It aims to give a complete description of the training and prediction procedures, as well as the options available for tuning. This guide is intended as an informal and practical reference; for a theoretical treatment of GRF, please consult the 'Generalized Random Forests' paper.

GRF extends the idea of a classic random forest to allow for estimating other statistical quantities besides the expected outcome. Each forest type, for example`quantile_forest`, trains a random forest targeted at a particular problem, like quantile estimation. The most common use of GRF is in estimating treatment effects through the function `causal_forest`.

## Table of Contents
* [General Algorithm](#general-algorithm)
  * [Training](#training)
  * [Prediction](#prediction)
  * [Out-of-bag Prediction](#out-of-bag-prediction)
  * [Training Options](#training-options)
  * [Variance Estimates](#variance-estimates)
* [Causal Forests](#causal-forests)
  * [Orthogonalization](#orthogonalization)
  * [Selecting Balanced Splits](#selecting-balanced-splits)
  * [Average Treatment Effects](#average-treatment-effects)
* [Additional Features](#additional-features)
  * [Parameter Tuning](#parameter-tuning)
  * [Cluster-Robust Estimation](#cluster-robust-estimation)
* [Troubleshooting](#troubleshooting)
* [References](#references)

## General Algorithm

In this section, we describe GRF's overall approach to training and prediction. The descriptions given in this section apply to all the available forest models. Specific details about the `causal_forest` method can be found in the 'Causal Forests' section below.

We begin with a simple example to illustrate the estimation process:

```
# Train a causal forest.
n = 2000; p = 10
X = matrix(rnorm(n*p), n, p)
W = rbinom(n, 1, 0.5)
Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
causal.forest = causal_forest(X, Y, W)

# Estimate causal effects on new test data.
X.test = matrix(0, 100, p)
X.test[,1] = seq(-2, 2, length.out = 100)
predictions = predict(causal.forest, X.test)$predictions

# Estimate causal effects for the training data using out-of-bag prediction.
oob.predictions = predict(causal.forest)$predictions
```

We now explore each of these steps in more detail.

### Training

A random forest is at its core an ensemble model, composed of a group of decision trees. During training, a number of trees are grown on random subsamples of the dataset. Individual trees are trained through the following steps:
- First, a random subsample is drawn by sampling without replacement from the full dataset. A single root node is created containing this random sample.
- The root node is split into child nodes, and child nodes are split recursively to form a tree. The procedure stops when no nodes can be split further. Each node is split using the following algorithm:
    - A random subset of variables are selected as candidates to split on.
    - For each of these variables `x`, we look at all of its possible values `v` and consider splitting it into two children based on this value. The goodness of a split (`x`, `v`) is determined by how much it increases heterogeneity in the quantity of interest.  Certain splits are not considered, because the resulting child nodes would be too small or too different in size.
    - All examples with values for the split variable`x` that are less than or equal to the split value `v` are placed in a new left child node, and all examples with values greater than the `v` are placed in a right child node.
    - If a node has no valid splits, or if splitting will not result in an improved fit, the node is not split further and forms a leaf of the final tree.

The main difference between GRF's approach to growing trees and that of classic random forests is in how the quality of a split is measured. Because the various forest types seek to estimate different statistical quantities like quantiles and treatment effects, splitting must be tailored to the particular task at hand. The approach taken in GRF is to maximize the heterogeneity in the quantity of interest across the child nodes. For example, with causal effect estimation, the goodness of a split relates to how different the treatment effect estimates are in each node. A theoretical motivation for this split criterion can be found in section 2 of the GRF paper.

The quality of a split must be calculated for each possible split variable `x` and value `v`, so it is critical for it to be fast to compute. Optimizing the heterogeneity criterion directly is therefore too expensive; instead, we take the gradient of the objective and optimize a linear approximation to the criterion. This approach allows us to reuse classic tree splitting algorithms that work in terms of cumulative sums. This gradient-based approach also has close ties to the concept of 'influence functions', as discussed in 2.3 of the GRF paper.

### Prediction

Given a test example, the GRF algorithm computes a prediction as follows:
- For each tree, the test example is 'pushed down' to determine what leaf it falls in.
- Given this information, we create a list of neighboring training examples, weighted by how many times the example fell in the same leaf as the test example.
- A prediction is made using this weighted list of neighbors, using the relevant approach for the type of forest. For regression forests, the prediction is equal to the average outcome of the test example's neighbors. In causal prediction, we calculate the treatment effect using the outcomes and treatment status of the neighbor examples.

Those familiar with classic random forests might note that this approach differs from the way forest prediction is usually described. The traditional view is that to predict for a test example, each tree makes a prediction on that example. To make a final prediction, the tree predictions are combined in some way, for example through averaging or through 'majority voting'. It's worth noting that for regression forests, the GRF algorithm described above is identical this 'ensemble' approach, where each tree predicts by averaging the outcomes in each leaf, and predictions are combined through a weighted average.

### Out-of-bag Prediction

If a dataset is provided to the `predict` method, then predictions are made for these new test example. When no dataset is provided, prediction proceeds on the training examples. In particular, for each training example, all the trees that did not use this example during training are identified (the example was 'out-of-bag', or OOB). Then, a prediction for the test example is made using only these trees. These out-of-bag predictions can be useful in understanding the model's goodness-of-fit, and are also used in several of the methods for causal effect estimation methods described later in this guide.

### Training Options

#### `sample.fraction`

The `sample.fraction` parameter is a number in the range (0, 1] that controls the fraction of examples that should be used in growing each tree. By default, `sample.fraction` is set to 0.5. As noted in the section on honest forests, the fractional subsample will be further split into halves when honesty is enabled.

#### `num.trees`

The parameter `num.trees` controls how many trees are grown during training, and defaults to 2000. Generally, obtaining high-quality confidence intervals requires growing more trees than are needed for accurate predictions.

Tree training is parallelized across several threads in an effort to improve performance. By default, all available cores are used, but the number of threads can be set directly through `num.threads`.

#### `honesty`

By default, 'honest' forests are trained. In a classic random forest, a single subsample is used both to choose a tree's splits, and for the leaf node examples used in making predictions. In contrast, honest forests randomly split this subsample in half, and use only the first half when performing splitting. The second half is then used to populate the tree's leaf nodes: each new example is 'pushed down' the tree, and added to the leaf in which it falls. In a sense, the leaf nodes are 'repopulated' after splitting using a fresh set of examples.

The motivation behind honesty is to reduce bias in tree predictions, by using different subsamples for constructing the tree and for making predictions. Honesty is a well-explored idea in the academic literature on random forests, but is not yet common in software implementations. For a more formal overview, please see section 2.4 of Wager and Athey (2018).

It's important to note that honesty may hurt performance when working with very small datasets. In this set-up, the subsample used to determine tree splits is already small, and honesty further cuts this subsample in half, so there may no longer be enough information to choose high-quality splits. To disable honesty during training, you can set the parameter `honesty` to `FALSE`.

#### `mtry`

The `mtry` parameter determines the number of variables considered during each split. The value of `mtry` is often tuned as a way to improve the runtime of the algorithm, but can also have an impact on statistical performance.

By default, `mtry` is taken as `max(sqrt(p + 20), p)`, where `p` is the number of variables (columns) in the dataset. This value can be adjusted by changing the parameter `mtry` during training. Selecting a tree split is often the most resource-intensive component of the algorithm. Setting a large value for `mtry` may therefore slow down training considerably.

To more closely match the theory in the GRF paper, the number of variables considered is actually drawn from a poisson distribution with mean equal to `mtry`. A new number is sampled from the distribution before every tree split.

#### `min.node.size`

The parameter `min.node.size` relates to the minimum size a leaf node is allowed to have. Given this parameter, if a node reaches too small of a size during splitting, it will not be split further.

There are several important caveats to this parameter:
- When honesty is enabled, the leaf nodes are 'repopulated' after splitting with a fresh subsample. This means that the final tree may contain leaf nodes smaller than the `min.node.size` setting.
- For regression forests, the splitting will only stop once a node has become smaller than `min.node.size`. Because of this, trees can have leaf nodes that violate the `min.node.size` setting. We initially chose this behavior to match that of other random forest packages like `randomForest` and `ranger`, but will likely be changed as it is misleading (see [#143](https://github.com/grf-labs/grf/issues/143)).
- When training a causal forest, `min.node.size` takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Split Penalization' section below, under the 'Causal Forests' heading.

#### `alpha`

The parameter `alpha` controls the maximum imbalance of a split. In particular, when splitting a parent node, the size of each child node is not allowed to be less than `size(parent) * alpha`. Its value must lie between (0, 0.25), and defaults to 0.05.

When training a causal forest, this parameter takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Split Penalization' section below, under the 'Causal Forests' heading.

#### `imbalance.penalty`

The `imbalance.penalty` parameter controls how harshly imbalanced splits are penalized. When determining which variable to split on, each split is assigned a 'goodness measure' related to how much it increases heterogeneity across the child nodes. The algorithm applies a penalty to this value to discourage child nodes from having very different sizes, specified by `imbalance.penalty * (1.0 / size(left.child) + 1.0 / size(right.child)`. This penalty can be seen as a complement to the hard restriction on splits provided by `alpha`.

This parameter is still experimental, and unless `imbalance.penalty`is explicitly specified, it defaults to 0 so that no split penalty is applied.

When training a causal forest, this parameter takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Selecting Balanced Splits' section below, under the 'Causal Forests' heading.

### Variance Estimates

By default, all forest models are trained in such a way as to support variance estimates. To calculate these estimates, the flag `estimate.variance` can be provided to prediction:

```
causal.forest = causal_forest(X, Y, W)
prediction.result = predict(causal.forest, X.test, estimate.variance=TRUE)
standard.error = sqrt(prediction.result$variance.estimates)
```

The procedure works by training trees in small groups, then comparing the predictions within and across groups to estimate variance. In more detail:
- In each training pass, we sample the full dataset to create a subsample of half its size. Then, a small group of trees in trained on this half-sample. In particular, for each tree we draw a subsample of the half-sample, and grow the tree using these examples.
- When predicting, a variance estimate is also computed by comparing the variance in predictions within groups to the total variance. More details on the method can be found in section 4 of the GRF paper, or by examining the implementations of the C++ method `PredictionStrategy::compute_variance`.

Note that although training starts by drawing a half-sample, the `sample.fraction` option still corresponds to a fraction of the full sample. This means that when variance estimates are requested, `sample.fraction` cannot be greater than 0.5.

The number of trees in each group is controlled through the `ci.group.size` parameter, and defaults to 2. If variance estimates are not needed, `ci.group.size` can be set to 1 during training to avoid growing trees in small groups.

## Causal Forests

The `causal.forest` method uses the same general training and prediction framework described above:
- When choosing a split, the algorithm seeks to maximize the difference in treatment effect between the two child nodes. For computational efficiency, we precompute the gradient of each observation, and optimize a linear approximation of this difference.
- When predicting on a test example, we gather a weighted list of the sample's neighbors based on what leaf nodes it falls in. We then calculate the treatment effect using the outcomes and treatment status of the neighbor examples. To speed up the algorithm, we precompute certain statistics in each leaf during training, such as the average value of the treatment.

For a technical treatment of causal forest splitting and prediction, please refer to section 6.2 of the GRF paper.

Beyond this core training procedure, causal forests incorporate some additions specific to treatment effect estimation. These additions are described below.

### Orthogonalization

Recall that causal forests assume that potential outcomes are independent of treatment assignment, but only after we condition on features `X`. In this setting, in order to consistently estimate conditional average treatment effects, a naive causal forest would need to split both on features that affect treatment effects and those that affect treatment propensities. This can be wasteful, as splits 'spent' on modelling treatment propensities may not be useful in estimating treatment heterogeneity.

In GRF, we avoid this difficulty by 'orthogonalizing' our forest using Robinson's transformation (Robinson, 1988). Before running `causal_forest`, we compute estimates of the propensity scores `e(x) = E[W|X=x]` and marginal outcomes `m(x) = E[Y|X=x]` by training separate regression forests and performing out-of-bag prediction. We then compute the residual treatment `W - e(x)` and outcome `Y - m(x)`, and finally train a causal forest on these residuals. If propensity scores or marginal outcomes are known through prior means (as might be the case in a randomized trial) they can be specified through the training parameters `W.hat` and `Y.hat`. In this case, `causal_forest` will use these estimates instead of training separate regression forests.

Empirically, we've found orthogonalization to be essential in obtaining accurate treatment effect estimates in observational studies. More details on the orthogonalization procedure in the context of forests can be found in section 6.1.1 of the GRF paper. For a broader discussion on Robinson's tranformation for conditional average treatment effect estimation, including formal results, please see Nie and Wager (2017). 

### Selecting Balanced Splits

In the sections above on `min.node.size`, `alpha`, and `imbalance.penalty`, we described how tree training protects against making splits that result in a large size imbalance between the children, or leaf nodes that are too small. In a causal setting, it is not sufficient to consider the number of examples in each node -- we must also take into account the number treatment vs. control examples. Without a reasonable balance of treated and control examples, there will not be enough information in the node to obtain a good estimate of treatment effect. In the worst case, we could end up with nodes composed entirely of control (or treatment) examples.

For this reason, causal splitting uses modified notions of each split balancing parameter:
- `min.node.size` usually determines the minimum number of examples a node should contain. In causal forests, the requirement is more stringent: a node must contain at least `min.node.size` treated samples, and also at least that many control samples.
- In regression forests, `alpha` and `imbalance.penalty` help ensure that the size difference between children is not too large. When applying `alpha` and `imbalance.penalty`, we use a modified measure of node size that tries to capture how much 'information content' it contains. The new size measure is given by `\sum_{i in node} (W_i - \bar{W})^2`.

### Average Treatment Effects

In addition to personalized treatment effects, causal forests can be used to estimate the average treatment effect across the training population. Naively, one might estimate the average treatment effect by averaging personalized treatment effects across training examples. However, a more accurate estimate can be obtained by plugging causal forest predictions into a doubly robust average treatment effect estimator. As discussed in Chernozhukov et al. (2018), such approaches can yield semiparametrically efficient average treatment effect estimates and accurate standard error estimates under considerable generality. GRF provides the dedicated function `average_treatment_effect` to compute these estimates.

The `average_treatment_effect` function implements two types of doubly robust average treatment effect estimations: augmented inverse-propensity weighting (Robins et al., 1994), and targeted maximum likelihood estimation (van der Laan and Rubin, 2006). Which method to use can be specified through the `method` parameter. The parameter `target.sample` controls which population the average treatment effect is taken over:
- `target.sample = "all"`: the ATE on the whole population, `sum_{i = 1}^n E[Y(1) - Y(0) | X = X_i] / n`.
- `target.sample = "treated"`: the ATE on the treated examples, `sum_{W_i = 1} E[Y(1) - Y(0) | X = X_i] / |{i : W_i = 1}|`.
- `target.sample = "control"`: the ATE on the control examples, `sum_{W_i = 0} E[Y(1) - Y(0) | X = X_i] / |{i : W_i = 0}|`.
- `target.sample = "overlap"`: the overlap-weighted ATE `sum_{i = 1}^n e(Xi) (1 - e(Xi)) E[Y(1) - Y(0) | X = Xi] / sum_{i = 1}^n e(Xi) (1 - e(Xi))`,
  where `e(x) = P[W_i = 1 | X_i = x]`. This last estimand is recommended by Li et al. (2017) in case of poor overlap (i.e., when the treatment propensities e(x) may be very close to 0 or 1), as it doesn't involve dividing by estimated propensities.

## Additional Features

The following sections describe other features of GRF that may be of interest.

### Parameter Tuning

The accuracy of a forest can be sensitive to the choice of `min.node.size`, `sample.fraction`, and `mtry`, as well as the split balance parameters `alpha` and `imbalance.penalty`. GRF provides a cross-validation procedure to select values of these parameters to use in training. To enable this parameter tuning in training, the option `tune.parameters = TRUE` can be passed to forest method. The cross-validation methods can also be called directly through `tune_regression_forest` and`tune_causal_forest`. Parameter tuning is currently disabled by default, as we are still finalizing some details of the algorithm.

The cross-validation procedure works as follows:
- Draw a number of random points in the space of possible parameter values. By default, 100 distinct sets of parameter values are chosen.
- For each set of parameter values, train a forest with these values and compute the out-of-bag error. There are a couple important points to note about this error measure, outlined below. The exact procedure for computing the error can be found in the C++ methods that implement`OptimizedPredictionStrategy#compute_debiased_error`.
  - For tuning to be computationally tractable, we only train 'mini forests' composed of 10 trees. With such a small number of trees, the out-of-bag error gives a biased estimate of the final forest error. We therefore debias the error through a simple variance decomposition.
  - While the notion of error is straightforward for regression forests, it can be more subtle in the context of treatment effect estimation. For causal forests, we use a measure of error developed in Nie and Wager (2017) motivated by residual-on-residual regression (Robinson, 1988).
- Finally, given the debiased error estimates for each set of parameters, we apply a smoothing function to determine the optimal parameter values.

### Cluster-Robust Estimation

For accurate predictions and variance estimates, it can be important to take into account for natural clusters in the data, as might occur if a dataset contains examples taken from the same household or small town. GRF provides support for cluster-robust forests by accounting for clusters in the subsampling process. To use this feature, the forest must be trained with the relevant cluster information by specifying the `clusters` and optionally `samples_per_cluster` parameters. Then, all subsequent calls to `predict` with will take clusters into account, including when estimating the variance of forest predictions. 

When clustering is enabled during training, all subsampling procedures operate on entire clusters as opposed to individual examples. Then, to determine the examples used for performing splitting and populating the leaves, `samples_per_cluster` examples are drawn from the selected clusters. By default, `sample_per_cluster` is equal to the size of the smallest cluster. Concretely, the cluster-robust training procedure proceeds as follows:
- For each 'CI group' used in estimating variance, sample half of the clusters. Each CI group is then associated with a list of cluster IDs.
- Within this half-sample of cluster IDs, sample `sample.fraction` of the clusters. Each tree is now associated with a list of cluster IDs.
- If honesty is enabled, split these cluster IDs in half, so that one half can be used for growing the tree, and the other half is used in repopulating the leaves.
- To grow the tree, draw `samples_per_cluster` examples from each of the cluster IDs, and do the same when repopulating the leaves for honesty.

Note that when clusters are provided, standard errors from `average_treatment_effect` and `average_partial_effect` estimation are also cluster-robust.

## Troubleshooting

### GRF isn't working well on a small dataset.

If you observe poor performance on a dataset with a small number of examples, it may be worth trying out two changes:
- Disabling honesty. As noted in the section on honesty above, when honesty is enabled, the training subsample is further split in half before performing splitting. This may not leave enough information for the algorithm to determine high-quality splits.
- Skipping the variance estimate computation, by setting `ci.group.size` to 1 during training, then increasing `sample.fraction`. Because of how variance estimation is implemented, `sample.fraction` cannot be greater than 0.5 when it is enabled. If variance estimates are not needed, it may help to disable this computation and use a larger subsample size for training.

### The variance estimates are jumpy or very large.

In this case, it would be good to try growing a larger number of trees. Obtaining good variance estimates often requires growing more trees than it takes to only obtain accurate predictions.

### The causal forest method is producing nonsensical results.

If the output of the `causal.forest` method doesn't pass a sanity check based on your knowledge of the data, it may be worth checking whether the overlap assumption is violated. In order for conditional average treatment effects to be properly identified, a dataset's propensity scores must be bounded away from 0 and 1. A simple way to validate this assumption is to calculate the propensity scores by regressing the treatment assignments W against X, and examining the out-of-bag predictions. Concretely, you can perform the following steps:

```
propensity.forest = regression.forest(X, W)
W.hat = predict(propensity.forest)$predictions
hist(W.hat, xlab = "propensity score")
```

If there is strong overlap, the histogram will be concentrated away from 0 and 1. If the data is instead concentrated at the extremes, the overlap assumption likely does not hold.

For further discussion of the overlap assumption, please see Imbens and Rubin (2015). In practice, this assumption is often violated due to incorrect modelling decision: for example one covariate may be a deterministic indicator that the example received treatment.

### Regression forest predictions differ from those of the randomForest and ranger packages.

While the algorithm in `regression_forest` is very similar to that of classic random forests, it has several notable differences, including 'honesty', group tree training for variance estimates, and restrictions during splitting to avoid imbalanced child nodes. These features can cause the predictions of the algorithm to be different, and also lead to a slower training procedure than other packages. We welcome GitHub issues that shows cases where GRF does notably worse than other packages (either in statistical or computational performance), as this will help us choose better defaults for the algorithm, or potentially point to a bug.

## References

Athey, Susan, Julie Tibshirani, and Stefan Wager. Generalized Random Forests. *Annals of Statistics (forthcoming)*, 2018. 

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, and James Robins. Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 2018.

Imbens, Guido W., and Donald B. Rubin. Causal inference in statistics, social, and biomedical sciences. *Cambridge University Press*, 2015. 

Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky. Balancing covariates via propensity score weighting. *Journal of the American Statistical Association*, 2018.

Nie, Xinkun, and Stefan Wager. Learning Objectives for Treatment Effect Estimation. *arXiv preprint arXiv:1712.04912*, 2017.

Robins, James M., Andrea Rotnitzky, and Lue Ping Zhao. Estimation of regression coefficients when some regressors are not always observed. *Journal of the American statistical Association*, 1994.

Robinson, Peter M. Root-n-consistent semiparametric regression. *Econometrica*, 1988.

Van Der Laan, Mark J., and Daniel Rubin. Targeted maximum likelihood learning. *The International Journal of Biostatistics*, 2006.

Wager, Stefan, and Susan Athey. Estimation and Inference of Heterogeneous Treatment Effects using Random Forests. *Journal of the American Statistical Association (forthcoming)*, 2018.
