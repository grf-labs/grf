# The GRF Algorithm

## Table of Contents
1. [General Algorithm](#general-algorithm)
2. [Causal Forests](#causal-forests)
3. [Additional Features](#additional-features)
4. [Troubleshooting](#troubleshooting)

The following guide gives an introduction to the generalized random forests algorithm as implemented in the `grf` package. It aims to give a complete description of the training and prediction procedures, as well as the options available for tuning. This guide is intended as an informal and practical reference; for a theoretical treatment of GRF, please consult the ['Generalized Random Forests' paper](https://arxiv.org/abs/1610.01271).

GRF extends the idea of a classic random forest to allow for estimating other statistical quantities besides the expected outcome. Each forest type, for example`quantile_forest`, trains a random forest targeted at a particular problem, like quantile estimation. The most common use of GRF is in estimating treatment effects through the function `causal_forest`.

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

The main difference between GRF's approach to growing trees and that of classic random forests is in how the quality of a split is measured. Because the various forest types seek to estimate different statistical quantities like quantiles and treatment effects, splitting must be tailored to the particular task at hand. The approach taken in GRF is to maximize the heterogeneity in quantity of interest across the child nodes. For example, with causal effect estimation, the goodness of a split relates to how different the treatment effect estimates are in each node. A theoretical motivation for this split criterion can be found in section 2 of the [GRF paper](Splitting to Maximize Heterogeneity).

The quality of a split must be calculated for each possible split variable `x` and value `v`, so it is critical for it to be fast to compute. Optimizing the heterogeneity criterion directly is therefore too expensive; instead, we take the gradient of the objective and optimize a first-order approximation to the criterion. This approach is much faster, and also allows us to reuse classic forest splitting algorithms that have proven to have good performance.

### Prediction

Given a test example, the GRF algorithm computes a prediction as follows:
- For each tree, the test example is 'pushed down' to determine what leaf it falls in.
- Given this information, we create a list of neighboring training examples, weighted by how many times the example fell in the same leaf as the test example.
- A prediction is made using this weighted list of neighbors, using the relevant approach for the type of forest. For regression forests, the prediction is equal to the average outcome of the test example's neighbors. In causal prediction, we calculate the treatment effect using the outcomes and treatment status of the neighbor examples.

Those familiar with classic random forests might note that this approach differs from the way forest prediction is usually described. The traditional view is that to predict for a test example, each tree makes a prediction on that example. To make a final prediction, the tree predictions are combined in some way, for example through averaging or through 'majority voting'. It's worth noting that for regression forests, the GRF algorithm described above is identical this 'ensemble' approach, where each tree predicts by averaging the outcomes in each leaf, and predictions are combined through a weighted average.

### Out-of-bag (OOB) Prediction

If a dataset is provided to the `predict` method, then predictions are made for these new test example. When no dataset is provided, prediction proceeds on the training examples. In particular, for each training example, all the trees that did not use this example during training are identified (the example was 'out-of-bag'). Then, a prediction for the test example is made using only these trees. These out-of-bag predictions can be useful in understanding the model's goodness-of-fit, and are also used in several of the methods for causal effect estimation methods described later in this guide.

### Training Options

#### `sample.fraction`

The `sample.fraction` parameter is a number in the range (0, 1] that controls the fraction of examples that should be used in growing each tree. By default, `sample.fraction` is set to 0.5. As noted in the section on honest forests, the fractional subsample will be further split into halves when honesty is enabled.

#### `num.trees`

The parameter `num.trees` controls how many trees are grown during training, and defaults to 2000. Generally, obtaining high-quality confidence intervals requires growing more trees than are needed for accurate predictions.

Tree training is parallelized across several threads in an effort to improve performance. By default, all available cores are used, but the number of threads can be set directly through `num.threads`.

### `honesty`

By default, 'honest' forests are trained. In a classic random forest, a single subsample is used both to choose a tree's splits, and for the leaf node examples used in making predictions. In contrast, honest forests randomly split this subsample in half, and use only the first half when performing splitting. The second half is then used to populate the tree's leaf nodes: each new example is 'pushed down' the tree, and added to the leaf in which it falls. In a sense, the leaf nodes are 'repopulated' after splitting using a fresh set of examples.

The motivation behind honesty is to reduce bias in tree predictions, by using different subsamples for constructing the tree and for making predictions. Honesty is a well-explored idea in the academic literature on random forests, but is not yet common in software implementations. For a more formal overview, please see section 2.4 of ['Estimation and Inference of Heterogeneous Treatment Effects using Random Forests'](https://arxiv.org/pdf/1510.04342.pdf).

It's important to note that honesty may hurt performance when working with very small datasets. In this set-up, the subsample used to determine tree splits is already small, and honesty further cuts this subsample in half, so there may no longer be enough information to choose high-quality splits. To disable honesty during training, you can set the parameter `honesty` to `FALSE`.

#### `mtry`

The `mtry` parameter determines the number of variables considered during each split. The value of `mtry` is often tuned as a way to improve the runtime of the algorithm, but can also have an impact on statistical performance.

By default, `mtry` is taken as `max(sqrt(p + 20), p)`, where `p` is the number of variables (columns) in the dataset. This value can be adjusted by changing the parameter `mtry` during training. Selecting a tree split is often the most resource-intensive component of the algorithm. Setting a large value for `mtry` may therefore slow down training considerably.

To more closely match the theory in the GRF paper, the number of variables considered is actually drawn from a poisson distribution with mean equal to `mtry`. A new number is sampled from the distribution before every tree split.

#### `min.node.size`

The parameter `min.node.size` relates to the minimum size a leaf node is allowed to have. Given this parameter, if a node reaches too small of a size during splitting, it will not be split further.

There are several important caveats to this parameter:
- When honesty is enabled, the leaf nodes are 'repopulated' after splitting with a fresh subsample. This means that the final tree may contain leaf nodes smaller than the `min.node.size` setting.
- For regression forests, the splitting will only stop once a node has become smaller than `min.node.size`. Because of this, trees can have leaf nodes that violate the `min.node.size` setting. We initially chose this behavior to match that of other random forest packages like `randomForest` and `ranger`, but will likely be changed as it is misleading (see [#143](https://github.com/swager/grf/issues/143)).
- When training a causal forest, `min.node.size` takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Split Penalization' section below, under the 'Causal Forests' heading.

#### `alpha`

The parameter `alpha` controls the maximum imbalance of a split. In particular, when splitting a parent node, the size of each child node is not allowed to be less than `size(parent) * alpha`. Its value must lie between (0, 0.25), and defaults to 0.05.

When training a causal forest, this parameter takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Split Penalization' section below, under the 'Causal Forests' heading.

#### `imbalance.penalty`

The `imbalance.penalty` parameter controls how harshly imbalanced splits are penalized. When determining which variable to split on, each split is assigned a 'goodness measure' related to how much it increases heterogeneity across the child nodes. The algorithm applies a penalty to this value to discourage child nodes from having very different sizes, specified by `imbalance.penalty * (1.0 / size(left.child) + 1.0 / size(right.child)`. This penalty can be seen as a complement to the hard restriction on splits provided by `alpha`.

This parameter is still experimental, and unless `imbalance.penalty`is explicitly specified, it defaults to 0 so that no split penalty is applied.

When training a causal forest, this parameter takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Split Penalization' section below, under the 'Causal Forests' heading.

### Variance Estimates

By default, all forest models are trained in such a way as to support variance estimates. To calculate these estimates, the flag `estimate.variance` can be provided to prediction:

```
causal.forest = causal_forest(X, Y, W)
prediction.result = predict(causal.forest, X.test, estimate.variance=TRUE)
variance.estimates = prediction.result$variance.estimates
```

The procedure works by training trees in small groups, then comparing the predictions within and across groups to estimate variance. In more detail:
- In each training pass, we sample the full dataset to create a subsample of half its size. Then, a small group of trees in trained on this half-sample. In particular, for each tree we draw a subsample of the half-sample, and grow the tree using these examples.
- When predicting, a variance estimate is also computed by comparing the variance in predictions within groups to the total variance. More details on the method can be found in section 4 of the [GRF paper](https://arxiv.org/pdf/1610.01271.pdf), or by examining the implementations of the C++ method `PredictionStrategy::compute_variance`.

Note that although training starts by drawing a half-sample, the `sample.fraction` option still corresponds to a fraction of the full sample. This means that when variance estimates are requested, `sample.fraction` cannot be greater than 0.5.

The number of trees in each group is controlled through the `ci.group.size` parameter, and defaults to 2. If variance estimates are not needed, `ci.group.size` can be set to 1 during training to avoid growing trees in small groups.

## Causal Forests

### Local Centering (TODO)

### Split Penalization (TODO)

### Average Treatment Effects (TODO)

## Additional Features

The following sections describe other features of GRF that may be of interest.

### Parameter Tuning (TODO)

### Clustering (TODO)

## Troubleshooting

### GRF isn't working well on a small dataset.

If you observe poor performance on a dataset with a small number of examples, it may be worth trying out two changes:
- Disabling honesty. As noted in the section on honesty above, when honesty is enabled, the training subsample is further split in half before performing splitting. This may not leave enough information for the algorithm to determine high-quality splits.
- Skipping the variance estimate computation, by setting `ci.group.size` to 1 during training, then increasing `sample.fraction`. Because of how variance estimation is implemented, `sample.fraction` cannot be greater than 0.5 when it is enabled. If variance estimates are not needed, it may help to disable this computation and use a larger subsample size for training.

### The variance estimates are jumpy or very large.

In this case, it would be good to try growing a larger number of trees. Obtaining good variance estimates often requires growing more trees than it takes to only obtain accurate predictions.

### The predictions for `regression_forest` differ from those of `randomForest` or `ranger`.

While the algorithm in `regression_forest` is very similar to that of classic random forests, it has several notable differences, including 'honesty', group tree training for variance estimates, and restrictions during splitting to avoid imbalanced child nodes. These features can cause the predictions of the algorithm to be different, and also lead to a slower training procedure than other packages. We welcome GitHub issues that shows cases where GRF does notably worse than other packages (either in statistical or computational performance), as this will help us choose better defaults for the algorithm, or potentially point to a bug.

