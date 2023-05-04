# The GRF Algorithm

<img src='https://raw.githubusercontent.com/grf-labs/grf/master/images/logo/grf_logo_wbg_cropped.png' align="right" height="120" />

The following guide gives an introduction to the generalized random forests algorithm as implemented in the `grf` package. It aims to give a complete description of the training and prediction procedures, as well as the options available for tuning. This guide is intended as an informal and practical reference; for a theoretical treatment of GRF, please consult the 'Generalized Random Forests' [paper](https://arxiv.org/abs/1610.01271).

GRF extends the idea of a classic random forest to allow for estimating other statistical quantities besides the expected outcome. Each forest type, for example `quantile_forest`, trains a random forest targeted at a particular problem, like quantile estimation. The most common use of GRF is in estimating treatment effects through the function `causal_forest`.

## Table of Contents
* [General Algorithm](#general-algorithm)
  * [Training](#training)
  * [Prediction](#prediction)
  * [Out-of-bag Prediction](#out-of-bag-prediction)
  * [Training Options](#training-options)
  * [Variance Estimates](#variance-estimates)
* [Causal Forests](#causal-forests)
  * [Orthogonalization (the R-learner)](#orthogonalization-the-r-learner)
  * [Selecting Balanced Splits](#selecting-balanced-splits)
  * [Average Treatment Effects](#average-treatment-effects)
  * [Rank-Weighted Average Treatment Effects](#rank-weighted-average-treatment-effects)
  * [Best Linear Projection of the CATE](#best-linear-projection-of-the-cate)
  * [Multiple Outcomes and Multiple Treatments](#multiple-outcomes-and-multiple-treatments)
  * [Right-Censored Survival Outcomes](#right-censored-survival-outcomes)
* [Additional Features](#additional-features)
  * [Parameter Tuning](#parameter-tuning)
  * [Merging Forests](#merging-forests)
  * [Boosted Regression Forests](#boosted-regression-forests)
  * [Cluster-Robust Estimation](#cluster-robust-estimation)
  * [Sample Weighting](#sample-weighting)
  * [Categorical Inputs](#categorical-inputs)
  * [Missing values](#missing-values)
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

The following table illustrates with a simple regression example. Given a forest trained on n = 4 training samples, we are given a new test point *x* we want to compute a prediction at. We consider all the training samples used to build the forest, and ask for each tree: for which training samples does the new *x* land in the same terminal node? As a table this would look like:

| Training sample 	| Does test sample  *x*  fall into the same leaf? 	|         	|     	|         	| "Count Yes"	| Similarity weight for *x*: &alpha;(*x*)| Y (outcome)      	| &alpha; * Y           |
|-----------------	|-------------------------------------------------	|---------	|-----	|---------	|-----------	|-------------------------------------	|--------------------	|---------------------	|
|                 	| *tree 1*                                         	| *tree 2* 	| ... 	| *tree B* 	|           	|                                     	|                    	|                     	|
| X<sub>1</sub>    	| **Yes**                                         	| No      	|     	| **Yes** 	| 1/3+0+1   	| 4/9                                 	| 4                  	| 16/9                 	|
| X<sub>2</sub>    	| No                                              	| No      	|     	| No      	| 0+0+0     	| 0                                   	| -0.5               	| 0                   	|
| X<sub>3</sub>    	| **Yes**                                         	| No      	|     	| No      	| 1/3+0+0   	| 1/9                                 	| 2                  	| 2/9                 	|
| X<sub>4</sub>    	| **Yes**                                         	| **Yes** 	|     	| No      	| 1/3+1+0   	| 4/9                                 	| 4                  	| 16/19                	|
|                 	|                                                 	|         	|     	|         	|           	|                                     	| *Final prediction:*	| 34/9                  |

where in the right column we sum each training sample's weighted contribution to a "Yes", then normalize these to sum to one. This effectively acts as a similarity measure telling us how "similar" each training sample is to *x*. Combining these forest weights &alpha;(*x*) with the estimation method given by the particular task at hand, is what delivers the final point predictions. For regression, the estimation method is a mean, giving us a forest weighted average (in a "traditional" random forest *tree 1* would use sample *1*, *3* and *4* to predict an average Y, etc., and these averages would be aggregated for the final prediction). For other statistical tasks, the final "aggregation" could be more complicated, for example, when using a random forest for instrumental variables estimation, GRF would effectively compute a forest-weighted two-stage least squares estimate for each target sample *x*.

Some readers may find this description similar to kernel weighted estimators, and that is right: GRF uses random forests as an adaptive nearest neighbor method. The GRF algorithm finds a weighted set of neighbors &alpha;(*x*) that are similar to a test point *x*, where the notion of similarity varies by the statistical task. For regression, it is how similar the outcomes are, for CATE estimation, it is how similar the estimated treatment effects are. An important consequence of this is that asking for individual tree predictions, as one would for a traditional ensemble random forest, is inappropriate here, as it is not individual trees that make up predictions, but rather the neighborhood they induce.

### Out-of-bag Prediction

If a dataset is provided to the `predict` method, then predictions are made for these new test example. When no dataset is provided, prediction proceeds on the training examples. In particular, for each training example, all the trees that did not use this example during training are identified (the example was 'out-of-bag', or OOB). Then, a prediction for the test example is made using only these trees. These out-of-bag predictions can be useful in understanding the model's goodness-of-fit, and are also used in several of the methods for causal effect estimation methods described later in this guide.

### Training Options

#### `sample.fraction`

The `sample.fraction` parameter is a number in the range (0, 1] that controls the fraction of examples that should be used in growing each tree. By default, `sample.fraction` is set to 0.5. As noted in the section on honest forests, the fractional subsample will be further split into halves when honesty is enabled.

#### `num.trees`

The `num.trees` parameter controls how many trees are grown during training, and defaults to 2000. Tree training is parallelized across several threads in an effort to improve performance. By default, all available cores are used, but the number of threads can be set directly through `num.threads`.

Forests are a randomized ensemble algorithm, and as such every forest grown with a different initial seed will produce slightly different estimates, even when fit on the same data. We call this sort of perturbation error the `excess.error`, because it does not come from the inherent sampling variability in the data. Large forests have smaller `excess.error`, so we recommend that users grow as many trees as necessary to ensure that the values in `excess.error` are negligible relative to `variance.estimates`.

In addition, obtaining tighter confidence intervals requires growing even more trees than are needed for accurate predictions. When the number of trees in a forest is small, the confidence intervals can be too wide, and therefore too conservative. We recommend that users grow trees in proportion to the number of observations.

If you are interested in checking the evolution of `excess.error` or confidence interval widths as the number of trees increases, you can grow forests iteratively using the function `merge_forests` (see [Merging Forests](#merging-forests) below).


#### `honesty`, `honesty.fraction`, `honesty.prune.leaves`

By default, 'honest' forests are trained. The motivation behind honesty is to reduce bias in tree predictions, by using different subsamples for constructing the tree and for making predictions. Honesty is a well-explored idea in the academic literature on random forests, but is not yet common in software implementations. This section gives an algorithmic explanation of honesty; for a more formal overview, please see section 2.4 of Wager and Athey (2018).

In a classic random forest, a single subsample is used both to choose a tree's splits, and for the leaf node examples used in making predictions. In contrast, honest forests randomly split this subsample in half, and use only the first half when performing splitting. The second half is then used to populate the tree's leaf nodes: each new example is 'pushed down' the tree, and added to the leaf in which it falls. In a sense, the leaf nodes are 'repopulated' after splitting using a fresh set of examples.

After repopulating a tree's leaves using the second half-sample, it is possible that some leaves end up empty. With empty leaves, a tree is not able to handle certain test examples and needs to be skipped when computing those predictions. By default, empty leaves are pruned away after training so that each tree is able to handle all test points. GRF's behavior with respect to empty leaves can be controlled through the parameter `honesty.prune.leaves`.

It's important to note that honesty may hurt performance when working with very small datasets. In this set-up, the subsample used to determine tree splits is already small, and honesty further cuts this subsample in half, so there may no longer be enough information to choose high-quality splits. There are a couple options for mitigating the cost of honesty in small samples:
- The parameter `honesty.fraction` allows for increasing the fraction of samples used in selecting tree splits. For example, an honesty fraction of 0.7 directs GRF to use 70% of the tree subsample for splitting, and the other 30% to populate the leaf nodes. If GRF is not working well on a small sample, we've found empirically that it can help to increase `honesty.fraction` and set `honesty.prune.leaves` to `FALSE`. With these settings, it may also be necessary to increase the number of trees grown in training through `num.trees`.
- Alternatively, honesty can be completely disabled during training by setting the parameter `honesty` to `FALSE`.

When automatically tuning parameters through `tune.parameters`, GRF will try varying `honesty.fraction` between 0.5 and 0.8, and consider both options for `honesty.prune.leaves`. More information on automatic parameter selection can be found in the 'Parameter Tuning' section below.

#### `mtry`

The `mtry` parameter determines the number of variables considered during each split. The value of `mtry` is often tuned as a way to improve the runtime of the algorithm, but can also have an impact on statistical performance.

By default, `mtry` is taken as `min(sqrt(p) + 20, p)`, where `p` is the number of variables (columns) in the dataset. This value can be adjusted by changing the parameter `mtry` during training. Selecting a tree split is often the most resource-intensive component of the algorithm. Setting a large value for `mtry` may therefore slow down training considerably.

To more closely match the theory in the GRF paper, the number of variables considered is actually drawn from a poisson distribution with mean equal to `mtry`. A new number is sampled from the distribution before every tree split.

#### `min.node.size`

The parameter `min.node.size` relates to the minimum size a leaf node is allowed to have. Given this parameter, if a node reaches too small of a size during splitting, it will not be split further.

There are several important caveats to this parameter:
- When honesty is enabled, the leaf nodes are 'repopulated' after splitting with a fresh subsample. This means that the final tree may contain leaf nodes smaller than the `min.node.size` setting.
- For regression forests, the splitting will only stop once a node has become smaller than `min.node.size`. Because of this, trees can have leaf nodes that violate the `min.node.size` setting. We initially chose this behavior to match that of other random forest packages like `randomForest` and `ranger`, but will likely be changed as it is misleading (see [#143](https://github.com/grf-labs/grf/issues/143)).
- When training a causal forest, `min.node.size` takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Selecting Balanced Splits' section below, under the 'Causal Forests' heading.

#### `alpha`

The parameter `alpha` controls the maximum imbalance of a split. In particular, when splitting a parent node, the size of each child node is not allowed to be less than `size(parent) * alpha`. Its value must lie between (0, 0.25), and defaults to 0.05.

When training a causal forest, this parameter takes on a slightly different notion related to the number of treatment and control samples. More detail can be found in the 'Selecting Balanced Splits' section below, under the 'Causal Forests' heading.

#### `imbalance.penalty`

The `imbalance.penalty` parameter controls how harshly imbalanced splits are penalized. When determining which variable to split on, each split is assigned a 'goodness measure' related to how much it increases heterogeneity across the child nodes. The algorithm applies a penalty to this value to discourage child nodes from having very different sizes, specified by `imbalance.penalty * (1.0 / size(left.child) + 1.0 / size(right.child)`. This penalty can be seen as a complement to the hard restriction on splits provided by `alpha`. It defaults to 0 so that no split penalty is applied.

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

### Orthogonalization (the R-learner)

Recall that causal forests assume that potential outcomes are independent of treatment assignment, but only after we condition on features `X`. In this setting, in order to consistently estimate conditional average treatment effects, a naive causal forest would need to split both on features that affect treatment effects and those that affect treatment propensities. This can be wasteful, as splits 'spent' on modeling treatment propensities may not be useful in estimating treatment heterogeneity.

In GRF, we avoid this difficulty by 'orthogonalizing' our forest using Robinson's transformation (Robinson, 1988). Before running `causal_forest`, we compute estimates of the propensity scores `e(x) = E[W|X=x]` and marginal outcomes `m(x) = E[Y|X=x]` by training separate regression forests and performing out-of-bag prediction. We then compute the residual treatment `W - e(x)` and outcome `Y - m(x)`, and finally train a causal forest on these residuals. If propensity scores or marginal outcomes are known through prior means (as might be the case in a randomized trial) they can be specified through the training parameters `W.hat` and `Y.hat`. In this case, `causal_forest` will use these estimates instead of training separate regression forests.

Empirically, we've found orthogonalization to be essential in obtaining accurate treatment effect estimates in observational studies. This first-stage residualization is conceptually different from the GRF splitting routine described in the previous sections. It can be interpreted as giving rise to a loss function tailored for heterogeneous treatment effects. In recognition of the work of Robinson (1988) and to emphasize the role of residualization, we refer to this approach as the 'R-learner'. More details on the orthogonalization procedure in the context of forests can be found in section 6.1.1 of the GRF paper. For a broader discussion on Robinson's transformation for conditional average treatment effect estimation, including formal results, please see Nie and Wager (2021).

### Selecting Balanced Splits

In the sections above on `min.node.size`, `alpha`, and `imbalance.penalty`, we described how tree training protects against making splits that result in a large size imbalance between the children, or leaf nodes that are too small. In a causal setting, it is not sufficient to consider the number of examples in each node -- we must also take into account the number treatment vs. control examples. Without a reasonable balance of treated and control examples, there will not be enough information in the node to obtain a good estimate of treatment effect. In the worst case, we could end up with nodes composed entirely of control (or treatment) examples.

For this reason, causal splitting uses modified notions of each split balancing parameter:
- `min.node.size` usually determines the minimum number of examples a node should contain. In causal forests, the requirement is more stringent: a node must contain at least `min.node.size` treated samples, and also at least that many control samples.
- In regression forests, `alpha` and `imbalance.penalty` help ensure that the size difference between children is not too large. When applying `alpha` and `imbalance.penalty` in causal forests, we use a modified measure of node size that tries to capture how much 'information content' it contains. The new size measure is given by `\sum_{i in node} (W_i - \bar{W})^2`.

The above description of `min.node.size` assumes that the treatment is binary, which in most cases is an oversimplification. The precise algorithm for enforcing `min.node.size` is as follows. Note that this same approach is used both when the treatment is binary or continuous.
- Take the average of the parent node's treatment values.
- When considering a split, require that each child node have `min.node.size` samples with treatment value less than the average, and at least that many samples with treatment value greater than or equal to the average.

### Average Treatment Effects

In addition to personalized treatment effects, causal forests can be used to estimate the average treatment effect across the training population. Naively, one might estimate the average treatment effect by averaging personalized treatment effects across training examples. However, a more accurate estimate can be obtained by plugging causal forest predictions into a doubly robust average treatment effect estimator. As discussed in Chernozhukov et al. (2018), such approaches can yield semiparametrically efficient average treatment effect estimates and accurate standard error estimates under considerable generality. GRF provides the dedicated function `average_treatment_effect` to compute these estimates.

The `average_treatment_effect` function implements two types of doubly robust average treatment effect estimations: augmented inverse-propensity weighting "AIPW" (Robins et al., 1994), and targeted maximum likelihood estimation (van der Laan and Rubin, 2006). Which method to use can be specified through the `method` parameter. The following estimates are available:
- The average treatment effect (`target.sample = all`): `E[Y(1) - Y(0)]`.
- The average treatment effect on the treated (`target.sample = treated`): `E[Y(1) - Y(0) | Wi = 1]`.
- The average treatment effect on the controls (`target.sample = control`): `E[Y(1) - Y(0) | Wi = 0]`.
- The overlap-weighted average treatment effect (`target.sample = overlap`): `E[e(X) (1 - e(X)) (Y(1) - Y(0))] / E[e(X) (1 - e(X))`, where `e(x) = P[Wi = 1 | Xi = x]`.

This last estimand is recommended by Li et al. (2018) in case of poor overlap (i.e., when the treatment propensities e(x) may be very close to 0 or 1), as it doesn't involve dividing by estimated propensities.

### Rank-Weighted Average Treatment Effects

Even though there is treatment effect heterogeneity across the population or a sub-population, the average treatment effect might still be zero. To assess if the estimated CATE function `tau(x) = E[Y(1) - Y(0) | X = x]` does well in identifying personalised effects, another summary measure would be useful. One such proposed measure is the Rank-Weighted Average Treatment Effect (RATE), which uses the relative ranking of the estimated CATEs to gauge if it can effectively target individuals with high treatment effects on a separate evaluation data set (and can thus be used to test for the presence of heterogeneous treatment effects).

This approach is implemented in the function `rank_average_treatment_effect`, which provides valid bootstrapped errors of the RATE, along with a Targeting Operator Characteristic (TOC) curve which can be used to visually inspect how well a CATE estimator performs in ordering observations according to treatment benefit. For more details on the RATE metric see Yadlowsky et al., 2021.

### Best Linear Projection of the CATE

Sometimes, it is helpful to have a more expressive summary of the CATE function `tau(x) = E[Y(1) - Y(0) | X = x]` than the average treatment effect. One useful summary of this type is the "best linear projection" of tau(x) onto a set of features A, i.e., the population solution to the linear model

`tau(X) ~ beta0 + A * beta`

Qualitatively, these betas can be used to assess the association of the CATE function with the features A. The features A could, for example, be a subset of the features X used to train the forest. Note that, if the features A are mean zero, the intercept beta0 corresponds to the average treatment effect.

The function `best_linear_projection` provides estimates of the coefficients beta in the above model. Note that this function does not simply regress the CATE estimates tau.hat(Xi) given by the causal forest against Ai. Instead, we use a doubly robust estimator that generalizes the augmented inverse-propensity weighted estimator for the average treatment effect described above.

### Multiple Outcomes and Multiple Treatments

As mentioned above, causal forest choose splits that maximize the difference in the treatment effect tau between two child nodes by a gradient based approximation. When the response `Yi` and treatment `Wi` are real-valued this results in a scalar valued approximation &rho;i which is used by the subsequent splitting routine (see section 2.3 of the GRF paper for more details).

There are however empirical applications where there are more than one primary outcome `Y` of interest, and a practitioner may wish to train a forest that jointly expresses heterogeneity in treatment effects across several outcomes `Yj`, j=1...M. The same algorithmic framework described above can be employed for this purpose (assuming all outcomes are on the same scale) by computing the gradient based approximation of tau_j across each outcome, then concatenating them to a M-sized pseudo response vector &rho; which is used by a subsequent CART regression splitting routine. The CART criterion GRF uses for vector-valued responses is the squared L2 norm.

Another closely related practical application are settings where there are several interventions, as for example in medical trials with multiple treatment arms. In the event there are K mutually exclusive treatment choices, we can use the same algorithmic principles described above to build a forest that jointly targets heterogeneity across the K-1 different treatment contrasts. In the software package this is done with (20) from the GRF paper, where Wi is a vector encoded as {0, 1}^(K-1), and &xi; selects the K-1 gradient approximations of the contrasts.

The functionality described above is available in `multi_arm_causal_forest`. The intended use-case for this function is for a "handful" of arms or outcomes. Statistical considerations aside (s.t. limited overlap with many arms), our implementation has a memory scaling that for each tree takes the form `O(number.of.nodes * M * K^2)`. The general case where all `Wk` are real-valued is implemented in the conditionally linear model forest: `lm_forest`.

### Right-Censored Survival Outcomes

Many applications featuring time-to-event data involve right-censored responses where instead of observing the survival time `Ti` we observe `Yi = min(Ti, Ci)` along with an event indicator `Di = 1{Ti <= Ci}`. To estimate treatment effects in such a setting we need to account for censoring in order to obtain unbiased estimates.

GRF supports this use-case in the function `causal_survival_forest` by incorporating an extension of the celebrated AIPW score that is robust to censoring to the orthogonalization approach described in the previous section (see Cui et al., 2020 for more details). In addition to estimates of the propensity score, this approach relies on estimates of the survival functions `S(t, x) = P(T > t | X = x)` and `C(t, x) = P(C > t | X = x)`, which GRF estimates using a variant of random survival forests (Ishwaran et al., 2008), where the notable difference is GRF uses honest splitting and forest weights to produce a kernel-weighted survival function, instead of aggregating terminal node survival curves.

The balanced split criterions described in the earlier paragraph are extended to take censoring into account. For both `causal_survival_forest` and `survival_forest` the `alpha` parameter controls the minimum number of non-censored samples each child node needs for splitting to proceed.

## Additional Features

The following sections describe other features of GRF that may be of interest.

### Parameter Tuning

The accuracy of a forest can be sensitive to several training parameters:
- the core options tree-growing options `min.node.size`, `sample.fraction`, and `mtry`
- the parameters that control honesty behavior `honesty.fraction` and `honesty.prune.leaves`
- the split balance parameters `alpha` and `imbalance.penalty`

GRF provides a cross-validation procedure to select values of these parameters to use in training. To enable this tuning during training, the option `tune.parameters = "all"` can be passed to main forest method. Parameter tuning is currently disabled by default.

The cross-validation procedure works as follows:
- Draw a number of random points in the space of possible parameter values. By default, 100 distinct sets of parameter values are chosen (`tune.num.reps`).
- For each set of parameter values, train a forest with these values and compute the out-of-bag error. There are a couple important points to note about this error measure, outlined below. The exact procedure for computing the error can be found in the C++ methods that implement `OptimizedPredictionStrategy#compute_debiased_error`.
  - For tuning to be computationally tractable, we only train 'mini forests' composed of a small number of trees (`tune.num.trees`). With such a small number of trees, the out-of-bag error gives a biased estimate of the final forest error. We therefore debias the error through a simple variance decomposition.
  - While the notion of error is straightforward for regression forests, it can be more subtle in the context of treatment effect estimation. For causal forests, we use a measure of error developed in Nie and Wager (2021) motivated by residual-on-residual regression (Robinson, 1988).
- Finally, given the debiased error estimates for each set of parameters, we apply a smoothing function to determine the optimal parameter values. The optimal parameters are the ones minimizing the predicted smoothed error on a new random draw of possible parameter values (of size `tune.num.draws`).

Note that `honesty.fraction` and `honesty.prune.leaves` are only considered for tuning when `honesty = TRUE` (its default value). Parameter tuning does not try different options of `honesty` itself.

Tuning can be sensitive. Tuning only a subset of the parameters (with for example `tune.parameters = c("min.node.size", "honesty.prune.leaves")` may give better performance than trying to find the optimum for all. On smaller data, increasing `tune.num.trees` may not necessarily increase training time drastically and could result in more stable results.

### Merging Forests

In order to ensure valid predictions and tight confidence intervals, users may need to grow a large number of trees. However, it is hard to know exactly how many trees to grow in advance. That is why
GRF allows users to grow their forests separately and then create a single forest using the function `merge_forests`. This functionality allows users to sequentially grow small forests, merge them into a large forest, and check if they have attained the desired level of `excess.error` or tight enough confidence intervals.


### Boosted Regression Forests

Ghosal and Hooker (2018) show how a boosting step can reduce bias in random forest predictions. We provide a `boosted_regression_forest` method which trains a series of forests, each of which are trained on the OOB residuals from the previous step. `boosted_regression_forest` includes the same parameters as `regression_forest`, which are passed directly to the forest trained in each step.

The `boosted_regression_forest` method also contains parameters to control the step selection procedure. By default, the number of boosting steps is automatically chosen through the following cross-validation procedure:
- First, a single `regression_forest` is trained and added to the boosted forest.
- To decide if another step should be taken, we train a small forest of size `boost.trees.tune` on the OOB residuals. If it's estimated OOB error reduces the OOB error from the previous step by more than `boost.error.reduction`, then we will prepare for another step.
- To take the next step, we train a full-sized `regression_forest` on the residuals and add it to the boosted forest.
- This process continues until `boost.error.reduction` cannot be met when training the small forest, or if the total number of steps exceeds the limit `boost.max.steps`.

Alternatively, you can to skip the cross-validation procedure and specify the number of steps directly through the parameter `boost.steps`.

Some additional notes about the behavior of boosted regression forests:
- For computational reasons, if `tune.parameters` is enabled, then parameters are chosen by the `regression_forest` procedure once in the first boosting step. The selected parameters are then applied to train the forests in any further steps.
- The `estimate.variance` parameter is not available for boosted forests.
- OOB predictions are available for the training data, which combine the OOB predictions for the forests in each boosting step.
- We have found that boosting improves out-of-bag forest predictions most in scenarios where there is a strong signal-to-noise ratio.

### Cluster-Robust Estimation

For accurate predictions and variance estimates, it can be important to take into account for natural clusters in the data, as might occur if a dataset contains examples taken from the same household or small town. GRF provides support for cluster-robust forests by accounting for clusters in the subsampling process. To use this feature, the forest must be trained with the relevant cluster information by specifying the `clusters` and optionally the `equalize.cluster.weights` parameter. Then, all subsequent calls to `predict` will take clusters into account, including when estimating the variance of forest predictions.

When clustering is enabled during training, all subsampling procedures operate on entire clusters as opposed to individual examples. Then, to determine the examples used for performing splitting and populating the leaves, `samples_per_cluster` examples are drawn from the selected clusters. By default, `samples_per_cluster` is all the observations in the clusters. If `equalize.cluster.weights = TRUE` `samples_per_cluster` is equal to the size of the smallest cluster.
Concretely, the cluster-robust training procedure proceeds as follows:
- If `estimate.variance` is enabled, sample half of the clusters. Each 'tree group' is associated with this half-sample of clusters (as described in the 'Variance Estimates' section). If we are not computing variance estimates, then keep the full sample instead.
- Within this set of cluster IDs, sample `sample.fraction` of the clusters. Each tree is now associated with a list of cluster IDs.
- If honesty is enabled, split these cluster IDs in half, so that one half can be used for growing the tree, and the other half is used in repopulating the leaves.
- To grow the tree, draw `samples_per_cluster` examples from each of the cluster IDs, and do the same when repopulating the leaves for honesty.

Note that when clusters are provided, standard errors from `average_treatment_effect` estimation are also cluster-robust. Moreover, If clusters are specified, then each unit gets equal weight by default. For example, if there are 10 clusters with 1 unit each and per-cluster ATE = 1, and there are 10 clusters with 19 units each and per-cluster ATE = 0, then the overall ATE is 0.05 (additional sample.weights allow for custom weighting). If equalize.cluster.weights = TRUE each cluster gets equal weight and the overall ATE is 0.5.

### Sample Weighting

When the distribution of data that you observe is not representative of the population you are interested in scientifically, it can be important to adjust for this. We can pass `sample.weights` to specify that in our population of interest, we observe Xi with probability proportional to `sample.weights[i]`. By default, these weights are constant, meaning that our population of interest is the population from which X1 ... Xn are sampled. For causal validity, the weights we use should not be confounded with the potential outcomes --- typically this is done by having them be a function of Xi. One common example is inverse probability of complete case weighting to adjust for missing data, which allows us to work with only the complete cases (the units with nothing missing) to estimate properties of the full data distribution (all units as if nothing were missing).

The impact these weights have depends on what we are estimating. When our estimand is a function of x, e.g. `r(x) = E[Y | X=x]` in `regression_forest` or `tau(x)=E[Y(1)-Y(0)| X=x]` in `causal_forest`, passing weights that are a function of x does not change our estimand. It instead prioritizes fit on our population of interest. In `regression_forest`, this means minimizing weighted mean squared error, i.e. mean squared error over the population specified by the sample weights. In `causal_forest`, this means minimizing weighted R-loss (Nie and Wager (2021), Equations 1/2). When our estimand is an average of such a function, as in `average_treatment_effect`, it does change the estimand. Our estimand will be the average treatment/partial effect over the population specified by our sample weights.

Currently most forest types takes `sample.weights` into account during splitting and prediction. The exception is local linear forests and `quantile_forest` where sample weighting is currently not a feature. `survival_forest` only takes `sample.weights` into account during prediction (due to properties of the log-rank splitting criterion not being "amendable" to sample weights).

### Categorical inputs

GRF can only handle numerical inputs. If your data has a column with categorical values, here are some ways in which you can proceed.

If there is a natural way to order the categories, then simply encode them as integers according to this ordering. For instance, an individual's educational attainment could be encoded as 1 for "primary", 2 for "secondary", 3 for "undergraduate", and so on. This representation will make it easy for forests to split individuals into lower and higher education groups. Moreover, even if your variable does not have an obvious ordering at first sight, you may want to be a little creative in imposing an artificial ordering that is relevant to your statistical problem. For example, numbering contiguous geographical areas sequentially can often help forests make splits that group nearby areas together. If geographical proximity is important, the forest performance will likely increase.

However, sometimes the category has no reasonable ordering. In that case, one common option is to create one-hot encoded vectors (also know as dummy vectors), replacing the column with several binary columns indicating the category to which each observation belongs. Unfortunately, this sort of representation can be wasteful, because each binary column will likely not carry a lot of information. When the number of categories is large, you might want to use a more sophisticated representation method such as feature hashing. Alternatively, we recommend users try our sister package [`sufrep`](https://github.com/grf-labs/sufrep), which implements several categorical variable representation methods. It is currently in beta version, so please feel free to submit github issues for bugs or additional feature requests.

Finally, please bear in mind that in general the statistical performance of any statistical algorithm will depend on how the data is represented. Therefore, it may be worthwhile to spend some time thinking about how to encode your variables before starting the actual analysis.

### Missing values

GRF can handle missing covariate values, meaning that a the input covariate matrix can contain `NA` in some cells instead of numerical values.

These missing values can be _informative_, in the sense that the fact that a covariate is missing can be predictive of the unobserved value itself or of the outcome Y, or it can be _non-informative_.
For instance, an individual's annual income could be missing because this individual is wealthy and prefers to conceal his/her income and this in turn is likely to be predictive, for instance, of the type of housing of this individual.

GRF handles missing values implicitly using the missing incorporated in attributes criterion (MIA; Twala et al., 2008), meaning that as soon as there is a missing value in a variable _j_, there are three candidate splits for a given threshold.
- split the samples with observed value for _j_ according to the threshold and send all samples with missing value for _j_ to the __left__,
- split the samples with observed value for _j_ according to the threshold and send all samples with missing value for _j_ to the __right__,
- send all samples with observed value for _j_ to the left and all samples with missing value for _j_ to the right.

Note that by default, missing values are sent to the left, which implies that if the covariates from the training set are completely observed but there are missing values occurring in a new covariate matrix from the test set, then the prediction functions will not produce an error but simply send the incomplete observations to the left side of the corresponding nodes.

Missing values occur almost inevitably in most applications and the way to handle them depends on the type of analysis one wishes to perform and on the mechanism that generated these missing values.
For a broader discussion on treatment effect estimation with missing attributes please see Mayer et al. (2019).

## Troubleshooting

### GRF isn't working well on a small dataset.

If you observe poor performance on a dataset with a small number of examples, there are two changes worth looking into:
- Adjusting honesty parameters during training. When honesty is enabled, the training subsample is further split in half before performing splitting. This may not leave enough information for the algorithm to determine high-quality splits. The section above on the `honesty` parameter gives guidance on how to mitigate the issue.
- Skipping the variance estimate computation by setting `ci.group.size` to 1 during training, then increasing `sample.fraction`. Because of how variance estimation is implemented, `sample.fraction` cannot be greater than 0.5 when it is enabled. If variance estimates are not needed, it may help to disable this computation and use a larger subsample size for training.

### GRF isn't working well on a massive dataset.

GRF is a package designed for forest-based estimation on datasets of 'moderate' size. The tree building process can be sped up significantly on a machine with many cores, but note that on very large datasets the problem quickly becomes memory bound, due to the information (sufficient statistics, etc.) stored by each tree. As a point of reference a `causal_forest` with 1 000 000 observations and 30 continuous covariates takes around 1 hour to train (using default settings on a machine with 24 cores and 150 GB RAM), and the final forest takes up ca. 25 GB of memory (this is expected to scale linearly in the number of trees). In order to make forest training on very large data more manageable, we provide the following recommendations for parameters to experiment with:
- Increase `min.node.size`, effectively growing shallower trees. Having a `min.node.size` in the thousands will make scaling to a data set in the millions perfectly feasible.
- Decrease the number of trees, `num.trees`. The default of 2000 is high in order to ensure good performance of the optionally estimated pointwise confidence intervals. If these are not needed, fewer trees might give very similar prediction performance.

In addition, the following might also be useful, depending on application:
- Decrease the subsampling rate `sample.fraction`, reducing the number of samples per tree.
- Split up training and take a union of forests on a machine with more RAM using `merge_forests` (note that merging is still a copy inducing operation).
- For a moderate gain in speed, one could also round down continuous covariates, effectively reducing the number of split points to search over.

### The variance estimates are jumpy or very large.

In this case, it would be good to try growing a larger number of trees. Obtaining good variance estimates often requires growing more trees than it takes to only obtain accurate predictions. See the discussion under [`num.trees`](#numtrees) and [Merging Forests](#merging-forests).

### The causal forest method is producing nonsensical results.

#### Overlap issues

If the output of the `causal_forest` method doesn't pass a sanity check based on your knowledge of the data, it may be worth checking whether the overlap assumption is violated. In order for conditional average treatment effects to be properly identified, a dataset's propensity scores must be bounded away from 0 and 1. A simple way to validate this assumption is to calculate the propensity scores by regressing the treatment assignments W against X, and examining the out-of-bag predictions. Concretely, you can perform the following steps:

```
propensity.forest = regression_forest(X, W)
W.hat = predict(propensity.forest)$predictions
hist(W.hat, xlab = "propensity score")
```

If there is strong overlap, the histogram will be concentrated away from 0 and 1. If the data is instead concentrated at the extremes, the overlap assumption likely does not hold.

For further discussion of the overlap assumption, please see Imbens and Rubin (2015). In practice, this assumption is often violated due to incorrect modeling decision: for example one covariate may be a deterministic indicator that the example received treatment.

#### Subgroups based on out-of-bag (OOB) estimates

When forming subgroups based on which quantile of the OOB treatment effect estimates a unit belongs to, it is possible to observe the counterintuitive result that the average treatment effect in the high group is small, and the average treatment effect in the low group is large.

This is unfortunately a known artifact of using OOB estimates for this kind of exercise and can be avoided by doing for example a train/test split where one sample is used to estimate subgroups and another to estimate average treatment effects. See for example the documentation for the function `rank_average_treatment_effect`.

### Regression forest predictions differ from those of the randomForest and ranger packages.

While the algorithm in `regression_forest` is very similar to that of classic random forests, it has several notable differences, including 'honesty', group tree training for variance estimates, and restrictions during splitting to avoid imbalanced child nodes. These features can cause the predictions of the algorithm to be different, and also lead to a slower training procedure than other packages. We welcome GitHub issues that shows cases where GRF does notably worse than other packages (either in statistical or computational performance), as this will help us choose better defaults for the algorithm, or potentially point to a bug.


### Forests predict different values depending on the platform even though the seed is the same

Overall, GRF is designed to produce the same estimates across platforms when using a consistent value for the random seed through the training option seed. However, there are still some cases where GRF can produce different estimates across platforms. When it comes to cross-platform predictions, the output of GRF will depend on a few factors beyond the forest seed.

One such factor is the compiler that was used to build GRF. Different compilers may have different default behavior around floating-point rounding, and these could lead to slightly different forest splits if the data requires numerical precision. Another factor is how the forest construction is distributed across different threads. Right now, our forest splitting algorithm can give different results depending on the number of threads that were used to build the forest.

Therefore, in order to ensure consistent results, we provide the following recommendations.
- Make sure arguments `seed` and `num.threads` are the same across platforms
- Round data to 8 significant digits

Also, please note that we have not done extensive testing on Windows platforms, although we do not expect random number generation issues there to be different from Linux/Mac. Regardless of the platform, if results are still not consistent please help us by submitting a Github issue.


## References

Athey, Susan, Julie Tibshirani and Stefan Wager. Generalized Random Forests, *Annals of Statistics*, 2019. ([arxiv](https://arxiv.org/abs/1610.01271))

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, and James Robins. Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 2018.

Cui, Yifan, Michael R. Kosorok, Erik Sverdrup, Stefan Wager, and Ruoqing Zhu. Estimating Heterogeneous Treatment Effects with Right-Censored Data via Causal Survival Forests. *Journal of the Royal Statistical Society: Series B*, 2023.

Ghosal, Indrayudh, and Giles Hooker. Boosting Random Forests to Reduce Bias; One-Step Boosted Forest and its Variance Estimate. *Journal of Computational and Graphical Statistics*, 2020.

Ishwaran, Hemant, Udaya B. Kogalur, Eugene H. Blackstone, and Michael S. Lauer. Random survival forests. *The Annals of Applied Statistics*, 2008.

Imbens, Guido W., and Donald B. Rubin. Causal inference in statistics, social, and biomedical sciences. *Cambridge University Press*, 2015.

Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky. Balancing covariates via propensity score weighting. *Journal of the American Statistical Association*, 2018.

Mayer, Imke, Erik Sverdrup, Tobias Gauss, Jean-Denis Moyer, Stefan Wager and Julie Josse. Doubly robust treatment effect estimation with missing attributes. *The Annals of Applied Statistics*, 2020.

Nie, Xinkun, and Stefan Wager. Quasi-oracle estimation of heterogeneous treatment effects. *Biometrika*, 2021.

Robins, James M., Andrea Rotnitzky, and Lue Ping Zhao. Estimation of regression coefficients when some regressors are not always observed. *Journal of the American statistical Association*, 1994.

Robinson, Peter M. Root-n-consistent semiparametric regression. *Econometrica*, 1988.

Twala, B. E. T. H., M. C. Jones, and David J. Hand. Good methods for coping with missing data in decision trees. *Pattern Recognition Letters 29*,2008.

Van Der Laan, Mark J., and Daniel Rubin. Targeted maximum likelihood learning. *The International Journal of Biostatistics*, 2006.

Wager, Stefan, and Susan Athey. Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 2018.

Yadlowsky, Steve, Scott Fleming, Nigam Shah, Emma Brunskill, and Stefan Wager. Evaluating Treatment Prioritization Rules via Rank-Weighted Average Treatment Effects. *arXiv preprint arXiv:2111.07966*, 2021.
