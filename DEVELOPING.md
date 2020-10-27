# Developing

In addition to providing out-of-the-box forests for quantile regression and instrumental variables, grf provides a framework for creating forests tailored to new statistical tasks. Certain components around splitting and prediction can be swapped out, within the general infrastructure for growing and predicting on trees.

## Contributing

This repository follows the standard open source protocol and setup with git where there is an abundance of existing resources to get up to speed (see for example the
contributing guidelines for well known packages in other languages, like Scikit-learn, Scipy, and pandas)

Condensed greatly, the workflow is to fork this repository, check out a branch, commit your changes (forming an ideally legible commit history),
then submitting a pull request explaining your contribution, ideally referring to the issue you created, or the issue you chose to work on.

## Working with the code

The core forest implementation is written in C++, with an R interface powered by Rcpp. We recommend using a full-powered C++ IDE such as CLion, Xcode, or Visual Studio when working with the core code.

## R package

To build the R package from source, cd into `r-package` and run `build_package.R`. Required development dependencies are listed there. This mimics the tests run when submitting a pull request. Additional online package documentation is built on Travis with [pkgdown](https://pkgdown.r-lib.org/). Usage examples in the form of R Markdown files under _grf\vignettes_ are built and rendered and the R method reference (along with the package articles) is displayed according to the layout defined in _pkgdown.yml_. To build the site locally run `pkgdown::build_site()` from the R package directory.

An alternative development workflow is to use the accompanying grf.Rproj and build and test the package with RStudio's build menu, which can be convenient
for quickly iterating C++/R code changes.

### Note for Windows users:

Symlinks in the src directory point to the core C++ and R bindings. On Windows one has to clone this repository with symlinks enabled: `git clone -c core.symlinks=true https://github.com/grf-labs/grf.git` (this command needs to be run as an administrator: right click _Command Prompt -> Run as administrator_). Caveat: the above RStudio workflow is not tested on Windows.

## Core C++

### Code structure

![GRF Architecture Diagram](https://raw.githubusercontent.com/grf-labs/grf/master/images/arch_diagram.png)

The forest implementation is composed of two top-level components, [ForestTrainer](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestTrainer.h) and [ForestPredictor](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestPredictor.h).

ForestTrainer drives the tree-growing process, and has two pluggable components.
* [RelabelingStrategy](https://github.com/grf-labs/grf/blob/master/core/src/relabeling/RelabelingStrategy.h) is applied before every split, and produces a set of relabelled outcomes given the observations for a group of samples. In the case of quantile forests, for example, this strategy computes the quantiles for the group of samples, then relabels them with a factor corresponding to the quantile they belong to.
* [SplittingRule](https://github.com/grf-labs/grf/blob/master/core/src/splitting/SplittingRule.h) is called to find the best split for a particular node, given a set of outcomes. There are currently implementations for standard regression and multinomial splitting.

The trained forest produces a [Forest](https://github.com/grf-labs/grf/blob/master/core/src/forest/Forest.h) object. This can then be passed to the ForestPredictor to predict on test samples. The predictor has a pluggable 'prediction strategy', which computes a prediction given a test sample. Prediction strategies can be one of two types:
* [DefaultPredictionStrategy](https://github.com/grf-labs/grf/blob/master/core/src/prediction/DefaultPredictionStrategy.h) computes a prediction given a weighted list of training sample IDs that share a leaf with the test sample. Taking quantile forests as an example, this strategy would compute the quantiles of the weighted leaf samples.
* [OptimizedPredictionStrategy](https://github.com/grf-labs/grf/blob/master/core/src/prediction/OptimizedPredictionStrategy.h) does not predict using a list of neighboring samples and weights, but instead precomputes summary values for each leaf during training, and uses these during prediction. This type of strategy will also be passed to ForestTrainer, so it can define how the summary values are computed.

Prediction strategies can also compute variance estimates for the predictions, given a forest trained with grouped trees. Because of performance constraints, only 'optimized' prediction strategies can provide variance estimates.

A particular type of forest is created by pulling together a set of pluggable components. As an example, a quantile forest is composed of a QuantileRelabelingStrategy, ProbabilitySplittingRule, and QuantilePredictionStrategy.
The factory classes [ForestTrainers](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestTrainers.h) and [ForestPredictors](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestPredictors.h) define the common types of forests like regression, quantile, and causal forests.

### Creating a custom forest

When designing a forest for a new statistical task, we suggest following the format suggested in the paper: provide a custom RelabelingStrategy and PredictionStrategy, but use the standard RegressionSplittingRule. If modifications to the splitting rule must be made, they should be done carefully, as it is an extremely performance-sensitive piece of code.

To avoid the overhead of setting up new classes and Rcpp bindings, we provide a template for a 'custom' forest. To get started, fill in the implementation for CustomRelabelingStrategy and CustomPredictionStrategy. By default, the custom forest uses the standard regression splitting rule -- if you'd like to change this, consult the ForestTrainers.custom_trainer method.

This forest is made available in R as `custom.forest` (with `predict.custom.forest`). You can find a starter template for a test exercising the forest in `testthat/test_custom_forest.R`. Note that you'll need to re-run `build_package.R` after making changes to the C++ source.

### Tree Splitting Algorithm

The follow section outlines pseudocode for some of the components listed above.

Each splitting rule follow the same pattern: given a set of samples, and a possible split variable, iterate over all the split points for that variable and compute an error criterion on both sides of the split, then select the split where the decrease in the summed error criterion is smallest. For `RegressionSplittingRule` the criterion is the mean squared error, and is usually referred to as standard CART splitting. `InstrumentalSplittingRule` use the same error criterion, but it incorporates more constraints on potential splits, such as requiring a certain number of treated and control observations (more details are in the [Algorithm Reference](https://grf-labs.github.io/grf/REFERENCE.html)).

Two added features to grf's splitting rule beyond basic CART is the addition of sample weights and support for missing values with Missingness Incorporated in Attributes (MIA) splitting. The algebra behind the decrease calculation is detailed below.

***Algorithm*** (`RegressionSplittingRule`): find the best split for variable `var` at n samples `samples = i...j`

_Input_: `X`: covariate matrix, `weights`: weight n-vector, `response`: the responses (pseudo-outcomes) n-vector, `min_child_size`: split constraint on number of samples

_Output_: The best split value and the direction to send missing values

```
let x = X[samples, var]

for each unique value u in x:
  let count_left[u] be the number of samples with x <= u
  let sum_left[u], weight_sum_left[u] be the sum of weight * response and weight with x <= u

  let count_right[u], sum_right[u], weight_sum_right[u] be the same for samples with x > u
  let count_missing, sum_missing, weight_sum_missing be the same for samples with x = missing

  decrease[u, missing = on_left] = (sum_left[u] + sum_missing)^2 / (weight_sum_left[u] + weight_sum_missing)
                                  + sum_right[u]^2 / weight_sums_right[u]

  decrease[u, missing = on_right] = sum_left[u]^2 / weight_sum_left[u]
                                  + (sum_right[u] + sum_missing)^2 / (weight_sums_right[u] + weight_sum_missing)

return (u, missing) that maximizes decrease such that                                 
count_left[u] and count_right[u] >= min_child_size
```

---

***Algorithm*** (Decrease calculation in weighted CART): find the best split `Sj` for variable `Xj`. The impurity criterion is the sum of weighted mean squared errors, which for a specific split `Sj` is minimized by the weighted sample average `Ybar` (see for example [Elements of Statistical Learning, Ch 9](https://web.stanford.edu/~hastie/ElemStatLearn/)).

_Input_: `Xj` : Covariate j, `w`: sample weight vector, `y`: response vector

_Output_: The best split criterion

```
minimize wrt. s: MSE(s, left) + MSE(s, right)

where

MSE(s, left) = \sum_{i \in L(s)} wi [yi - Ybar_L(s)]^2
MSE[s, right] = \sum_{i \in R(s)} wi [yi - Ybar_R(s)]^2

L(s) are all points Xj such that Xj <= s, and R(s) are all points such that Xj > s.
Ybar_L(s) and Ybar_R(s) are the weighted sample averages for the respective partitions.

Expand the child impurity and write out the averages (minimize wrt. s):
= \sum_{i \in L(s)} wi yi^2 + W_L(s) [Ybar_L(s)]^2
 + \sum_{i \in R(s)} wi yi^2 + W_R(s) [Ybar_R(s)]^2
= \sum_i wi yi^2
  + 1/W_L(s) [\sum_{i \in L(s)} wi yi]^2
  + 1/W_R(s) [\sum_{i \in R(s)} wi yi]^2

Where W is the sum of sample weights.

The parent impurity is (minimize wrt. s):
= \sum_i wi [yi - Ybar]^2
= \sum_i wi yi^2 + W Ybar^2
= \sum_i wi yi^2 + 1/W [\sum_i wi yi]^2

Child impurity <= parent impurity then reduces to

1/W_L(s) [\sum_{i \in L(s)} wi yi]^2 + 1/W_R(s) [\sum_{i \in R(s)} wi yi]^2 <= 1/W [\sum_i wi yi]^2

Or sum_left^2 / weight_sum_left + sum_right^2 / weight_sum_right
```

For multivariate CART `sum_left^2` and `sum_right^2` becomes the squared L2 norm.

---

***Algorithm*** (`SurvivalSplittingRule`): find the best split for variable `var` at n samples `samples = i...j`.

_Input_: `X`: covariate matrix, `response`: the failure times, `min_child_size`: split constraint on number of failures, `m`: the number of failure times in this node.

_Output_: The best split value and the direction to send missing values

grf performs survival splits by finding the partition that maximizes the [log-rank](https://en.wikipedia.org/wiki/Logrank_test) test for the comparison of the survival distribution in the left and right node. Using the notation from [Ishwaran et. al (2008), equation 1](https://kogalur.github.io/randomForestSRC/theory.html) the statistic is:

```
logrank(x, split.value) =  sum over all times k up to m
  (dk,l - Yk,l * dk/Yk)  /
  (Yk,l / Yk * (1 - Yk,l / Yk) * (Yk - dk) / (Yk - 1) dk)

dk,l: the count of failures in the left node at time k
Yk,l: the count of observations at risk in the left node at time k (At risk: all i such that Ti >= k)
dk: the count of failures in the parent node at time k
Yk: the count of observations at risk in the parent node at time k
```

Note that the absolute value of the response variable (the failure time) is never used, only counts matter, constructed from the relative ordering between them. Therefore, the response variable is always relabeled to range from consecutive integers starting at zero and ending at the number of failures. This can be done quickly with binary search:

```
let failure_times be a sorted vector of the original failure values (T1, T2, .. TM)
initialize response_relabeled = zeros(n)
for each i=1:n
  relabeled_failure_time = binary_search(response[i], failure_times)
  response_relabeled[i] = relabeled_failure_time
```

A label of 0 is given to all samples with response less than the first failure time, a value of 1 is given to all samples with failure value greater or equal to the first failure value and less than the second failure value, etc.

The algorithm proceeds in the same manner as outlined above for RegressionSplitting: iterate over all possible splits and calculate the logrank statistic, one with all missing values on the left, and one with all missing on the right. Then select the split that yielded the maximum logrank test, subject to a constraint that there are a sufficient amount of failures on both sides of the split.
