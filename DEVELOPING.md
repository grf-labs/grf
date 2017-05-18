# Developing

## Working with the code

The core forest implementation is written in C++, with an R interface powered by Rcpp. We recommend using a full-powered C++ IDE such as CLion, Xcode, or Visual Studio when working with the core code. To build the R package from source, cd into `r-package` and run the `build_package.R` script.

## Code structure

The forest implementation is composed of two top-level components, [ForestTrainer](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestTrainer.h) and [ForestPredictor](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestPredictor.h).

ForestTrainer drives the tree-growing process, and has two pluggable components.
* [RelabelingStrategy](https://github.com/swager/gradient-forest/blob/master/core/src/relabeling/RelabelingStrategy.h) is applied before every split, and produces a set of relabelled outcomes given the observations for a group of samples. In the case of quantile forests, this strategy computes the quantiles for the group of samples, then relabels them with a factor corresponding to the quantile they belong to.
* [SplittingRule](https://github.com/swager/gradient-forest/blob/master/core/src/splitting/SplittingRule.h) is called to find the best split for a particular node. There are currently implementations to standard regression and multinomial splitting.

The trained forest produces a [Forest](https://github.com/swager/gradient-forest/blob/master/core/src/forest/Forest.h) object. This can then be passed to the ForestPredictor to predict on test samples. The predictor has a pluggable 'prediction strategy', which computes a prediction given a test sample. Prediction strategies can be one of two types:
* [DefaultPredictionStrategy](https://github.com/swager/gradient-forest/blob/master/core/src/prediction/DefaultPredictionStrategy.h) computes a prediction given a weighted list of training sample IDs that share a leaf with the test sample. For quantile forests, this strategy computes the quantiles of the weighted leaf samples.
* [OptimizedPredictionStrategy](https://github.com/swager/gradient-forest/blob/master/core/src/prediction/OptimizedPredictionStrategy) does not predict using a list of neighboring samples and weights, but instead precomputes summary values for each leaf during training, and uses these during prediction. This type of strategy will also be passed to ForestTrainer, so it can define how the summary values are computed.

Prediction strategies can also compute variance estimates for the predictions, given a forest trained with grouped trees. Because of performance constraints, only 'optimized' prediction strategies can provide variance estimates.

A particular type of forest is created by pulling together a set of pluggable components. As an example, a quantile forest is composed of a QuantileRelabelingStrategy, ProbabilitySplittingRule, and QuantilePredictionStrategy. The factory classes [ForestTrainers](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestTrainers.h) and [ForestPredictors](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestPredictors.h) define the common types of forests like regression, quantile, and causal forests.

## Creating a custom forest

If no relabeling is needed (as in the case of a vanilla regression forest), the NoopRelabelingStrategy should be used.

