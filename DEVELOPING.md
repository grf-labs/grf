# Developing

## Working with the code

The core forest implementation is written in C++, with an R interface powered by Rcpp. We recommend using a full-powered C++ IDE such as CLion, Xcode, or Visual Studio when working with the core code. To build the R package from source, cd into `r-package` and run the `build_package.R` script.

## Code structure

The forest implementation is composed of two top-level components, [ForestTrainer](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestTrainer.h) and [ForestPredictor](https://github.com/swager/gradient-forest/blob/master/core/src/forest/ForestPredictor.h).

ForestTrainer drives the tree-growing process, and has two pluggable components.
* [RelabelingStrategy](https://github.com/swager/gradient-forest/blob/master/core/src/relabeling/RelabelingStrategy.h) is applied before every split, and produces a set of relabelled outcomes given the observations for a group of samples. In the case of quantile forests, this strategy computes the quantiles for the group of samples, then relabels them with a factor corresponding to the quantile they belong to.
* [SplittingRule](https://github.com/swager/gradient-forest/blob/master/core/src/splitting/SplittingRule.h) is called to find the best split for a particular node. The existing implementations perform standard regression and multinomial splitting.

The trained forest produces a [Forest](https://github.com/swager/gradient-forest/blob/master/core/src/forest/Forest.h) object. This can then be passed to the ForestPredictor to predict on test samples.

## Creating a custom forest

If no relabeling is needed (as in the case of a vanilla regression forest), the NoopRelabelingStrategy should be used.

