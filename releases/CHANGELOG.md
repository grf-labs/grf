# Changelog
All notable changes to `grf` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.3.0] - 2023-05-10

### Added
- Add support for a continuous treatment `W` in `causal_survival_forest`, mimicking the average partial effect estimand in `causal_forest`. [#1280](https://github.com/grf-labs/grf/pull/1280)
- Add `target.sample = "overlap"` option to `best_linear_projection` allowing for estimation with weights equal to e(X)(1 - e(X)). [#1258](https://github.com/grf-labs/grf/pull/1258)
- Add optional `drop` argument to `predict.multi_arm_causal_forest`, `predict.multi_regression_forest`, and `predict.lm_forest` which drops the singleton dimension in the prediction array in case the forest is fit with only a single outcome. [#1271](https://github.com/grf-labs/grf/pull/1271), [#1281](https://github.com/grf-labs/grf/pull/1282)
- Add column names to the `variance.estimates` matrix when calling `predict.multi_arm_causal_forest`, `predict.probability_forest`, and `predict.lm_forest`. [#1272](https://github.com/grf-labs/grf/pull/1272), [#1283](https://github.com/grf-labs/grf/pull/1283)

### Fixed
- Fix missing `num.threads` argument in internal `causal_survival_forest` nuisance estimation that would cause some components to utilize all available threads even though the `num.threads` arguments specified less. [#1267](https://github.com/grf-labs/grf/pull/1267)
- Drop C++11 compiler flag in src/Makevars per latest CRAN guidelines. [#1305](https://github.com/grf-labs/grf/pull/1305)

## [2.2.1] - 2022-12-14

### Added
- Add `rank_average_treatment_effect.fit`, an optional interface to RATE which allows for user-supplied doubly robust evaluation set scores. [#1242](https://github.com/grf-labs/grf/pull/1242)

### Fixed
- Fix minor docstring example typo in `rank_average_treatment_effect`. [#1209](https://github.com/grf-labs/grf/pull/1209)
- Move a missing values input check in `rank_average_treatment_effect`. [#1212](https://github.com/grf-labs/grf/pull/1212)
- Allow factor variables in the covariate matrix `A` passed to `best_linear_projection`. [#1215](https://github.com/grf-labs/grf/pull/1215)
- Fix the "ai1" and "ai2" example DGPs included in `dgps.R`. [#1241](https://github.com/grf-labs/grf/pull/1241)

## [2.2.0] - 2022-08-06

### Added
- Add new "linear model" forest `lm_forest` allowing for forest-based estimation of the conditionally linear model Y = c(x) + h_1(x)W_1 + ... + h_K(x)W_K at X = x (this generalization nests causal and multi-arm causal forest). As an example, one could use this forest to estimate heterogeneous treatment effects in a regression discontinuity design by setting W_1 to 1{Z >= cutoff} and W_2 = Z, where Z is a running variable and "cutoff" a treatment threshold. [#1138](https://github.com/grf-labs/grf/pull/1138)
- Add optional `prediction.times` argument to `predict.survival_forest`. The new option "times" allows for easier computation of IPCW-type weights. [#1139](https://github.com/grf-labs/grf/pull/1139)
- Add support for `instrumental_forest` to `best_linear_projection`. [#1158](https://github.com/grf-labs/grf/pull/1158)

### Fixed
- Fix optional `failure.times` argument in `predict.survival_forest`. This only affects `survival_forest` predictions in cases the optional `failure.times` contains time points outside the training grid. [#1132](https://github.com/grf-labs/grf/pull/1132)
- Fix `ll_regression_forest` such that all allowable input types `X` work when `enable.ll.split = TRUE`. [#1198](https://github.com/grf-labs/grf/pull/1198)

## [2.1.0] - 2022-03-17

### Changed (breaking)
**IMPORTANT** Some of these changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- The optional event grid (`failure.times`) in `survival_forest` is now required to be strictly increasing. [#1034](https://github.com/grf-labs/grf/pull/1034)
- The time values `Y` in `survival_forest` are now required to be non-negative. [#1058](https://github.com/grf-labs/grf/pull/1058)

### Added
- Add new feature `rank_average_treatment_effect` ("RATE") which calculates a scalar metric together with a Targeting Operating Characteristic curve which can be used to assess how well a CATE estimator does in identifying subpopulations which benefit from treatment. [#1086](https://github.com/grf-labs/grf/pull/1086)
- Add support for survival probability difference estimation with `causal_survival_forest` using new required arguments `target` and `horizon`. **NOTE**: this change breaks the experimental API from the previous release and marks a new stable interface. [#1055](https://github.com/grf-labs/grf/pull/1055)
- Add warning if `multi_arm_causal_forest` is trained with a treatment factor vector `W` which contains missing levels. [#1039](https://github.com/grf-labs/grf/pull/1039)
- Add documentation to `multi_arm_causal_forest` emphasizing that in case multiple outcomes `Y` are supplied they should be on the same scale. [#1041](https://github.com/grf-labs/grf/pull/1041)
- Add documentation note to `survival_forest` emphasizing that the `alpha` parameter works as a split constraint on the number of events in child nodes, thus suggesting that if the forest does not split on very low event-rate data lowering this may be useful. [#1072](https://github.com/grf-labs/grf/pull/1072)

### Fixed
- Fix doubly robust scores for `multi_arm_causal_forest`. For forests trained with more than two treatments the doubly robust score construction used for ATE estimation had a typo which could lead to slight under-coverage. [#1021](https://github.com/grf-labs/grf/pull/1021)
- Allow forest tuning to continue with tuning even when some random parameter draws are inadmissible. This prevents the message "Could not tune forest because some small forest error estimates were NA." from appearing in some cases. [#1067](https://github.com/grf-labs/grf/pull/1067)
- For forests requiring estimates of additional nuisance components in summary functions (such as a causal forest with continuous treatment), use the original forest seed when estimating these with a new nuisance forest, avoiding confusion for users who expect the ATE summary functions to produce the exact same result across different invocations. [#1070](https://github.com/grf-labs/grf/pull/1070)
- Fix sample weighted `average_treatment_effect` for `target.sample = c("treated", "control")` in the case the forest is trained without clusters. [#1103](https://github.com/grf-labs/grf/pull/1103)
- Disable sample weighted `average_treatment_effect` for `method = "TMLE"` as these were silently ignored anyways. [#1102](https://github.com/grf-labs/grf/pull/1102)
- Fix sample weighted `average_treatment_effect` for the rare case when some weights are zero. [#1104](https://github.com/grf-labs/grf/pull/1104), [#1105](https://github.com/grf-labs/grf/pull/1105)
- Remove undocumented legacy "sentinel" `seed = 0` which meant a new random seed was drawn every time. [#1093](https://github.com/grf-labs/grf/pull/1093)
- Fix documentation notation in `survival_forest` where `Y` appeared instead of `T`. [#1056](https://github.com/grf-labs/grf/pull/1056)

## [2.0.2] - 2021-07-14

### Fixed
- Minor patch release for CRAN Solaris compatibility. [#1011](https://github.com/grf-labs/grf/pull/1011)

## [2.0.0] - 2021-06-22

### Changed (breaking)
**IMPORTANT** Some of these changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- Unify the interface to ATE-type estimators: 1) `average_treatment_effect` is the new entry point for all ATE summaries, meaning `average_late` and `average_partial_effect` is removed. 2) This function now targets population-type quantities for all forests, meaning some confidence intervals may be slightly wider than before. 3) Some ad-hoc normalization schemes are removed, but can be manually specified through the `debiasing.weights` argument. [#723](https://github.com/grf-labs/grf/pull/723)
- Remove all `tune_***_forest` functions, restricting the tuning interface to the pre-existing `tune.parameters` argument in all tuning-compatible forests. [#790](https://github.com/grf-labs/grf/pull/790)
- Remove the optional `orthog.boosting` argument in `causal_forest`, since tailored `m(x)` estimates can be passed through the existing `Y.hat` argument. [#892](https://github.com/grf-labs/grf/pull/892)
- Remove support for sparse `X` in forest training (as the internal C++ implementation did not leverage sparsity in `X` beyond storage mode). To train a forest with sparse data do `forest(as.matrix(X), Y)`. [#939](https://github.com/grf-labs/grf/pull/939)
- Remove `custom_forest`. For a template for getting started with a custom GRF estimator, consider using an existing simple forest, like `regression_forest` as a scaffold. For more details see the GRF [developing document](https://grf-labs.github.io/grf/DEVELOPING.html). [#870](https://github.com/grf-labs/grf/pull/870)
- Rename `get_sample_weights` to `get_forest_weights`. [#894](https://github.com/grf-labs/grf/pull/894)
- Return `quantile_forest` predictions in the `predictions` attribute of a new output list in order to conform with the GRF convention of returning point predictions as `predict(forest)$predictions`. [#822](https://github.com/grf-labs/grf/pull/822)
- Change the way optional sample weights (passed through the `sample.weights` argument) interact with the GRF forest weights. When forming estimates according to (2) and (3) in the [GRF paper](https://arxiv.org/pdf/1610.01271.pdf), sample weights now enter through alpha_i(x)' = alpha_i(x) * sample.weight_i. In addition, `causal_forest` and `instrumental_forest` now take sample weights into account in the relabeling step ([#752](https://github.com/grf-labs/grf/pull/752)). Sample weights are also explicitly disabled for local linear forests ([#841](https://github.com/grf-labs/grf/pull/841)). [#796](https://github.com/grf-labs/grf/issues/796)

### Added
- Add `causal_survival_forest` for estimating conditional average treatment effects with right-censored data. [#660](https://github.com/grf-labs/grf/pull/660)
- Add `multi_arm_causal_forest`, an extension of `causal_forest` to multiple categorial treatments `W`, and optionally multiple responses `Y`. [#748](https://github.com/grf-labs/grf/pull/748)
- Add `multi_regression_forest` for estimating several conditional mean functions mu_i(x) = E[Y_i | X = x]. [#742](https://github.com/grf-labs/grf/pull/742)
- Add `probability_forest` for estimating conditional class probabilities P[Y = k | X = x]. [#711](https://github.com/grf-labs/grf/pull/711)
- Add `get_scores` returning doubly robust scores for a number of estimands. [#732](https://github.com/grf-labs/grf/pull/732)
- Add `get_leaf_node` utility function which given a GRF tree object returns the leaf node a test sample falls into. [#739](https://github.com/grf-labs/grf/pull/739)
- Add a `vcov.type` standard error option to `test_calibration` and `best_linear_projection`. On large datasets with clusters, setting this option to `"HC0"` or `"HC1"` will significantly speed up the computation.
- Add optional Nelson-Aalen estimates of the survival function. [#685](https://github.com/grf-labs/grf/pull/685)
- Add a docstring example to `survival_forest` on how to calculate concordance with the optional `survival` package. [#956](https://github.com/grf-labs/grf/pull/956)

### Fixed
- Fix the output name in `average_treatment_effect` when `method = "TMLE"`. [#864](https://github.com/grf-labs/grf/pull/864)
- Fix pointwise variance estimates in (the very unlikely) zero variance case. [#907](https://github.com/grf-labs/grf/pull/907)
- Fix `survival_forest` test set predictions with sample weights. [#969](https://github.com/grf-labs/grf/pull/969)
- Make forest tuning respect the `seed` argument when drawing a random grid of parameter values, allowing reproducibility without an explicit `set.seed` before training. [#704](https://github.com/grf-labs/grf/pull/704)

## [1.2.0] - 2020-06-04

### Added
- Add Survival Forest functionality. [#647](https://github.com/grf-labs/grf/pull/647)
- Add optional argument `debiasing.weights` to `average_partial_effect`. [#637](https://github.com/grf-labs/grf/pull/637)
- Add optional `compute.oob.predictions` argument to Quantile Forest. [#665](https://github.com/grf-labs/grf/pull/665)

### Fixed
- Fix a performance regression in DefaultPredictionCollector. This improves prediction speed for forests such as Quantile Forest. [#650](https://github.com/grf-labs/grf/pull/650)
- Predict with training quantiles by default. [#668](https://github.com/grf-labs/grf/pull/668)

## [1.1.0] - 2020-03-12

### Changed (breaking)
**IMPORTANT** These changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- Performance improvement: remove an unnecessary splitting rule loop. Note: this may cause very small differences from earlier versions because it changes the order in which potential splits are evaluated. [#592](https://github.com/grf-labs/grf/pull/592)

### Added
- Add support for missing values in the covariates X with [MIA](https://github.com/grf-labs/grf/issues/457) splitting. [#612](https://github.com/grf-labs/grf/pull/612)
- Add local linear splitting. An experimental option `enable.ll.split` fits a forest with splits based on ridge residuals as opposed to standard CART splits. Note: local linear tuning does not take the new splits into account. [#603](https://github.com/grf-labs/grf/pull/603)
- Add sample weighted splitting. Previously, if a user passed `sample.weights`, they would only be used for prediction. Now they are used in splitting as well. Note: this will make results fitted with sample weights different from previous versions. [#590](https://github.com/grf-labs/grf/pull/590)

### Fixed
- Remove a superfluous predict call in tuning. [#597](https://github.com/grf-labs/grf/pull/597)
- Fix `average_partial_effect` calibration in case of low variation W.hat. [#611](https://github.com/grf-labs/grf/pull/611)
- Update `best_linear_projection` to handle non-binary treatment. [#615](https://github.com/grf-labs/grf/pull/615)
- Add an error message in case summary functions are passed a subset that refers to too few distinct units. [#629](https://github.com/grf-labs/grf/pull/629)

## [1.0.1] - 2019-12-05

### Fixed
- Fix a bug where the nodes of the printed trees would be in the wrong order. [#587](https://github.com/grf-labs/grf/pull/587)

## [1.0.0] - 2019-11-29

### Changed (breaking)
**IMPORTANT** These changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- Rename `prune.empty.leaves` to `honesty.prune.leaves`. [#529](https://github.com/grf-labs/grf/pull/529)
- Simplify the forest tuning API. Previously, tuning was enabled during forest training by setting the option `tune.parameters=TRUE`. All relevant parameters were tuned by default, except for those that were explicitly passed to the forest (like `min.node.size=100`). Now the option `tune.parameters` directly takes a list of parameters to tune, for example `tune.parameters=c("min.node.size", "mtry")`, or `tune.parameters="all"`. [#534](https://github.com/grf-labs/grf/pull/534)
- Change how data points are weighted in cluster-robust estimation. Previously, each cluster was given equal weight when training the forest and computing estimates. Now, each point is weighted equally regardless of its cluster size. This behavior can be controlled through a new option `equalize.cluster.weights`, which defaults to `FALSE` but can be set to `TRUE` to match the old behavior of weighting clusters equally. The old option `samples.per.cluster` has been removed. [#545](https://github.com/grf-labs/grf/pull/545).

### Added
- Improve the performance of `get_tree`. [#528](https://github.com/grf-labs/grf/pull/528)
- Add support for tuning instrumental forests (currently marked 'experimental'). [#547](https://github.com/grf-labs/grf/pull/547)
- Introduce optimizations to tree splitting. These improvements lead to a small speed-up in forest training. [#560](https://github.com/grf-labs/grf/pull/560), [#561](https://github.com/grf-labs/grf/pull/561)
- Add `best_linear_projection`, a doubly robust estimate of the best linear projection of the conditional average treatment effect onto a set of covariates. [#574](https://github.com/grf-labs/grf/pull/574)
- Speed up forest prediction by introducing additional parallelization. [#566](https://github.com/grf-labs/grf/pull/566), [#576](https://github.com/grf-labs/grf/pull/576)

### Fixed
- Allow the data matrix `X` to be a data frame. [#540](https://github.com/grf-labs/grf/pull/540)
- When merging forests, validate that all forests were trained on the same data. [#543](https://github.com/grf-labs/grf/pull/543)
- Fix a major performance issue in `get_sample_weights`. [#578](https://github.com/grf-labs/grf/pull/578)

## [0.10.4] - 2019-09-01

### Changed (breaking)
**IMPORTANT** These changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- Ensure forest estimates are consistent across platforms. [#469](https://github.com/grf-labs/grf/pull/469), [#492](https://github.com/grf-labs/grf/pull/492)
- The number of trees used for orthogonalization was changed from `min(500, num.trees)` to `max(50, num.trees / 4)`. [#439](https://github.com/grf-labs/grf/pull/439)
- Solidify the parameter tuning procedure. If the optimization procedure fails, or if the selected parameters perform worse than defaults, we now return default parameters instead. [#455](https://github.com/grf-labs/grf/pull/455)
- Introduce parameters `honesty.fraction` and `prune.empty.leaves` to help mitigate the effect of honesty on small datasets, and tune over them when `tune.parameters=TRUE`. [#456](https://github.com/grf-labs/grf/pull/456), [#484](https://github.com/grf-labs/grf/pull/484)

### Added
- Add variance estimates for local linear forests. [#442](https://github.com/grf-labs/grf/pull/442)
- Include information about leaf samples in plotting and printing. [#460](https://github.com/grf-labs/grf/pull/460)
- Add example of saving a plot with DiagrammeRsvg. [#478](https://github.com/grf-labs/grf/pull/478)
- Support average effect estimates for instrumental forests (ACLATE). [#490](https://github.com/grf-labs/grf/pull/490)

### Fixed
- Performance improvements to forest training. [#514](https://github.com/grf-labs/grf/pull/514)

## [0.10.3] - 2019-06-01

### Changed (breaking)
**IMPORTANT** These changes might cause small differences in results compared to previous releases, even if the same random seed is used.
- Fix two bugs in the termination criterion for tree splitting.
  - Remove the purity condition on outcomes during splitting. For all tree types, we used to stop splitting if all outcomes in a leaf are the same. This behavior does not make sense for causal forests (which incorporates other observations besides the outcome), so it was removed. [#362](https://github.com/grf-labs/grf/pull/362)
  - Stop splitting if the objective can no longer be improved. With this change, `causal_forest` may split slightly less aggressively. [#415](https://github.com/grf-labs/grf/pull/415)

### Added
- In out-of-bag prediction, return the Monte Carlo error alongside the debiased error. [#327](https://github.com/grf-labs/grf/pull/327)
- Allow for passing a factor for the `cluster` parameter. [#329](https://github.com/grf-labs/grf/pull/329)
- Support taking a union of forests through the `merge_forests` method. [#347](https://github.com/grf-labs/grf/pull/347)
- Include a summary of the parameter tuning procedure in the forest object. [#419](https://github.com/grf-labs/grf/pull/419)
- Add experimental support for sample weighting to regression, causal, and instrumental forests. [#376](https://github.com/grf-labs/grf/pull/376), [#418](https://github.com/grf-labs/grf/pull/418)
- Add an experimental new forest type `boosted_regression_forest`, which applies boosting to regression forests. Allow boosting to be used during orthogonalization through the `orthog.boosting` parameter. [#388](https://github.com/grf-labs/grf/pull/388)

### Fixed
- Improve input data validation. [#354](https://github.com/grf-labs/grf/pull/354), [#378](https://github.com/grf-labs/grf/pull/378), [#430](https://github.com/grf-labs/grf/pull/430)
- Improve the `test_calibration` function by switching to one-sided p-values. [#370](https://github.com/grf-labs/grf/pull/370)
- For custom forests, fix a bug in OOB prediction where the train and tests datasets were switched. [#372](https://github.com/grf-labs/grf/pull/372)
- Decrease memory usage during training and out-of-bag prediction. [#408](https://github.com/grf-labs/grf/pull/408), [#412](https://github.com/grf-labs/grf/pull/412)
- Allow roxygen to autogenerate the `NAMESPACE` file. [#423](https://github.com/grf-labs/grf/pull/423), [#428](https://github.com/grf-labs/grf/pull/428)

## [0.10.2] - 2018-11-23
### Added
- Add support for confidence intervals in local linear regression forests.

### Changed
- Allow samples_per_cluster to be larger than smallest cluster size.

### Fixed
- Make sure average effect estimation doesn't error on data with a single feature.
- Fix a bug in local linear prediction where the penalty wasn't properly calculated.
- Fix two issues in causal forest tuning that could lead to unstable results.
- Ensure that the ATE and APE functions correctly account for cluster membership.

## [0.10.1] - 2018-09-23
### Added
- Add basic support for tree plotting (through the `plot` method).
- Add the method `test_calibration`, which performs an omnibus test for presence of heterogeneity via calibration.
- For local linear regression forests, add support for selecting the value of `ll.lambda` through cross-validation.
- Introduce a training option `honesty.fraction` that can be used to specify the fraction of data that should be used in selecting splits vs. performing estimation. Note that this parameter is only relevant when honesty is enabled (the default).
- Start a practical guide to the `grf` algorithm (https://github.com/grf-labs/grf/blob/master/REFERENCE.md).
- In `average_treatment_effect` and `average_partial_effect`, add an option `subset` to support estimating the treatment effect over a subsample of the data.

### Fixed
- Fix a bug in random sampling where features listed earlier in the data matrix were more likely to be selected for splitting.
- Make sure that the sample indices returned in `get_tree` are 1-indexed, as opposed to 0-indexed.

## [0.10.0] - 2018-05-08
### Added
- Replace the local.linear option with linear.correction.variables, which allows for a subset of
  variables to be considered during local linear regression.
- Update causal_forest interface to allow user to specify Y.hat and W.hat. Note that these options supercede
  the precompute.nuisance parameter, which has been removed. To recreate the behavior of
  precompute.nuisance = TRUE, NULL can be provided for Y.hat and W.hat, and for precompute.nuisance = FALSE,
  Y.hat and W.hat should be 0.
- Add a causal forest example with variable selection and parameter tuning.

### Changed
- Adjust the defaults for the causal forest tuning algorithm.

### Fixed
- Prevent tuning down to a min.node.size of 0.

## [0.9.6] - 2018-04-13

### Added
- Debiased error criterion for measuring the out-of-bag accuracy of a forest using only a few trees.
- Automated tuning via cross-validation for regression and causal forests.
- Estimation of average partial effects with a continuous treatment.
- Overlap-weighted average treatment effects.
- Cluster-robust standard errors for regression and causal forests, and average effect estimates (contributed by @lminer).
- Locally linear prediction in regression forests (contributed by @rinafriedberg).
- Regularize splits in causal/instrumental forests via a variance penalty.

### Changed
- Avoid causal forest leaves with all treated or all control samples (controlled via stabilize.splits = TRUE).
- Store in-bag rather than out-of-bag samples to save memory.
- Only support sampling with replacement (as some features are ambiguously defined with bootstrapping).

### Fixed
- Bug in IV CI construction (#209).

## [0.9.5] - 2018-01-04

### Added
- Create a simple method for variable importance based on split frequency and depth.
- Add support for sparse data matrices of type 'dgCMatrix'.

### Changed
- Use RcppEigen for the R package, as opposed to the eigen source.

### Fixed
- Fix a few places where we were still using the old default for mtry. This issue was
causing poor performance for even moderately large numbers of features.

## [0.9.4] - 2017-11-25

### Changed
- Update the default for mtry to sqrt(p) + 20.

### Fixed
- Fix an issue where split_frequencies fails when p = 1.
- Use a Solaris-compatible version of std::sqrt.

## [0.9.3] - 2017-07-20

### Fixed
- Fix an out of bounds error when there are fewer trees than threads.
- Fix a bug in the get_tree function where the same tree was always returned.

## [0.9.2] - 2017-07-06

### Added
- Add an experimental regularized version of the regression splitting rule.

### Fixed
- Several bugfixes for CRAN compatibility.

## [0.9.1] - 2017-06-27

### Added
- Add a regression.splits option to quantile forests to allow emulating the approach in Meinhausen (2006).

## [0.9.0] - 2017-06-25
First official (beta) release. The package currently supports
- standard regression forests
- causal, instrumental, and quantile forests
- confidence intervals for causal, instrumental, and regression forests
- training 'honest' versions of the above forests
