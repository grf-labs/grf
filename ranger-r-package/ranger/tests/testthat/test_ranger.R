library(ranger)
library(survival)

context("ranger")

rg.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE)
rg.reg <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE)
rg.surv <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, write.forest = TRUE)

#TODO: This would require the GenABEL package in suggests
# dat.gwaa <- readRDS("../test_gwaa.Rds")
# rg.gwaa <- ranger(CHD ~ ., data = dat.gwaa, verbose = FALSE, write.forest = TRUE)

test_that("classification result is of class ranger with 14 elements", {
  expect_that(rg.class, is_a("ranger"))
  expect_that(length(rg.class), equals(14))
})

test_that("regression result is of class ranger with 14 elements", {
  expect_that(rg.reg, is_a("ranger"))
  expect_that(length(rg.reg), equals(14))
})

test_that("survival result is of class ranger with 15 elements", {
  expect_that(rg.surv, is_a("ranger"))
  expect_that(length(rg.surv), equals(15))
})

#TODO: This would require the GenABEL package in suggests
# test_that("classification gwaa rf is of class ranger with 14 elements", {
#   expect_that(rg.gwaa, is_a("ranger"))
#   expect_that(length(rg.gwaa), equals(14))
# })

test_that("results have 500 trees", {
  expect_that(rg.class$num.trees, equals(500))
  expect_that(rg.reg$num.trees, equals(500))
  expect_that(rg.surv$num.trees, equals(500))
})

test_that("results have right number of independent variables", {
  expect_that(rg.class$num.independent.variables, equals(ncol(iris) - 1))
  expect_that(rg.reg$num.independent.variables, equals(ncol(iris) - 1))
  expect_that(rg.surv$num.independent.variables, equals(ncol(veteran) - 2))
})

test_that("unique death times in survival result is right", {
  expect_that(rg.surv$unique.death.times, equals(sort(unique(veteran$time))))
})

test_that("importance measures work", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "impurity")
  expect_that(rg.imp$variable.importance, is_a("numeric"))
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation")
  expect_that(rg.imp$variable.importance, is_a("numeric"))
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = TRUE)
  expect_that(rg.imp$variable.importance, is_a("numeric"))
})

test_that("gini importance is larger than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "impurity")
  expect_that(rg.imp$variable.importance[1], is_more_than(1))
})

test_that("unscaled importance is smaller than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = FALSE)
  expect_that(rg.imp$variable.importance[1], is_less_than(1))
})

test_that("scaled importance is larger than 1", {
  rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                    importance = "permutation", scale.permutation.importance = TRUE)
  expect_that(rg.imp$variable.importance[1], is_more_than(1))
})

test_that("probability estimations are a matrix with correct size", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]

  rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
  prob <- predict(rg.prob, iris.test)

  expect_that(prob$predictions, is_a("matrix"))
  expect_that(nrow(prob$predictions), equals(nrow(iris.test)))
  expect_that(ncol(prob$predictions), equals(length(rg.prob$forest$levels)))
})

test_that("probability estimations are between 0 and 1 and sum to 1", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]

  rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
  prob <- predict(rg.prob, iris.test)

  expect_that(all(prob$predictions > -1e-5 & prob$predictions <= 1 + 1e-5), is_true())
  expect_that(rowSums(prob$predictions), equals(rep(1, nrow(prob$predictions))))
})

test_that("predict returns good prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE)
  pred <- predict(rf, iris)
  expect_that(mean(iris$Species == predictions(pred)), is_more_than(0.9))
})

test_that("Alternative interface works for classification", {
  rf <- ranger(dependent.variable.name = "Species", data = iris)
  expect_that(rf$treetype, equals("Classification"))
})

test_that("Alternative interface works for regression", {
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = iris)
  expect_that(rf$treetype, equals("Regression"))
})

test_that("Alternative interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran)
  expect_that(rf$treetype, equals("Survival"))
})

test_that("Matrix interface works for classification", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), write.forest = TRUE, classification = TRUE)
  expect_that(rf$treetype, equals("Classification"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(iris)[1:4]))
})

test_that("Matrix interface prediction works for classification", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, classification = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("Matrix interface works for regression", {
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = data.matrix(iris), write.forest = TRUE)
  expect_that(rf$treetype, equals("Regression"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(iris)[2:5]))
})

test_that("Matrix interface prediction works for regression", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Sepal.Length", data = dat, write.forest = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("Matrix interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = data.matrix(veteran), write.forest = TRUE)
  expect_that(rf$treetype, equals("Survival"))
  expect_that(rf$forest$independent.variable.names, equals(colnames(veteran)[c(1:2, 5:8)]))
})

test_that("Matrix interface prediction works for survival", {
  dat <- data.matrix(veteran)
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = dat, write.forest = TRUE)
  expect_that(predict(rf, dat), not(throws_error()))
})

test_that("save.memory option works for classification", {
  rf <- ranger(Species ~ ., data = iris, save.memory = TRUE)
  expect_that(rf$treetype, equals("Classification"))
})

test_that("save.memory option works for regression", {
  rf <- ranger(Sepal.Length ~ ., data = iris, save.memory = TRUE)
  expect_that(rf$treetype, equals("Regression"))
})

test_that("save.memory option works for probability", {
  rf <- ranger(Species ~ ., data = iris, probability = TRUE, save.memory = TRUE)
  expect_that(rf$treetype, equals("Probability estimation"))
})

test_that("C-index splitting works", {
  rf <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, 
               splitrule = "C")
  expect_that(rf$treetype, equals("Survival"))
})

test_that("predict works for single observations, classification", {
  pred <- predict(rg.class, head(iris, 1))
  expect_that(pred$predictions, equals(iris[1,"Species"]))
})

test_that("predict works for single observations, probability prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, head(iris, 1))
  expect_that(names(which.max(pred$predictions)), equals(as.character(iris[1,"Species"])))
})

test_that("predict works for single observations, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, write.forest = TRUE)
  pred <- predict(rf, head(veteran, 1))
  expect_that(length(pred$survival), equals(length(rf$unique.death.times)))
})

test_that("same result with same seed", {
  ind = 1:150 %in% sample(150, 100)
  
  set.seed(2)
  mod1 = ranger(Species ~ ., data = iris[ind, ], write.forest = TRUE, num.trees = 50)
  pred1 = predict(mod1, data = iris[!ind, ])
  
  set.seed(2)
  mod2 = ranger(Species ~ ., data = iris[ind, ], write.forest = TRUE, num.trees = 50)
  pred2 = predict(mod2, data = iris[!ind, ])
  
  expect_that(pred1$predictions, equals(pred2$predictions))
})

test_that("same result with same seed, different interface", {
  ind = 1:150 %in% sample(150, 100)
  
  set.seed(2)
  mod1 = ranger(Species ~ ., data = iris[ind, ], write.forest = TRUE, num.trees = 50)
  pred1 = predict(mod1, data = iris[!ind, ])
  
  set.seed(2)
  mod3 = ranger(dependent.variable.name = "Species", data = iris[ind, ], write.forest = TRUE, num.trees = 50)
  pred3 = predict(mod3, data = iris[!ind, ])
  
  expect_that(pred1$predictions, equals(pred3$predictions))
})

test_that("no warning if data.frame has to classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_that(ranger(Species ~ ., data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("no warning if character vector in data", {
  dat <- iris
  dat$Test <- paste0("AA",as.character(1:nrow(dat)))
  expect_that(ranger(Species ~ ., data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("no error if character vector in data, prediction", {
  dat <- iris
  dat$Test <- paste0("AA",as.character(1:nrow(dat)))
  rf <- ranger(Species~., dat, write.forest = TRUE)
  expect_that(predict(rf, dat),
              not(throws_error()))
})

test_that("no warning if character vector in data, alternative interface", {
  dat <- iris
  dat$Test <- paste0("AA",as.character(1:nrow(dat)))
  expect_that(ranger(dependent.variable.name = "Species", data = dat, verbose = FALSE), 
              not(gives_warning()))
})

test_that("no error if character vector in data, alternative interface, prediction", {
  dat <- iris
  dat$Test <- paste0("AA",as.character(1:nrow(dat)))
  rf <- ranger(dependent.variable.name = "Species", data = dat, verbose = FALSE, write.forest = TRUE)
  expect_that(predict(rf, dat),
              not(throws_error()))
})

test_that("confusion matrix is of right dimension", {
  expect_that(dim(rg.class$confusion.matrix), 
              equals(rep(nlevels(iris$Species), 2)))
})

test_that("confusion matrix rows are the true classes", {
  expect_that(as.numeric(rowSums(rg.class$confusion.matrix)), 
              equals(as.numeric(table(iris$Species))))
})

test_that("case weights work", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, case.weights = rep(1, nrow(iris))), 
              not(throws_error()))
  
  ## Should only predict setosa now
  weights <- c(rep(1, 50), rep(0, 100))
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = weights, write.forest = TRUE)
  pred <- predict(rf, iris)$predictions
  expect_that(all(pred == "setosa"), is_true())
})

test_that("predict.all for classification returns numeric matrix of size trees x n", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_that(pred$predictions, is_a("matrix"))
  expect_that(dim(pred$predictions), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("predict.all for regression returns numeric matrix of size n x trees", {
  rf <- ranger(Petal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_that(pred$predictions, is_a("matrix"))
  expect_that(dim(pred$predictions), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Majority vote of predict.all for classification is equal to forest prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  
  ## Majority vote
  pred_num <- apply(pred_trees$predictions, 1, function(x) {
    which(tabulate(x) == max(tabulate(x)))
  })
  pred <- factor(pred_num, levels = 1:length(rf$forest$levels),
                 labels = rf$forest$levels)
  
  expect_that(pred, equals(pred_forest$predictions))
})

test_that("Mean of predict.all for regression is equal to forest prediction", {
  rf <- ranger(Petal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  expect_that(rowMeans(pred_trees$predictions), equals(pred_forest$predictions))
})

test_that("split select weights work", {
  expect_that(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3, 0.4)), 
              not(throws_error()))
  expect_that(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3)), 
              throws_error())
})

test_that("Tree-wise split select weights work", {
  num.trees <- 5
  weights <- replicate(num.trees, runif(ncol(iris)-1), simplify = FALSE)
  expect_that(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights), 
              not(throws_error()))
  
  weights <- replicate(num.trees+1, runif(ncol(iris)-1), simplify = FALSE)
  expect_that(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights), 
              throws_error())
})

test_that("Inbag count matrix if of right size, with replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, without replacement", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, with replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})

test_that("Inbag count matrix if of right size, without replacement, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, replace = FALSE, case.weights = runif(nrow(iris)), keep.inbag = TRUE)
  expect_that(dim(data.frame(rf$inbag.counts)), 
              equals(c(nrow(iris), rf$num.trees)))
})
