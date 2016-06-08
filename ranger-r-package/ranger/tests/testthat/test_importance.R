## Tests for importance measures

library(ranger)
context("ranger_imp")

## Initialize the random forests
rg.imp <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                 importance = "impurity")
rg.perm <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                 importance = "permutation")
rg.scale.perm <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                 importance = "permutation", scale.permutation.importance = TRUE)

## Tests
test_that("importance measures work", {
  expect_is(rg.imp$variable.importance, "numeric")
  expect_is(rg.perm$variable.importance, "numeric")
  expect_is(rg.scale.perm$variable.importance, "numeric")
})

test_that("gini importance is larger than 1", {
  expect_gt(rg.imp$variable.importance[1], 1)
})

test_that("unscaled importance is smaller than 1", {
  expect_lt(rg.perm$variable.importance[1], 1)
})

test_that("scaled importance is larger than 1", {
  expect_gt(rg.scale.perm$variable.importance[1], 1)
})

