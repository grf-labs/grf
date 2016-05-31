## Tests for print function

library(ranger)
context("ranger_print")

## Initialize the random forest
rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)

## Test print ranger function
expect_that(print(rf), prints_text("Ranger result"))

## Test print forest function
expect_that(print(rf$forest), prints_text("Ranger forest object"))

## Test print prediction function
expect_that(print(predict(rf, iris)), prints_text("Ranger prediction"))

## Test str ranger function
expect_that(str(rf), prints_text("List of 13"))

## Test str forest function
expect_that(str(rf$forest), prints_text("List of 10"))
