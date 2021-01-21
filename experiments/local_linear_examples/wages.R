rm(list = ls())
set.seed(123)

library(grf)
library(glmnet)
library(BART)
library(xgboost)
library(ggplot2)
library(dplyr)
library(splines)
library(rlearner)

# load the data
load("cps1976_2018.RData")
data = data.frame(cps1976_2018)

# clean data: remove NA and missingg outcomes
data = data[!is.na(data$incwage),]
data = data[data$incwage != 9999999,] # remove missing
data = data[data$incwage != 9999998,] # remove missing

# extract 2018 data
data = data[data$year == 2018,]

# add age^2, educ^2 covariates
data$agesq = data$age**2
data$educsq = data$educ**2

covariates = c("age", "agesq", "educ", "educsq", "occ2010", "occ10ly", "sex", "race",
               "marst", "labforce", "ind1950", "classwkr", "wkstat", "uhrswork1", "metro", "famsize")

continuous.covariates = which(covariates %in% c("age", "educ", "agesq", "educsq", "uhrswork1", "famsize"))
outcome = "incwage"

data = data[,c(covariates, outcome)]
data = data[complete.cases(data),]

# transform outcome (standard for wage regressions)
data$incwage = log(data$incwage + 1)

######################
## ERROR EVALUATION ##
######################

num.reps = 100
size.test = 10000
sample.sizes = c(2000, 5000, 10000, 50000)

mse.sample.sizes = data.frame(t(sapply(sample.sizes, function(size){
  index.train = sample(1:nrow(data), size = size, replace = FALSE)

  X = data[index.train, covariates]
  Y = data$incwage[index.train]

  X = X[complete.cases(X),]
  Y = Y[complete.cases(X)]

  results = data.frame(t(sapply(1:num.reps, function(i){
    print(i)

    index.test = sample((1:nrow(data))[-index.train], size = size.test, replace = FALSE)

    X.test = data[index.test, covariates]
    truth = data$incwage[index.test]

    X.test = X.test[complete.cases(X.test),]
    truth = truth[complete.cases(X.test)]

    forest = regression_forest(as.matrix(X), Y, honesty = TRUE, tune.parameters = "all")

    llf.preds = predict(forest, as.matrix(X.test),
                              linear.correction.variables = continuous.covariates,
                              ll.weight.penalty = T)$predictions
    llf.mse = mean((llf.preds - truth)**2)

    rf.preds = predict(forest, as.matrix(X.test))$predictions
    rf.mse = mean((rf.preds - truth)**2)

    ols.form = as.formula(paste("Y", paste(covariates, collapse = "+"), sep = "~"))
    dd.ols = cbind(Y, X)
    ols.fit = lm(ols.form, dd.ols)
    ols.preds = predict(ols.fit, X.test)
    ols.mse = mean((ols.preds - truth)**2)

    mm = model.matrix( ~.^2, data = X)
    lasso.mod = cv.glmnet(mm, Y, alpha = 1)
    mmtest = model.matrix( ~.^2, data = X.test)
    lasso.preds = predict(lasso.mod, newx = mmtest, lambda= lasso.mod$lambda.min)
    lasso.mse = mean((lasso.preds - truth)**2)

    bart.mod = wbart(X, Y, X.test)
    bart.preds = bart.mod$yhat.test.mean
    bart.mse = mean((bart.preds-truth)**2)

    boost.cv.fit = rlearner::cvboost(as.matrix(X), Y)
    xgb.preds = predict(boost.cv.fit, as.matrix(X.test))
    xg.mse = mean((xgb.preds - truth)**2)

    return(c(ols.mse, llf.mse, rf.mse, lasso.mse, xg.mse, bart.mse))
  })))

  mses = colMeans(results)
  sds = apply(results, MARGIN = 2, FUN = sd)

  as.numeric(c(mses, sds))
})))

colnames(mse.sample.sizes) = c("OLS", "LLF", "RF", "Lasso", "XG", "BART",
                               "OLS.sd", "LLF.sd", "RF.sd", "Lasso.sd", "XG.sd", "BART.sd")

mse.sample.sizes$size = sample.sizes

write.csv(mse.sample.sizes,"wages_sample_sizes.csv", row.names = FALSE)
