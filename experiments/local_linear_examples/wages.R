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

  results = data.frame(t(sapply(1:num.reps, function(i){
    index.train = sample(1:nrow(data), size = size, replace = FALSE)
    index.test = sample((1:nrow(data))[-index.train], size = size.test, replace = FALSE)
    
    X = data[index.train, covariates]
    Y = data$incwage[index.train]
    
    X = X[complete.cases(X),]
    Y = Y[complete.cases(X)]
    
    forest = regression_forest(as.matrix(X), Y, honesty = TRUE, tune.parameters = TRUE)
    
    X.test = data[index.test, covariates]
    truth = data$incwage[index.test]
    
    X.test = X.test[complete.cases(X.test),]
    truth = truth[complete.cases(X.test)]
    
    ll.lambda = tune_ll_regression_forest(forest, linear.correction.variables = continuous.covariates, 
                                          ll.weight.penalty = TRUE)$lambda.min
    llf.lasso.preds = predict(forest, as.matrix(X.test), 
                              linear.correction.variables = continuous.covariates,
                              ll.lambda = ll.lambda, 
                              ll.weight.penalty = TRUE)$predictions
    llf.mse = mean((llf.lasso.preds - truth)**2)
    
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
    
    boost.cv.fit = cvboost(as.matrix(X), Y)
    xgb.preds = predict(boost.cv.fit, as.matrix(X.test))
    xg.mse = mean((xgb.preds - truth)**2)
    
    return(c(ols.mse, llf.mse, rf.mse, lasso.mse, xg.mse, bart.mse))
  })))
  
  mses = colMeans(results)
  sds = apply(results, MARGIN = 2, FUN = sd)
  
  as.numeric(c(mses, sds))
})))

colnames(mse.sample.sizes) = c("OLS", "Lasso", "RF", "LLF", "XG", "BART",
                               "OLS.sd", "Lasso.sd", "RF.sd", "LLF.sd", "XG.sd", "BART.sd")

mse.sample.sizes$size = sample.sizes

write.csv(mse.sample.sizes,"wages_sample_sizes.csv", row.names = FALSE)

############################
## CALIBRATION AND IMAGES ##
############################

train.index = sample(nrow(data), size = floor(nrow(data)/2))
train = data[train.index,]
test = data[-train.index,]

subsets = c("age", "race", "educ", "famsize")
range.dummies = c("age.r", "race.r", "educ.r", "famsize.r")
test$age.r = test$age < 20 | test$age > 85
test$race.r = test$race > 200
test$educ.r = test$educ > 111
test$famsize.r = test$famsize > 6

methods = c("lasso", "llf", "ols", "rf")

forest.fit = regression_forest(train[,covariates], train$incwage, tune.parameters = T)
ols.mod = lm(incwage ~ ., data = train[,c(covariates, "incwage")])
mm = model.matrix( ~.^2, data = train[,covariates])
lasso.mod = cv.glmnet(mm, train$incwage, alpha = 1)

# generate plots for sparse regions and find corresponding t-statistics
t_statistics = data.frame(t(sapply(1:length(subsets), function(i){
  
  test.subset = subsets[i]
  range.subset = range.dummies[i]
  testdata = test[test[,range.subset] == T,]
  
  names = sapply(methods, function(method){
    paste(paste(method, test.subset, sep = "_"), "jpeg", sep = ".")
  })
  
  rf.preds = predict(forest.fit, testdata[,covariates])$predictions
  
  tl = tune_ll_regression_forest(forest.fit, linear.correction.variables = continuous.covariates,
                                 ll.weight.penalty = TRUE)$lambda.min
  llf.preds = predict(forest.fit, testdata[,covariates], linear.correction.variables = continuous.covariates,
                      ll.lambda = tl, ll.weight.penalty = TRUE)$predictions
  
  ols.preds = predict(ols.mod, newdata = testdata)
  
  mmtest = model.matrix( ~.^2, data=testdata[,covariates])
  lasso.preds = predict(lasso.mod, newx = mmtest, lambda= lasso.mod$lambda.min)
  
  dt = data.frame("Wages" = testdata$incwage,
                  "LLF" = llf.preds,
                  "RF" = rf.preds,
                  "OLS" = ols.preds,
                  "Lasso" = lasso.preds)
  colnames(dt)[5] = "Lasso"
  
  spline.fit = lm(dt$Wages ~ ns(dt$Lasso, df = 3))
  dt$spline = spline.fit$fitted.values
  ggplot(dt, aes(x = Lasso, y = Wages)) + geom_point(cex = 1) + 
    ylab("Observed log wages") + xlab("Predicted log wages") + 
    geom_line(aes(x = Lasso, y = spline), color = "red", lty = 1, lwd = 1.5) + 
    geom_abline(color = "Darkgray", intercept = 0, slope = 1, lty = 2, lwd = 1.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
    xlim(c(range(dt$Lasso, dt$Wages))) + ylim(c(range(dt$Lasso, dt$Wages)))
  ggsave(names[1])
  
  spline.fit = lm(dt$Wages ~ ns(dt$LLF, df = 3))
  dt$spline = spline.fit$fitted.values
  ggplot(dt, aes(x = LLF, y = Wages)) + geom_point(cex = 1) + 
    ylab("Observed log wages") + xlab("Predicted log wages") + 
    geom_line(aes(x = LLF, y = spline), color = "red", lty = 1, lwd = 1.5) + 
    geom_abline(color = "Darkgray", intercept = 0, slope = 1, lty = 2, lwd = 1.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
    xlim(c(range(dt$LLF, dt$Wages))) + ylim(c(range(dt$LLF, dt$Wages)))
  ggsave(names[2])
  
  spline.fit = lm(dt$Wages ~ ns(dt$RF, df = 3))
  dt$spline = spline.fit$fitted.values
  ggplot(dt, aes(x = RF, y = Wages)) + geom_point(cex = 1) + 
    ylab("Observed log wages") + xlab("Predicted log wages") + 
    geom_line(aes(x = RF, y = spline), color = "red", lty = 1, lwd = 1.5) + 
    geom_abline(color = "Darkgray", intercept = 0, slope = 1, lty = 2, lwd = 1.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
    xlim(c(range(dt$RF, dt$Wages))) + ylim(c(range(dt$RF, dt$Wages)))
  ggsave(names[3])
  
  spline.fit = lm(dt$Wages ~ ns(dt$OLS, df = 3))
  dt$spline = spline.fit$fitted.values
  ggplot(dt, aes(x = OLS, y = Wages)) + geom_point(cex = 1) + 
    ylab("Observed log wages") + xlab("Predicted log wages") + 
    geom_line(aes(x = OLS, y = spline), color = "red", lty = 1, lwd = 1.5) + 
    geom_abline(color = "Darkgray", intercept = 0, slope = 1, lty = 2, lwd = 1.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
    xlim(c(range(dt$OLS, dt$Wages))) + ylim(c(range(dt$OLS, dt$Wages)))
  ggsave(names[4])
  
  # squared differences: 
  ols.sds = (ols.preds - testdata$incwage)**2
  llf.sds = (llf.preds- testdata$incwage)**2
  lasso.sds = (lasso.preds - testdata$incwage)**2
  rf.sds = (rf.preds - testdata$incwage)**2
  
  # corresponding paired t-statistics: 
  t.ols = mean(ols.sds - llf.sds)/(sd(ols.sds - llf.sds)/sqrt(nrow(testdata)))
  t.lasso = mean(lasso.sds - llf.sds)/(sd(lasso.sds - llf.sds)/sqrt(nrow(testdata)))
  t.rf = mean(rf.sds - llf.sds)/(sd(rf.sds - llf.sds)/sqrt(nrow(testdata)))
  
  c(t.ols, t.lasso, t.rf)
})))

colnames(t_statistics) = c("OLS", "Lasso", "RF")
rownames(t_statistics) = subsets

write.csv(t_statistics,"t_statistics.csv", row.names = TRUE)
