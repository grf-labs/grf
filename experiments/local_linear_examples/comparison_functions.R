library(grf)
library(glmnet)
library(BayesTree)

run.comparison = function(data) {
  print("Beginning run")
  
  # Extract train and test data
  X = data$X
  Y = data$Y
  X.test = data$X.test
  truth = data$truth
  
  # Run random forest
  print("RF")
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = FALSE, honesty = TRUE)
  preds.grf = predict(forest, X.test)$predictions
  mse.grf = mean((preds.grf - truth)**2)
  
  # Run random forest-cv
  print("RFCV")
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = TRUE, honesty = TRUE)
  preds.grf = predict(forest, X.test)$predictions
  mse.grf.cv = mean((preds.grf - truth)**2)
  
  # Run local linear frorest
  print("LLF")
  preds.llf = predict(forest, X.test, linear.correction.variables = 1:ncol(X.test), lambda = 0.1)$predictions
  mse.llf = mean((preds.llf - truth)**2)
  
  # Local linear forest - cv
  print("LLF CV")
  preds.llf.cv = predict(forest, X.test, linear.correction.variables = 1:5, tune.lambda = FALSE)$predictions
  #preds.llf.cv = predict(forest, X.test, linear.correction.variables = 1:ncol(X.test), tune.lambda = TRUE)$predictions
  mse.llf.cv = mean((preds.llf.cv - truth)**2)
  
  # Run adaptive (non-honest) random forest 
  forest.adaptive = regression_forest(X, Y, num.trees = 1000, tune.parameters = TRUE, honesty = FALSE)
  preds.grf.adaptive = predict(forest.adaptive, X.test)$predictions
  mse.grf.adaptive = mean((preds.grf.adaptive - truth)**2)
  
  # Run BART
  bart.mod = bart(X, Y, x.test = X.test, verbose = FALSE)
  bart.preds = colMeans(bart.mod$yhat.test)
  mse.bart = mean((bart.preds - truth)**2)
  
  # Run lasso/random forest baseline
  n = nrow(X)
  inds = sample(n, size = n/2, replace = FALSE)
  # fold 1
  X1 = X[inds,]
  Y1 = Y[inds]
  # fold 2
  X2 = X[-inds,]
  Y2 = Y[-inds]
  
  # train lasso with fold 1
  lasso.mod = cv.glmnet(X1, Y1, alpha=1)
  lasso.preds = predict(lasso.mod, newx = X2, s ="lambda.min")
  lasso.resids = Y2 - lasso.preds
  
  # train random forest with fold 2
  rf = regression_forest(X2, lasso.resids, num.trees = 1000, tune.parameters = TRUE, honesty = TRUE)
  lassorf.preds = predict(rf, newdata = X.test)$predictions + predict(lasso.mod, newx = X.test, s = "lambda.min")
  mse.lasso.rf = mean((lassorf.preds - truth)**2)
  
  results = c(LLF = mse.llf,
              LLF.cv = mse.llf.cv, 
              RF.honest = mse.grf,
              RF.honest.cv = mse.grf.cv,
              RF.adaptive = mse.grf.adaptive,
              Lasso.RF = mse.lasso.rf,
              BART = mse.bart)
  
  results
  
}

check.timing = function(data){
  # Add LLF-CV
  # Extract train and test data
  X = data$X
  Y = data$Y
  X.test = data$X.test
  truth = data$truth
  
  print("RF")
  
  # Run random forest
  ptm = proc.time()
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = FALSE, honesty = TRUE)
  preds.grf = predict(forest, X.test)$predictions
  time.grf = proc.time() - ptm 
  
  print("RFCV")
  
  # Run random forest-CV
  ptm = proc.time()
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = TRUE, honesty = TRUE)
  preds.grf = predict(forest, X.test)$predictions
  time.grf.cv = proc.time() - ptm 
  
  print("LLF")
  
  # Run local linear forest 
  ptm = proc.time()
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = FALSE, honesty = TRUE)
  preds.llf = predict(forest, X.test, linear.correction.variables = 1:ncol(X), lambda = 0.1)$predictions
  time.llf = proc.time() - ptm 
  
  print("LLF CV")
  
  # Run local linear forest - cv
  ptm = proc.time()
  forest = regression_forest(X, Y, num.trees = 1000, tune.parameters = FALSE, honesty = TRUE)
  preds.llf.xc = predict(forest, X.test, linear.correction.variables = 1:ncol(X), tune.lambda = TRUE)$predictions
  time.llf.cv = proc.time() - ptm 
  
  # Run adaptive (non-honest) random forest
  ptm = proc.time()
  forest.adaptive = regression_forest(X, Y, num.trees = 1000, tune.parameters = TRUE, honesty = FALSE)
  preds.grf.adaptive = predict(forest.adaptive, X.test)$predictions
  time.grf.adaptive = proc.time() - ptm 
  
  # Run BART
  ptm = proc.time()
  bart.mod = bart(X, Y, x.test = X.test, verbose = FALSE)
  bart.preds = colMeans(bart.mod$yhat.test)
  time.bart = proc.time() - ptm 
  
  # Run lasso/random forest baseline
  ptm = proc.time()
  n = nrow(X)
  inds = sample(n, size = n/2, replace = FALSE)
  # fold 1
  X1 = X[inds,]
  Y1 = Y[inds]
  # fold 2
  X2 = X[-inds,]
  Y2 = Y[-inds]
  
  # train lasso with fold 1
  lasso.mod = cv.glmnet(X1, Y1, alpha=1)
  lasso.preds = predict(lasso.mod, newx = X2, s ="lambda.min")
  lasso.resids = Y2 - lasso.preds
  
  # train random forest with fold 2
  rf = regression_forest(X2, lasso.resids, num.trees = 1000, tune.parameters = TRUE, honesty = TRUE)
  lassorf.preds = predict(rf, newdata = X.test)$predictions + predict(lasso.mod, newx = X.test, s = "lambda.min")
  time.lasso.rf = proc.time() - ptm
  
  results = c(time.llf[3],
              time.llf.cv[3],
              time.grf[3],
              time.grf.cv[3],
              time.grf.adaptive[3],
              time.lasso.rf[3],
              time.bart[3])
  
  names(results) = c("LLF", "LLF.cv", "RF.honest", "RF.honest.cv", "RF.adaptive", "Lasso.RF", "BART")
  results 
}