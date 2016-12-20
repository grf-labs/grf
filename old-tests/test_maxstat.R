library(ranger)
library(maxstatRF)
library(survival)
library(reshape2)
library(ggplot2)

# Compare survival curves -------------------------------------------------
idx <- sample(nrow(veteran), 2/3*nrow(veteran), replace = FALSE)
train_data <- veteran[idx, ]
test_data <- veteran[-idx, ]

num_trees <- 50
mtry <- 3
minprop <- 0.1
alpha <- 0.5
min_node_size <- 5
replace <- TRUE
pmethod <- "minLau"

rf_maxstatRF <- maxstatRF(Surv(time, status) ~., train_data, 
                          num_trees = num_trees, mtry = mtry, minprop = minprop, 
                          alpha = alpha, min_node_size = min_node_size, replace = replace, 
                          pmethod = pmethod)

rf_ranger <- ranger(Surv(time, status) ~., train_data, splitrule = "maxstat",
                    write.forest = TRUE,
                    num.trees = num_trees, mtry = mtry, minprop = minprop, 
                    alpha = alpha, min.node.size = min_node_size, replace = replace)

pred_maxstatRF <- rf_maxstatRF$predict(test_data)
pred_ranger <- predict(rf_ranger, test_data)$survival

colnames(pred_maxstatRF) <- rf_maxstatRF$timepoints
colnames(pred_ranger) <- timepoints(rf_ranger)
df <- rbind(data.frame(Package = "maxstatRF", melt(pred_maxstatRF, value.name = "Survival", varnames = c("ID", "Time"))), 
            data.frame(Package = "ranger", melt(pred_ranger, value.name = "Survival", varnames = c("ID", "Time"))))

ggplot(df, aes(x = Time, y = Survival, color = Package)) + 
  geom_line() + 
  facet_wrap(~ID)

# Compare Brier score -------------------------------------------------
library(ranger)
library(maxstatRF)
library(survival)
library(prodlim)
library(pec)

## Pec function for maxstatRF
predictSurvProb.ForestSurvival <- function(object, newdata, times, ...) {
  pred <- object$predict(newdata)
  pos <- sapply(times, function(x) {
    max(c(0, which(x >= object$timepoints)))
  })    
  p <- cbind(1, pred)[, pos+1, drop = FALSE]
  return(p)
}

## Pec function for ranger
predictSurvProb.ranger <- function(object, newdata, times, ...) {
  pred <- predict(object, newdata)
  pos <- sapply(times, function(x) {
    max(c(0, which(x >= pred$unique.death.times)))
  })    
  p <- cbind(1, predictions(pred))[, pos+1, drop = FALSE]
  return(p)
}

## Parameters
num_trees <- 50
mtry <- 3
minprop <- 0.1
alpha <- 0.5
min_node_size <- 5
replace <- TRUE
pmethod <- "minLau"
formula <- Surv(time, status) ~.

compare_packages <- function() {
  ## Data
  idx <- sample(nrow(veteran), 2/3*nrow(veteran), replace = FALSE)
  train_data <- veteran[idx, ]
  test_data <- veteran[-idx, ]
  time <- sort(unique(train_data$time))
  
  ## Run
  rf_maxstatRF <- maxstatRF(formula, train_data, 
                            num_trees = num_trees, mtry = mtry, minprop = minprop, 
                            alpha = alpha, min_node_size = min_node_size, replace = replace, 
                            pmethod = pmethod)
  
  rf_ranger <- ranger(formula, train_data, splitrule = "maxstat",
                      write.forest = TRUE,
                      num.trees = num_trees, mtry = mtry, minprop = minprop, 
                      alpha = alpha, min.node.size = min_node_size, replace = replace)
  
  ## Compare
  fit.pec <- pec(list(maxstatRF = rf_maxstatRF, ranger = rf_ranger), 
                 formula = formula, data = test_data, 
                 times = time, cens.model = "marginal", reference = FALSE)
  crps(fit.pec)[, 1]
}

res <- replicate(100, compare_packages())
boxplot(t(res))
