
library(survival)
library(ranger)
library(randomForestSRC)

n <- 20
num.trees <- 500

formula <- Surv(time, status) ~ .
  
src <- replicate(n, {
  rfsrc(formula, veteran, ntree = num.trees)$err.rate[num.trees]
})

std <- replicate(n, {
  ranger(formula, veteran, num.trees = num.trees)$prediction.error
})

auc <- replicate(n, {
  ranger(formula, veteran, splitrule = "auc", num.trees = num.trees)$prediction.error
})

# auc.ignore.ties <- replicate(n, {
#   ranger(formula, veteran, splitrule = "auc_ignore_ties", num.trees = num.trees)$prediction.error
# })

result <- data.frame(src, std, auc)

boxplot(result)

##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon)
##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon, splitrule = "auc")

# x <- seq(min(veteran$time),  max(veteran$time), length = 100)
# hist(veteran$time, 50, probability = TRUE)
# lines(density(veteran$time), col = "red")
# lines(x, dweibull(x, shape = 1, scale = 100), col = "blue")

