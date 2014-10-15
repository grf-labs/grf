
library(survival)
library(ranger)

n <- 20
num.trees <- 500

std <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, num.trees = num.trees)$prediction.error
})

auc <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc", num.trees = num.trees)$prediction.error
})

auc.ignore.ties <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc_ignore_ties", num.trees = num.trees)$prediction.error
})

result <- data.frame(std, auc, auc.ignore.ties)

boxplot(result)

##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon)
##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon, splitrule = "auc")
