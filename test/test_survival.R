
library(survival)
library(ranger)

n <- 100

std <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, num.trees = 2000)$prediction.error
})

auc <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc", num.trees = 2000)$prediction.error
})

auc.ignore.ties <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc_ignore_ties", num.trees = 2000)$prediction.error
})

result <- data.frame(std, auc, auc.ignore.ties)

boxplot(result)

##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon)
##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon, splitrule = "auc")
