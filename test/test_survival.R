
library(survival)
library(ranger)

n <- 20

std <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran)$prediction.error
})

auc <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc")$prediction.error
})

auc.ignore.ties <- replicate(n, {
  ranger(Surv(time, status) ~ ., veteran, splitrule = "auc_ignore_ties")$prediction.error
})

result <- data.frame(std, auc, auc.ignore.ties)

boxplot(result)

##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon)
##ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon, splitrule = "auc")
