
library(survival)
library(ranger)

ranger(Surv(time, status) ~ ., veteran)
ranger(Surv(time, status) ~ ., veteran, splitrule = "auc")

ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon)
ranger(Surv(time, status) ~ age + sex + perfor + surg + etype, colon, splitrule = "auc")
