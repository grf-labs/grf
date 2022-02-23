# The following script reproduces the HIV application from
# the manuscript.
rm(list = ls())
library(ggplot2)
library(texreg)
library(speff2trial) # "ACTG175" data set.
library(grf)
library(xtable)
set.seed(123)

data = ACTG175[ACTG175$arms == 1 | ACTG175$arms == 3, ]
data$preanti = as.numeric(data$preanti > 0)
X = data[c("age", "wtkg", "karnof", "cd40", "cd80", "gender", "homo",
           "race", "symptom", "drugs", "hemo", "preanti", "cd420", "cd820")]
Y = data$days
W = as.numeric(data$arms == 1) # W = 0 : ddI, W = 1: ZDV+ddI
D = data$cens

# Overlaid histogram with T.max
ggplot(data.frame(Y, Censored = factor(D, labels = c("Yes", "No"))), aes(x = Y, fill = Censored)) +
        geom_histogram(alpha = 0.5) +
        geom_vline(xintercept = 1000, linetype = 1, col = "red") +
        xlab("Survival time (days)") +
        ylab("Frequency") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_grey()
ggsave("HIV_histogram.pdf", width = 6, height = 5)

# Truncate Y at Y.max
Y.max = 1000

cs.forest = causal_survival_forest(X, Y, W, D, horizon = Y.max, num.trees = 10000, ci.group.size = 12)

# Estimates and SEs for a random subset of individuals
idx = sample(nrow(X), 10)
pp = predict(cs.forest, estimate.variance = TRUE)
vimp = variable_importance(cs.forest)
colnames(X)[order(vimp)[1:4]]
df = data.frame(
        CATE = pp$predictions[idx],
        CATE.se = sqrt(pp$variance.estimates[idx]),
        hemophilia = ifelse(X[idx, "hemo"] == 1, "Yes", "No"),
        gender = ifelse(X[idx, "gender"] == 1, "Male", "Female"),
        homosexual.activity = ifelse(X[idx, "homo"] == 1, "Yes", "No"),
        antiretroviral.history = ifelse(X[idx, "preanti"] == 1, "Experienced", "Naive")
)
print(xtable(df[order(df$CATE), ]
             ), include.rownames = FALSE)

# BLP
full = best_linear_projection(cs.forest, X)
age = best_linear_projection(cs.forest, X[, "age", drop = F])

# Same names as in paper
varnames = c("Constant", "Age", "Weight", "Karnofsky score",
             "CD4 count", "CD8 count", "Gender", "Homosexual activity",
             "Race", "Symptomatic status", "Intravenous drug use",
             "Hemophilia", "Antiretroviral history",
             "CD4 count 20+/-5 weeks", "CD8 count 20+/-5 weeks"
             )

# BLP Table
texreg(list(full, age),
       custom.model.names = c("All covariates", "Age only"),
       table = FALSE,
       use.packages = FALSE,
       dcolumn = TRUE,
       single.row = TRUE,
       custom.coef.names = varnames
)

# CATE plot
X.median <- apply(X, 2, median)
age.test = seq(min(X$age), max(X$age))
X.test = matrix(rep(X.median, length(age.test)), length(age.test), byrow = TRUE)
X.test[, 1] = age.test
cs.pred = predict(cs.forest, X.test, estimate.variance = TRUE)
pt = cs.pred$predictions
ub = pt + sqrt(cs.pred$variance.estimates) * qnorm(0.975)
lb = pt - sqrt(cs.pred$variance.estimates) * qnorm(0.975)
pdf("HIV_data.pdf")
plot(X.test[, 1], pt, type = 'l', xlab = "Age (years)", ylab = "CATE (days)", ylim = c(min(lb), max(ub)))
lines(X.test[, 1], ub, lty = 2)
lines(X.test[, 1], lb, lty = 2)
grid()
dev.off()
