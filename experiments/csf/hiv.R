rm(list = ls())
library(ggplot2)
library(texreg)
library(speff2trial) # "ACTG175" data set.
library(lmtest)
library(grf)
set.seed(123)

data = ACTG175[ACTG175$arms == 1 | ACTG175$arms == 3, ]
data$preanti = as.numeric(data$preanti > 0)
X = data[c("age", "wtkg", "karnof", "cd40", "cd80", "gender", "homo",
           "race", "symptom", "drugs", "hemo", "preanti", "cd420", "cd820")]
Y = data$days
W = as.numeric(data$arms == 1) # W = 0 : ddI, W = 1: ZDV+ddI
D = data$cens

# Figure 3 - histogram overlaid
ggplot(data.frame(Y, Censored = factor(D, labels = c("Yes", "No"))), aes(x = Y, fill = Censored)) +
        geom_histogram(alpha = 0.5) +
        geom_vline(xintercept = 1000, linetype = 2, col = "red") +
        xlab("Survival time (days)") +
        ylab("Frequency") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_grey()
ggsave("HIV_histogram.pdf", width = 6, height = 5)

# Truncate Y at Y.max
Y.max = 1000
D[Y >= Y.max] = 1
Y[Y >= Y.max] = Y.max

cs.forest = causal_survival_forest(X, Y, W, D)

# BLP
full = best_linear_projection(cs.forest, X)
age = best_linear_projection(cs.forest, X[, "age", drop = F])

# Same names as in paper
varnames = c("Constant", "age", "weight", "Karnofsky score",
             "CD4 count", "CD8 count", "gender", "homosexual activity",
             "race", "symptomatic status", "intravenous drug use",
             "hemophilia", "antiretroviral history",
             "CD4 count 20+/-5 weeks", "CD8 count 20+/-5 weeks"
             )

# Table 4
texreg(list(full, age),
       custom.model.names = c("All covariates", "Age only"),
       table = FALSE,
       use.packages = FALSE,
       dcolumn = TRUE,
       single.row = TRUE,
       custom.coef.names = varnames
)

# Figure 4
X.median <- apply(X, 2, median)
age.test = seq(min(X$age), max(X$age))
X.test = matrix(rep(X.median, length(age.test)), length(age.test), byrow = TRUE)
X.test[, 1] = age.test
cs.pred = predict(cs.forest, X.test)
pt = cs.pred$predictions
pdf("HIV_data.pdf")
plot(X.test[, 1], pt, type = 'l', xlab = "Age (years)", ylab ="CATE (days)")
grid()
dev.off()
