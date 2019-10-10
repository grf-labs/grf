########################
## basic image script ##
########################

library(grf)
library(glmnet)
library(ggplot2)

mu = function(x){
  log(1+exp(6*x))
}

n = 600
p = 20
sigma = sqrt(20)
X = matrix(runif(n*p,-1,1), nrow = n)

ticks = sort(X[,1])
truth = mu(ticks)
X[,1] = ticks
Y = truth + sigma*rnorm(n)

forest = regression_forest(X,Y, num.trees = 1000,honesty=TRUE,tune.parameters="all")
preds = predict(forest, linear.correction.variables = 1, tune.lambda = TRUE,
                estimate.variance = TRUE)
df = data.frame(ticks = X[,1], predictions = preds$predictions,
                upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
ggplot(df, aes(ticks)) +
  geom_point(aes(y=predictions, color="RF-H"), show.legend=F, size=0.6) +
  geom_line(aes(y=truth)) +
  geom_line(aes(y=upper),color="gray",lty=2) +
  geom_line(aes(y=lower),color="gray",lty=2) +
  xlab("x") + ylab("y") + theme_bw()
ggsave("llf-preds.pdf")

preds = predict(forest, estimate.variance = TRUE)
df = data.frame(ticks = X[,1], predictions = preds$predictions,
                upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
ggplot(df, aes(ticks)) +
  geom_point(aes(y=predictions, color="RF-H"), show.legend=F, size=0.6) +
  geom_line(aes(y=truth)) +
  geom_line(aes(y=upper),color="gray",lty=2) +
  geom_line(aes(y=lower),color="gray",lty=2) +
  xlab("x") + ylab("y") + theme_bw()
ggsave("grf-preds.pdf")

################################
## coverage, length, and RMSE ##
################################

p = 5
sigma = sqrt(20)
ns = c(500, 2000)
num_reps = 50

linear.mod = function(x){
  return(10 / (1 + exp(-10 * (x[1] - 0.5))) + 5 / (1 + exp(-10 * (x[2] - 0.5))))
}

mu = function(x){
  log(1+exp(6*x))
}

friedman = function(x){
  return(10*sin(pi*x[1]*x[2]) + 20*((x[3] - 0.5)**2) + 10*x[4] + 5*x[5])
}

full_results = sapply(ns, function(n){
  results = replicate(num_reps, {
    X = matrix(runif(n*p,0,1), nrow = n)
    truth = apply(X, MARGIN = 1, FUN = friedman)
    Y = truth + sigma*rnorm(n)

    percent_llf = 0
    percent_grf = 0
    avg_grf = 0
    avg_llf = 0

    lasso.mod = cv.glmnet(X,Y,alpha=1)
    selected = which(coef(lasso.mod)!=0)
    if(length(selected) < 2){
      selected = 1:p
    }else{
      selected = selected[-1] - 1 # remove intercept
    }

    forest = regression_forest(X, Y, num.trees = 2000, ci.group.size = 4,
                                 sample.fraction = 0.5,
                                 honesty = TRUE,
                                 tune.parameters = TRUE)
    preds = predict(forest, linear.correction.variables = selected, tune.lambda=TRUE, estimate.variance = TRUE)
    mse.llf = mean((preds$predictions - truth)**2)
    df = data.frame(ticks = X[,1], predictions= preds$predictions,
                    upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                    lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
    for(i in 1:n){
      xtick = df$ticks[i]
      xlow = df$lower[i]
      xup = df$upper[i]
      truthi = truth[i]
      if(xlow <= truthi && truthi <= xup){
        percent_llf = percent_llf + 1;
      }
      avg_llf = avg_llf + abs(xup - xlow)
    }
    percent_llf = percent_llf/n
    avg_llf = avg_llf/n

    forest = regression_forest(X, Y, num.trees = 2000, ci.group.size = 4,
                               sample.fraction = 0.5,
                               honesty = TRUE,
                               tune.parameters = TRUE)
    preds = predict(forest, estimate.variance = TRUE)
    mse.grf = mean((preds$predictions - truth)**2)

    df = data.frame(ticks = X[,1], predictions= preds$predictions,
                    upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                    lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
    for(i in 1:n){
      xtick = df$ticks[i]
      xlow = df$lower[i]
      xup = df$upper[i]
      truthi = truth[i]
      if(xlow <= truthi && truthi <= xup){
        percent_grf = percent_grf + 1;
      }
      avg_grf = avg_grf + abs(xup - xlow)
    }
    percent_grf = percent_grf/n
    avg_grf = avg_grf/n

    return(c(sqrt(mse.grf), sqrt(mse.llf), percent_grf, percent_llf, avg_grf, avg_llf))
  })
  results = data.frame(t(results))
  round(colMeans(results),3)
})
full_results = data.frame(t(full_results))
colnames(full_results) = c("rmse.grf", "rmse.llf", "percent.grf", "percent.llf", "avg.grf", "avg.llf")
full_results$n = ns
full_results

write.table(full_results,"friedman_confidence_p5.csv", row.names = FALSE)

p=20
full_results = sapply(ns, function(n){
  results = replicate(num_reps, {
    X = matrix(runif(n*p,0,1), nrow = n)
    truth = apply(X, MARGIN = 1, FUN = friedman)
    Y = truth + sigma*rnorm(n)

    percent_llf = 0
    percent_grf = 0
    avg_grf = 0
    avg_llf = 0

    lasso.mod = cv.glmnet(X,Y,alpha=1)
    selected = which(coef(lasso.mod)!=0)
    if(length(selected) < 2){
      selected = 1:p
    }else{
      selected = selected[-1] - 1 # remove intercept
    }

    forest = regression_forest(X, Y, num.trees = 2000, ci.group.size = 2,
                                 sample.fraction = 0.5,
                                 honesty = TRUE,
                                 tune.parameters = TRUE)
    preds = predict(forest, linear.correction.variables = selected, tune.lambda=TRUE, estimate.variance = TRUE)
    mse.llf = mean((preds$predictions - truth)**2)
    df = data.frame(ticks = X[,1], predictions= preds$predictions,
                    upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                    lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
    for(i in 1:n){
      xtick = df$ticks[i]
      xlow = df$lower[i]
      xup = df$upper[i]
      truthi = truth[i]
      if(xlow <= truthi && truthi <= xup){
        percent_llf = percent_llf + 1;
      }
      avg_llf = avg_llf + abs(xup - xlow)
    }
    percent_llf = percent_llf/n
    avg_llf = avg_llf/n

    forest = regression_forest(X, Y, num.trees = 2000, ci.group.size = 2,
                               sample.fraction = 0.5,
                               honesty = TRUE,
                               tune.parameters = TRUE)
    preds = predict(forest, estimate.variance = TRUE)
    mse.grf = mean((preds$predictions - truth)**2)

    df = data.frame(ticks = X[,1], predictions= preds$predictions,
                    upper = preds$predictions + 1.96*sqrt(preds$variance.estimates),
                    lower = preds$predictions - 1.96*sqrt(preds$variance.estimates))
    for(i in 1:n){
      xtick = df$ticks[i]
      xlow = df$lower[i]
      xup = df$upper[i]
      truthi = truth[i]
      if(xlow <= truthi && truthi <= xup){
        percent_grf = percent_grf + 1;
      }
      avg_grf = avg_grf + abs(xup - xlow)
    }
    percent_grf = percent_grf/n
    avg_grf = avg_grf/n

    return(c(mse.grf, mse.llf, percent_grf, percent_llf, avg_grf, avg_llf))
  })
  results = data.frame(t(results))
  round(colMeans(results),3)
})
full_results = data.frame(t(full_results))
colnames(full_results) = c("mse.grf", "mse.llf", "percent.grf", "percent.llf", "avg.grf", "avg.llf")
full_results$n = ns

write.csv(full_results,"friedman_confidence_p20", row.names = FALSE)
