library(grf)
library(ggplot2)

mu = function(x){log(1+exp(6*x))}

n = 600
p = 20
sigma = sqrt(20)

ticks = seq(-1,1,length=n)
truth = mu(ticks)

X = matrix(runif(n*p,-1,1), nrow = n)
Y = mu(X[,1]) + sigma*rnorm(n)
  
X.test = matrix(runif(n*p, -1, 1), nrow=n)
X.test[,1] = ticks 

forest = regression_forest(X, Y, num.trees = 500)
preds.llf = predict(forest, X.test, linear.correction.variables = 1, lambda = 0.1)$predictions
preds.grf = predict(forest, X.test)$predictions
df = data.frame(cbind(ticks, truth, preds.llf))

ggplot(df, aes(ticks)) + 
  geom_point(aes(y = preds.llf, color = "LLF"), show.legend=F, size=0.6) +
  geom_line(aes(y = truth)) + 
  xlab("x") + ylab("y") + theme_bw()
ggsave("llf-bias.pdf")

ggplot(df, aes(ticks)) + 
  geom_point(aes(y = preds.grf, color = "RF-H"), show.legend=F, size=0.6) +
  geom_line(aes(y = truth)) + 
  xlab("x") + ylab("y") + theme_bw()
ggsave("grf-bias.pdf")