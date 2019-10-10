library(grf)
library(glmnet)
library(hte)
packageVersion("hte")

simulation <- function(n = 5000, d = 2, version = 2){
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  noise <- rnorm(n=n)
  main_effect <- rep(0, times = n)
  treatment_propensity <- 0.5
  W <- matrix(rbinom(n = n,
                     size = 1,
                     prob = treatment_propensity),
              nrow = n, ncol = 1)
  if(version == 2){
    zeta1 <- 1 + 1/(1+exp(-20*(X[,1]-(1/3))))
    zeta2 <- 1 + 1/(1+exp(-20*(X[,2]-(1/3))))
  }
  else{
    zeta1 <- 2/(1+exp(-20*(X[,1]-(1/3))))
    zeta2 <- 2/(1+exp(-20*(X[,2]-(1/3))))
  }
  Tau <- matrix(zeta1 * zeta2,nrow=n, ncol=1)
  Y <- matrix(main_effect + W * Tau + noise, nrow=n, ncol=1)
  output <- list(X= X, W = W, Y = Y, Tau = Tau)
  return(output)
}

num.reps = 50
lambdas = c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5)
ns = seq(200, 1200, by = 200)

run_simulation = function(ns, lambdas, num.reps, version){
  results = sapply(ns, function(n){
    basic.results = replicate(num.reps, {
      dat = simulation(n = n, d = 20, version = version)
      dat.test = simulation(n = 1000, d = 20, version = version)

      xl_bart <- X_BART(feat = list(dat$X), tr = as.numeric(dat$W), yobs = as.numeric(dat$Y))
      cate_esti_bart = EstimateCate(xl_bart, dat.test$X)
      err = mean((cate_esti_bart - dat.test$Tau)**2)

      forest = causal_forest(as.matrix(dat$X), as.numeric(dat$Y), as.numeric(dat$W), tune.parameters = "all")
      preds = predict(forest, as.matrix(dat.test$X))$predictions
      grf.err = mean((preds - dat.test$Tau)**2)

      lasso.mod = cv.glmnet(as.matrix(dat$X), as.numeric(dat$Y), alpha = 1)
      selected = which(coef(lasso.mod) != 0)
      if(length(selected) < 2){
        selected = 1:ncol(dat$X)
      }else{
        selected = selected[-1] - 1 # remove intercept
      }

      # predict: tuning done automatically unless lambda is specified
      preds.llf = predict(forest, as.matrix(dat.test$X), linear.correction.variables = selected, ll.weight.penalty = TRUE)$predictions
      err.llf = mean((preds.llf - dat.test$Tau)**2)

      return(c(sqrt(err), sqrt(grf.err), sqrt(err.llf)))
    })
    basic.results = data.frame(t(basic.results))
    colMeans(basic.results)
  })
  results = data.frame(t(results))
  colnames(results) = c("X-BART", "CF", "LLF")
  results$n = ns
  results
}

results_version1 = run_simulation(ns, lambdas, num.reps, version = 1)
results_version2 = run_simulation(ns, lambdas, num.reps, version = 2)

write.csv(results_version1,"version1.csv", row.names = FALSE)
write.csv(results_version2,"version2.csv", row.names = FALSE)
