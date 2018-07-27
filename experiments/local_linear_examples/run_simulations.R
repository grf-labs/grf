# Read inputs to simulation
source("comparison_functions.R")

main = function(args){
  print("running args")
  print(args)
  
  n = as.numeric(args[1])
  p = as.numeric(args[2])
  sigma = as.numeric(args[3])
  setup = as.numeric(args[4])
  experiment = as.numeric(args[5])
  num_reps = as.numeric(args[6])
  num_test = as.numeric(args[7])
  
  print("Defining setup")
  
  if(setup == 1){
    gen.data = function() {
      friedman = function(x){
        return(10*sin(pi*x[1]*x[2]) + 20*((x[3] - 0.5)**2) + 10*x[4] + 5*x[5])
      }
      
      X = matrix(runif(n*p, 0, 1), nrow = n)
      Y = apply(X, MARGIN = 1, FUN = friedman) + sigma*rnorm(n)
      
      X.test = matrix(runif(num_test*p, 0, 1), nrow = num_test)
      truth = apply(X.test, MARGIN = 1, FUN = friedman)
      
      data = list(X = X, Y = Y, X.test = X.test, truth = truth)
    }
    
  } else if(setup == 2){
    gen.data = function() {
      mu = function(x){log(1 + exp(6 * x))}
      
      X = matrix(runif(n*p, -1, 1), nrow = n)
      Y = mu(X[,1]) + sigma*rnorm(n)
      
      X.test = matrix(runif(num_test*p, -1, 1), nrow = num_test)
      truth = mu(X.test[,1])
      list(X = X, Y = Y, X.test = X.test, truth = truth)
    }
    
  }
  
  print("running expt")
  
  if (experiment == 1) {
    # MSE vs. other methods
    res = replicate(num_reps, {
      data = gen.data()
      run.comparison(data)
    })
    results = rowMeans(res)
    
  } else if (experiment == 2) {
    # timing 
    res = replicate(num_reps, {
      data = gen.data()
      check.timing(data)
    })
    results = rowMeans(res)
  }
  
  write.csv(results, file = paste0(paste("results/res", n, p, sigma, setup, experiment, sep = "-"), ".csv")) 
}

ns = 1000
ps.friedman = c(10, 50)
sigmas.friedman = c(5, 20)
setups = c(1, 2)
experiments = c(1, 2)
num_reps = 1
num_test = 100

args_friedman = expand.grid(n = ns, p = ps.friedman, sigma = sigmas.friedman, setup = 1, experiment = experiments, 
                            num_reps = num_reps, num_test = num_test)

ps.boundary = c(5, 10, 20)
sigmas.boundary = c(0.1, 1, 2)
args_boundary = expand.grid(n = ns, p = ps.boundary, sigma = sigmas.boundary, setup = 2, experiment = experiments, 
                            num_reps = num_reps, num_test = num_test)

args_full = rbind(args_friedman, args_boundary)
num_experiments = nrow(args_full)

for(i in 1:num_experiments){
  args = args_full[i,]
  main(args)
}

processed.files = data.frame(sapply(1:num_experiments, function(i){
  args = args_full[i,]
  n = as.numeric(args[1])
  p = as.numeric(args[2])
  sigma = as.numeric(args[3])
  setup = as.numeric(args[4])
  experiment = as.numeric(args[5])
  num_reps = as.numeric(args[6])
  num_test = as.numeric(args[7])
  
  fnm = paste0(paste("results/res", n, p, sigma, setup, experiment, sep = "-"), ".csv")
  errors = read.csv(fnm)
  
  c(n, p, sigma, setup, experiment, errors[,2])
}))
processed.files = t(processed.files)
colnames(processed.files) = c("n", "p", "sigma", "setup", "expt", "LLF", "LLF.cv", "RF.honest", "RF.honest.cv", "RF.adaptive", "Lasso.RF", "BART")

write.csv(processed.files, file = "summarized.results.csv", row.names = FALSE)