###########################################
## Main file to call all other functions ##
###########################################

rm(list = ls())
set.seed(1234)

source("bias_image.R")
source("friedman_table.R") 
source("confidence.R")  
source("causal_table.R") 
source("boundary_bias_table.R")
source("wages.R")