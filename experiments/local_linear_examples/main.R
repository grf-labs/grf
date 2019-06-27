###########################################
## Main file to call all other functions ##
###########################################

# Running this file will create the tables and images in the "Local Linear Forests" paper.

rm(list = ls())
set.seed(1234)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("bias_image.R")
source("friedman_table.R") 
source("confidence.R")  
source("causal_table.R") 
source("boundary_bias_table.R") 