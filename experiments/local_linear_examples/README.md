This folder has replication files for the paper "Local Linear Forests", by Friedberg, Tibshirani, Athey, and Wager, as well as a general introduction for practitioners.

The R markdown introduction to local linear forests can be found at llf_guide.Rmd. 
Required R packages: ggplot2, glmnet, grf.

All selected settings for the algorithms here are described in the paper.
Required R packages: BART, ggplot2, glmnet, grf, hte, xgboost.

Contained in this folder are:
* main.R, which runs all simulations found in the paper. Note that running time
(especially with tuning and the full replications) may be long. This file calls each of the files listed below.
* bias_image.R, replicates Figure 1 of the paper.
* boundary_bias_table.R, replicates errors for equation 1 (Table 5).
* causal_table.R, replicates MSE on treatment effects (Table 3).
* confidence.R, replicates confidence interval coverage and length (Table 2).
* friedman_table.R, replicates errors on Friedman's simulation (Equation 7; Table 4 in the Appendix) and Figure 4.
