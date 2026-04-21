_This folder has replication files for the paper "Efficient Log-Rank Updates for Random Survival Forests" by Sverdrup, Yang, and LeBlanc._

* Table 1 and Table 3: `timings.R`

* Figure 1 & Table 2: `simulation_benchmark.R`

* Figure 2 & Table 2: `simulation_rmse.R`

These scripts rely on the packages
```
grf, survival, survex, microbenchmark, ggplot2, xtable
```

* Appendix Table 4: `survset_benchmark.R` (requires the `reticulate` package along with the Python package `SurvSet`; the script uses a default [mamba](https://github.com/mamba-org/mamba) environment)
