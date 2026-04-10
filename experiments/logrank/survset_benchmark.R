# Update these paths accordingly:
Sys.setenv(RETICULATE_PYTHON = "~/mamba/bin/python")
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste0("~/mamba/bin:", old_path))

library(grf)
library(reticulate)
library(xtable)

surv_set_data = import("SurvSet.data")
loader = surv_set_data$SurvLoader()

tmp = tempfile(fileext = ".csv")
loader$df_ds$to_csv(tmp, index = FALSE)
df_ds = read.csv(tmp)
datasets = df_ds |>
  subset(n >= 1000) |>
  subset(is_td == "False")
datasets = datasets$ds
length(datasets) # 23

get_data = function(name) {
  options(na.action = "na.pass")
  pd = loader$load_dataset(ds_name=name)$df
  tmp = tempfile(fileext = ".csv")
  pd$to_csv(tmp, index = FALSE)
  df = read.csv(tmp)
  df$pid = NULL
  is.char = sapply(df, is.character)
  df[is.char] = lapply(df[is.char], as.factor)
  Y = df$time; df$time = NULL
  D = df$event; df$event = NULL
  X = model.matrix(~ . - 1, df)

  list(Y = Y, D = D, X = X)
}


get_cerror = function(forest) {
  1 - survival::concordance(Surv(forest$Y.orig, forest$D.orig) ~
                              rowSums(-log(forest$predictions)), reverse = TRUE)$concordance
}

get_IBS = function(forest) {
  pred = predict(forest)
  survex::integrated_brier_score(
    survival::Surv(forest$Y.orig, forest$D.orig),
    surv = pred$predictions,
    times = pred$failure.times)
}

# Appendix survset bench
out = list()
for (name in datasets) {
  data = get_data(name); cat(name, "\n")
  Y = data$Y
  D = data$D
  X = data$X

  sf.exact = survival_forest(X, Y, D, fast.logrank = FALSE, seed = 42, prediction.type = "Nelson-Aalen", num.trees = 500)
  sf.approx = survival_forest(X, Y, D, fast.logrank = TRUE, seed = 42, prediction.type = "Nelson-Aalen", num.trees = 500)

  sf.exact0 = survival_forest(X, Y, D, fast.logrank = FALSE, seed = 42, prediction.type = "Nelson-Aalen", num.trees = 500, alpha = 0)
  sf.approx0 = survival_forest(X, Y, D, fast.logrank = TRUE, seed = 42, prediction.type = "Nelson-Aalen", num.trees = 500, alpha = 0)

  df = data.frame(
    Data.set = paste0("$\\emph{", name, "}$"),
    n = nrow(X),
    p = ncol(X),
    M = length(unique(Y[D==1])),
    metric = c("$\\Delta PE_C$", "$\\Delta PE_{IBS}$"),
    value = c(get_cerror(sf.exact) - get_cerror(sf.approx), get_IBS(sf.exact) - get_IBS(sf.approx)),
    value0 = c(get_cerror(sf.exact0) - get_cerror(sf.approx0), get_IBS(sf.exact0) - get_IBS(sf.approx0))
    )
  out = c(out, list(df))
}
out.df = do.call(rbind, out)
out.df
write.csv(out.df, "survset_benchmark.csv", row.names = FALSE)

print(xtable(out.df, digits = 5),
      sanitize.text.function = identity,
      include.rownames = FALSE,
      format.args = list(big.mark = " ", decimal.mark = "."))
