set.seed(42)
rm(list = ls())
library(grf)

get_data = function(name = c("lung", "veteran", "pbc", "heart", "rotterdam"),
                    drop.na = TRUE) {
  name = match.arg(name)
  if (name == "lung")  {
    df = survival::lung
    if (drop.na) {
      ix = complete.cases(df)
      df = df[ix, ]
    }
    Y = df$time
    D = as.integer(df$status > 1)
    X = as.matrix(df[ ,!(names(df) %in% c("time", "status"))])
  }
  if (name == "veteran") {
    df = survival::veteran
    if (drop.na) {
      ix = complete.cases(df)
      df = df[ix, ]
    }
    Y = df$time
    D = df$status
    dfX = df[ ,!(names(df) %in% c("time", "status"))]
    X = model.matrix(~ . -1 , dfX)
  }
  if (name == "pbc") {
    df = survival::pbc
    if (drop.na) {
      ix = complete.cases(df)
      df = df[ix, ]
    }
    Y = df$time
    D = as.integer(df$status > 0)
    df$sex = ifelse(df$sex == "f", 1, 0)
    X = as.matrix(df[ ,!(names(df) %in% c("time", "status"))])
  }
  if (name == "heart") {
    df = survival::heart
    if (drop.na) {
      ix = complete.cases(df)
      df = df[ix, ]
    }
    Y = df$stop - df$start
    D = df$event
    df$transplant = as.integer(df$transplant)
    X = as.matrix(df[ ,!(names(df) %in% c("start", "stop", "event", "id"))])
  }
  if (name == "rotterdam") {
    df = survival::rotterdam
    if (drop.na) {
      ix = complete.cases(df)
      df = df[ix, ]
    }
    Y = df$dtime
    D = df$death
    df$size = as.integer(df$size)
    X = as.matrix(df[ ,!(names(df) %in% c("dtime", "death"))])
  }

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

datasets = c("lung", "veteran", "pbc", "heart", "rotterdam")
nsim = 250

out = list()
for (name in datasets) {
  for (sim in 1:nsim) {
    data = get_data(name)
    Y = data$Y
    D = data$D
    X = data$X

    seed = sim
    sf.exact = survival_forest(X, Y, D, fast.logrank = FALSE, seed = seed, prediction.type = "Nelson-Aalen", num.trees = 500)
    sf.approx = survival_forest(X, Y, D, fast.logrank = TRUE, seed = seed, prediction.type = "Nelson-Aalen", num.trees = 500)

    df = data.frame(
      metric = c("concordance", "IBS"),
      error = c(get_cerror(sf.exact) - get_cerror(sf.approx), get_IBS(sf.exact) - get_IBS(sf.approx)),
      sim = sim,
      name = name)
    out = c(out, list(df))
  }
}
out.df = do.call(rbind, out)

aggregate(
  list(value = out.df$error),
  by = list(
    metric = out.df$metric, name = out.df$name
  ),
  FUN = mean
)

# Result
library(ggplot2)

ggplot(out.df[out.df$metric == "concordance", ], aes(y = error)) +
  geom_boxplot() +
  facet_wrap(~ name, nrow = 1) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  ylab(NULL) +
  ggtitle("Difference in concordance") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  )
ggsave("concordance.pdf", width = 5, height = 3, scale = 1)

ggplot(out.df[out.df$metric == "IBS", ], aes(y = error)) +
  geom_boxplot() +
  facet_wrap(~ name, nrow = 1) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  ylab(NULL) +
  ggtitle("Difference in integrated Brier score") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  )
ggsave("IBS.pdf", width = 5, height = 3, scale = 1)

write.csv(out.df, "benchmark.csv", row.names = FALSE)


# Basic timings repeated
library(microbenchmark)
library(xtable)

fbench = function(data, fast.logrank) {
  survival_forest(data$X, data$Y, data$D, fast.logrank = fast.logrank,
                  compute.oob.predictions = FALSE,
                  num.trees = 1, sample.fraction = 1, mtry = ncol(data$X))
}

times = 100
out = list()
for (name in datasets) {
  data = get_data(name)

  bench.exact = microbenchmark(fbench(data, FALSE), times = times, unit = "milliseconds")
  time.exact = summary(bench.exact)$mean
  bench.approx = microbenchmark(fbench(data, TRUE), times = times, unit = "milliseconds")
  time.approx = summary(bench.approx)$mean

  diff = time.exact - time.approx
  ratio = time.exact / time.approx
  df = data.frame(
    metric = c("exact(ms)", "approx(ms)", "difference(ms)", "speedup.factor"),
    value = c(time.exact, time.approx, diff, ratio),
    name = name)
  out = c(out, list(df))
}
out.df = do.call(rbind, out)

tab.df = reshape(
  out.df,
  idvar = c("name"),
  timevar = "metric",
  direction = "wide"
)

print(xtable(tab.df[-4]), include.rownames = FALSE)
