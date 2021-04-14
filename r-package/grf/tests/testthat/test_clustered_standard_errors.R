library(grf)

set.seed(4)

test_that("Clustered standard errors are greater than unclustered", {
  n <- 200
  p <- 4
  cluster_size <- 10
  # data sim
  X <- matrix(rnorm(n * p), n, p)
  MU <- 5 + 2 * X[, 1]
  Y <- MU + rnorm(n)
  no_clusters <- (1:n)

  X_cluster <- Reduce(rbind, lapply(1:cluster_size, function(cc) X))
  Y_cluster <- rep(Y, cluster_size)
  clusters <- rep(1:n, cluster_size)
  MU_clusters <- rep(MU, cluster_size)


  forest_corrected <- regression_forest(X_cluster,
    Y_cluster,
    ci.group.size = 4,
    clusters = clusters
  )
  preds_corrected.oob <- predict(forest_corrected, estimate.variance = TRUE)

  forest_no_cluster <- regression_forest(X,
    Y,
    ci.group.size = 4
  )
  preds_no_cluster.oob <- predict(forest_no_cluster, estimate.variance = TRUE)

  forest_uncorrected <- regression_forest(X_cluster,
    Y_cluster,
    ci.group.size = 4
  )
  preds_uncorrected.oob <- predict(forest_uncorrected, estimate.variance = TRUE)

  forest_corrected_no_clusters <- regression_forest(X,
    Y,
    ci.group.size = 4,
    clusters = no_clusters
  )
  preds_corrected_no_cluster.oob <- predict(forest_corrected_no_clusters, estimate.variance = TRUE)

  mean_no_cluster <- mean(preds_no_cluster.oob$variance.estimates)
  mean_uncorrected <- mean(preds_uncorrected.oob$variance.estimates)
  mean_corrected <- mean(preds_corrected.oob$variance.estimates)
  mean_corrected_no_cluster <- mean(preds_corrected_no_cluster.oob$variance.estimates)

  expect_equal(mean_no_cluster, mean_corrected, tolerance = 0.2 * mean_no_cluster)
  expect_equal(mean_no_cluster, mean_corrected_no_cluster, tolerance = 0.2 * mean_no_cluster)
  expect_gt(mean_no_cluster, 2 * mean_uncorrected)

  mse_no_cluster <- mean((preds_no_cluster.oob$predictions - MU)^2)
  mse_uncorrected <- mean((preds_uncorrected.oob$predictions - MU_clusters)^2)
  mse_corrected <- mean((preds_corrected.oob$predictions - MU_clusters)^2)
  mse_corrected_no_cluster <- mean((preds_corrected_no_cluster.oob$predictions - MU)^2)

  expect_equal(mse_no_cluster, mse_corrected, tolerance = 0.16 * mse_no_cluster)
  expect_equal(mse_no_cluster, mse_corrected_no_cluster, tolerance = 0.1 * mse_no_cluster)
  expect_lt(mse_no_cluster, 2 * mse_uncorrected)
})

test_that("Clustered predictions are reasonable with unevenly sized clusters", {
  n <- 200
  p <- 4
  cluster_size <- 10

  # data sim
  X <- matrix(rnorm(n * p), n, p)
  MU <- 5 + 2 * X[, 1]
  Y <- MU + 4 * ((1:n) > n / 2) + rnorm(n)
  no_clusters <- (1:n)

  X_cluster <- rbind(
    Reduce(rbind, lapply(1:cluster_size, function(cc) X[1:(n / 2), ])),
    X[(n / 2 + 1):n, ]
  )
  Y_cluster <- c(rep(Y[1:(n / 2)], cluster_size), Y[(n / 2 + 1):n])
  clusters <- c(rep(1:(n / 2), cluster_size), (n / 2 + 1):n)
  MU_clusters <- c(rep(MU[1:(n / 2)], cluster_size), MU[(n / 2 + 1):n])

  forest_corrected <- regression_forest(X_cluster,
    Y_cluster,
    ci.group.size = 1,
    clusters = clusters,
    equalize.cluster.weights = TRUE
  )
  preds_corrected.oob <- predict(forest_corrected, estimate.variance = FALSE)

  forest_no_cluster <- regression_forest(X,
    Y,
    ci.group.size = 1
  )
  preds_no_cluster.oob <- predict(forest_no_cluster, estimate.variance = FALSE)

  forest_uncorrected <- regression_forest(X_cluster,
    Y_cluster,
    ci.group.size = 1
  )
  preds_uncorrected.oob <- predict(forest_uncorrected, estimate.variance = FALSE)

  forest_corrected_no_clusters <- regression_forest(X,
    Y,
    ci.group.size = 1,
    clusters = no_clusters,
    equalize.cluster.weights = TRUE
  )
  preds_corrected_no_cluster.oob <- predict(forest_corrected_no_clusters, estimate.variance = FALSE)

  mse_no_cluster <- mean((preds_no_cluster.oob$predictions - MU)^2)
  mse_uncorrected <- mean((preds_uncorrected.oob$predictions - MU_clusters)[n / 2 * (cluster_size - 1) + 1:n]^2)
  mse_corrected <- mean((preds_corrected.oob$predictions - MU_clusters)[n / 2 * (cluster_size - 1) + 1:n]^2)
  mse_corrected_no_cluster <- mean((preds_corrected_no_cluster.oob$predictions - MU)^2)

  expect_equal(mse_no_cluster, mse_corrected, tolerance = 0.1 * mse_no_cluster)
  expect_equal(mse_no_cluster, mse_corrected_no_cluster, tolerance = 0.1 * mse_no_cluster)
  expect_gt(mse_no_cluster, 2 * mse_uncorrected)

  meanp_no_cluster <- mean(preds_no_cluster.oob$predictions)
  meanp_uncorrected <- mean(preds_uncorrected.oob$predictions[n / 2 * (cluster_size - 1) + 1:n])
  meanp_corrected <- mean(preds_corrected.oob$predictions[n / 2 * (cluster_size - 1) + 1:n])
  meanp_corrected_no_cluster <- mean(preds_corrected_no_cluster.oob$predictions)

  expect_equal(meanp_no_cluster, meanp_corrected, tolerance = 0.2)
  expect_equal(meanp_no_cluster, meanp_corrected_no_cluster, tolerance = 0.2)
  expect_gt(meanp_no_cluster, meanp_uncorrected + 1)
})
