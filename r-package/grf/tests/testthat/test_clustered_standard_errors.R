library(grf)

set.seed(4321)

test_that("Clustered standard errors are greater than unclustered", {

    cluster_size <- 20
    # data sim
    X <- rnorm(1000)
    Y <- 5 + 2 * X + rnorm(1000)
    no_clusters <- (1:1000)

    X_cluster <- rep(X, cluster_size)
    Y_cluster <- rep(Y, cluster_size)
    clusters <- rep(1:1000, cluster_size)


    forest_corrected <- regression_forest(matrix(X_cluster),
                                          Y_cluster,
                                          num.trees = 1000,
                                          ci.group.size = 4,
                                          clusters = clusters)
    preds_corrected.oob <- predict(forest_corrected, estimate.variance = TRUE)

    forest_no_cluster <- regression_forest(matrix(X),
                                           Y,
                                           num.trees = 1000,
                                           ci.group.size = 4)
    preds_no_cluster.oob <- predict(forest_no_cluster, estimate.variance = TRUE)

    forest_uncorrected <- regression_forest(matrix(X_cluster),
                                            Y_cluster,
                                            num.trees = 1000,
                                            ci.group.size = 4)
    preds_uncorrected.oob <- predict(forest_uncorrected, estimate.variance = TRUE)
    
    forest_corrected_no_clusters <- regression_forest(matrix(X),
                                                      Y,
                                                      num.trees = 1000,
                                                      ci.group.size = 4,
                                                      clusters = no_clusters)
    preds_corrected_no_cluster.oob <- predict(forest_corrected_no_clusters, estimate.variance = TRUE)

    mean_no_cluster <- mean(preds_no_cluster.oob$variance.estimates)
    mean_uncorrected <- mean(preds_uncorrected.oob$variance.estimates)
    mean_corrected <- mean(preds_corrected.oob$variance.estimates)
    mean_corrected_no_cluster <- mean(preds_corrected_no_cluster.oob$variance.estimates)

    expect_true(mean_uncorrected < mean_corrected)
})
