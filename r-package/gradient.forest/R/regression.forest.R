regression.forest <- function(X, Y,
	sample.fraction=0.5,
	mtry=ceiling(ncol(X)/3),
	num.trees=500,
	num.threads=NULL,
	min.node.size=NULL,
	keep.inbag = FALSE,
	honesty=TRUE,
	seed=NULL) {

	sparse.data <- as.matrix(0)

	if (is.null(mtry)) {
        mtry <- 0
     } else if (!is.numeric(mtry) | mtry < 0) {
        stop("Error: Invalid value for mtry")
     }
     
     verbose=FALSE
    
     if (is.null(num.threads)) {
        num.threads <- 0
     } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
     }
    
     if (is.null(min.node.size)) {
        min.node.size <- 0
     } else if (!is.numeric(min.node.size) | min.node.size < 0) {
        stop("Error: Invalid value for min.node.size")
     }
     
     sample.with.replacement <- FALSE
     
     if (!is.logical(keep.inbag)) {
        stop("Error: Invalid value for keep.inbag")
     }
     
     if (!is.numeric(sample.fraction) | sample.fraction <= 0 | 
        sample.fraction > 1) {
        stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
     }

     if (is.null(seed)) {
        seed <- runif(1, 0, .Machine$integer.max)
     }

	input.data <- as.matrix(cbind(X, Y))	
	variable.names <- c(colnames(X), "outcome")
	outcome.index <- ncol(input.data)
	outcome.index.zeroindexed <- outcome.index - 1
	no.split.variables <- numeric(0)

	forest <- regression_train(input.data,
						outcome.index.zeroindexed,
						sparse.data,
						variable.names,
						mtry,
						num.trees,
						verbose,
						num.threads,
						min.node.size,
						sample.with.replacement,
						keep.inbag,
						sample.fraction,
						no.split.variables,
						seed,
						honesty)
						
	forest[["original.data"]] <- input.data
	class(forest) <- "regression.forest"
	forest
}

predict.regression.forest <- function(forest, newdata = NULL, num.threads = NULL) {

	if (is.null(num.threads)) {
    	num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
    	stop("Error: Invalid value for num.threads")
    }
     
	sparse.data <- as.matrix(0)
	variable.names <- character(0)
	
	forest.short <- forest[-which(names(forest) == "original.data")]
	
	if (!is.null(newdata)) {
		input.data <- as.matrix(cbind(newdata, NA))
		regression_predict(forest.short,
	    				   input.data,
	     				   sparse.data,
	     				   variable.names,
	     				   num.threads)
	 } else {
	 	input.data <- forest[["original.data"]]
		regression_predict_oob(forest.short,
	     				       input.data,
	     				       sparse.data,
	     				       variable.names,
	     				       num.threads)
	 }
} 
