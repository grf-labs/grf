instrumental.forest <- function(X, Y, W, Z,
	sample.fraction=0.5,
	mtry=ceiling(ncol(X)/3),
	num.trees=500,
	num.threads=NULL,
	min.node.size=NULL,
	keep.inbag = FALSE,
	honesty=TRUE,
	precompute.nuisance=TRUE,
	split.regularization=0,
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

     if (!is.numeric(split.regularization) | split.regularization < 0 | 
        split.regularization > 1) {
        stop("Error: Invalid value for split.regularization. Please give a value in [0,1].")
     }
     
     if (!precompute.nuisance) {
     	
     	input.data <- as.matrix(cbind(X, Y, W, Z))	
     	
     } else {
     	
     	forest.Y <- regression.forest(X, Y,
			sample.fraction=sample.fraction,
			mtry=mtry,
			num.trees=min(500, num.trees),
			num.threads=num.threads,
			min.node.size=NULL,
			keep.inbag = FALSE,
			honesty=TRUE,
			seed=seed)
		Y.hat = predict(forest.Y)
		
		forest.W <- regression.forest(X, W,
			sample.fraction=sample.fraction,
			mtry=mtry,
			num.trees=min(500, num.trees),
			num.threads=num.threads,
			min.node.size=NULL,
			keep.inbag = FALSE,
			honesty=TRUE,
			seed=seed)
		W.hat = predict(forest.W)
		
		forest.Z <- regression.forest(X, Z,
			sample.fraction=sample.fraction,
			mtry=mtry,
			num.trees=min(500, num.trees),
			num.threads=num.threads,
			min.node.size=NULL,
			keep.inbag = FALSE,
			honesty=TRUE,
			seed=seed)
		Z.hat = predict(forest.Z)
		
		input.data <- as.matrix(cbind(X, Y - Y.hat, W - W.hat, Z - Z.hat))

     }

	
	variable.names <- c(colnames(X), "outcome", "treatment", "instrument")
	outcome.index.zeroindexed <- ncol(X)
	treatment.index.zeroindexed <- ncol(X) + 1
	instrument.index.zeroindexed <- ncol(X) + 2
	
	no.split.variables <- numeric(0)

	forest <- instrumental_train(input.data,
						outcome.index.zeroindexed,
						treatment.index.zeroindexed,
						instrument.index.zeroindexed,
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
						honesty,
						split.regularization)
						
						
	class(forest) <- "instrumental.forest"
	forest
}

predict.instrumental.forest <- function(forest, newdata, num.threads = NULL) {

	if (is.null(num.threads)) {
        num.threads <- 0
     } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
     }
     
	 sparse.data <- as.matrix(0)
	 variable.names <- character(0)
	
	 input.data <- as.matrix(cbind(newdata, NA))	
	
     instrumental_predict(forest,
     				  input.data,
     				  sparse.data,
     				  variable.names,
     				  num.threads)    				 
} 
