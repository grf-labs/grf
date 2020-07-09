# For series
library(AER)
library(splines)

# For KNN
library(FNN)

iv.series = function(X, Y, W, Z, X.test, df = 3, interact = FALSE) {
	
	X.all = rbind(X, X.test)
	X.spl = Reduce(cbind , lapply(data.frame(X.all), function(xx) ns(xx, df)))
	
	if (interact) {
		X.reg = model.matrix( ~ . * . + 0, data = data.frame(X.spl))
	} else {
		X.reg = X.spl
	}

	if(ncol(X.reg) >= length(Y)) return(rep(NA, nrow(X.test)))

	series = ivreg(Y ~ X.reg[1:length(Y),] * W | X.reg[1:length(Y),] * Z)
	beta = coef(series)
	
	tau.hat = beta[ncol(X.reg) + 2] + X.reg[length(Y) + 1:nrow(X.test),] %*% beta[ncol(X.reg) + 2 + 1:ncol(X.reg)]
	
	tau.hat
	
}

iv.knn = function(X, Y, W, Z, X.test, k = 100) {

	neighbors = get.knnx(X, X.test, k)$nn.index
	
	apply(neighbors, 1, function(nn) {
		y.hat = mean(Y[nn])
		w.hat = mean(W[nn])
		z.hat = mean(Z[nn])
		yz.hat = mean(Y[nn] * Z[nn])
		wz.hat = mean(W[nn] * Z[nn])
		
		(yz.hat - y.hat * z.hat) / (wz.hat - w.hat * z.hat)
	})
}


