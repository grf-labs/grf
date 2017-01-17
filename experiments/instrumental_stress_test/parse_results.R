rm(list = ls())

library(xtable)

setwd("~/git_local/split-relabel/experiments/instrumental_stress_test")

filenames = list.files("output", pattern="*", full.names=TRUE)

res.names = c("cent. forest", "plain forest", "kNN", "series")
param.names = c("nuisance", "confounding", "additive", "sparsity", "$p$", "$n$")

raw = data.frame(t(sapply(filenames, function(fnm) {
	
	res.all = read.table(fnm, sep=",", header=FALSE)
	res = rowMeans(res.all)
	res.keep = c(res[c(1, 2, 3)], min(res[5:8], na.rm=TRUE))
	
	best.idx = which.min(res.keep)
	res.print = sprintf("%.2f", round(res.keep, 2))
	res.print[best.idx] = paste("\\bf", res.print[best.idx])
	
	params = strsplit(fnm, "-")[[1]][c(3, 5, 6, 7, 9, 10)]
	params[1] = params[1] != 0
	params[2] = params[2] != 0
	
	params.print = as.character(params)
	params.print[params.print=="FALSE"] = "no"
	params.print[params.print=="TRUE"] = "yes"
	
	ret = c(params.print, res.print)
		
})))

names(raw) = c(param.names, res.names)
rownames(raw) = 1:nrow(raw)