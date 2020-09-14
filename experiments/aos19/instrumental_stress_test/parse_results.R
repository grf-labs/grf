rm(list = ls())

library(xtable)

setwd("~/git/grf/experiments/instrumental_stress_test")

filenames = list.files("output", pattern="*", full.names=TRUE)

res.names = c("cent. forest", "plain forest", "series", "kNN")
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

raw = raw[order(raw[,3], decreasing=TRUE),]
raw = raw[order(raw[,1], decreasing=FALSE),]

raw = raw[,c(1, 3, 2, 4:10)]

raw2 = raw[,c(1:6, 10:7)]

collapsed = cbind(raw2[1:32, 2:10], raw2[33:64, 7:10])

xtab = xtable(collapsed, align=c("r", "|", "|", rep("c", 2), "|", rep("r", 3), "|", "|", rep("c", 4), "|", rep("c", 4), "|", "|"))
print(xtab, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, hline.after = c(-1, -1, 0, 0, 4, 8, 12, 16, 16, 20, 24, 28, 32, 32), file = "simulation_results.tex")
