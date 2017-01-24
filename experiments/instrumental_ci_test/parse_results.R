rm(list = ls())

library(xtable)

setwd("~/git/split-relabel/experiments/instrumental_ci_test")

filenames = list.files("output", pattern="*", full.names=TRUE)

param.names = c("nuisance", "confounding", "additive", "sparsity", "$p$", "$n$")

raw = data.frame(t(sapply(filenames, function(fnm) {
	
	res.all = read.table(fnm, sep=",", header=FALSE)
	res.mean = mean(res.all[,1])
	res = sprintf("%.2f", round(res.mean, 2))
	
	params = strsplit(fnm, "-")[[1]][c(3, 5, 6, 7, 9, 10)]
	params[1] = params[1] != 0
	params[2] = params[2] != 0
	
	params.print = as.character(params)
	params.print[params.print=="FALSE"] = "no"
	params.print[params.print=="TRUE"] = "yes"
	
	ret = c(params.print, res)
		
})))

names(raw) = c(param.names, "coverage")
rownames(raw) = 1:nrow(raw)

raw = raw[,-c(1, 2)]
raw = raw[order(raw[,1], decreasing=TRUE),]

collapsed = cbind(raw[1:8, 2:5], "coverage"=raw[9:16, 5])

xtab = xtable(collapsed, align=c("r", "|", rep("c", 3), "|", rep("c", 2), "|"))
print(xtab, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, hline.after = c(-1, 0, 4, 8), file = "simulation_results.tex")