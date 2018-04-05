rm(list = ls())

library(xtable)

setwd("~/git/grf/experiments/instrumental_ci_test")

filenames = list.files("output", pattern="*", full.names=TRUE)

param.names = c("nuisance", "confounding", "additive", "sparsity", "$p$", "$n$")

raw = data.frame(t(sapply(filenames, function(fnm) {
	
	res.all = read.table(fnm, sep=",", header=FALSE)[,1]
	res.mean = mean(res.all)
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
raw = raw[order(as.numeric(as.character(raw[,4]))),]
raw = raw[order(as.numeric(as.character(raw[,3]))),]
raw = raw[order(as.numeric(as.character(raw[,2]))),]
raw = raw[order(raw[,1], decreasing=TRUE),]

collapsed = cbind(raw[1:24, 2:5], "coverage"=raw[25:48, 5])

xtab = xtable(collapsed, align=c("r", "|", rep("c", 3), "|", rep("c", 2), "|"))
print(xtab, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity,
      hline.after = c(-1, -1, 0, 0, 4, 8, 12, 12, 16, 20, 24, 24),
      file = "simulation_results_true.tex")


raw = data.frame(t(sapply(filenames, function(fnm) {
  
  res.all = read.table(fnm, sep=",", header=FALSE)[,2]
  res.mean = mean(res.all)
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
raw = raw[order(as.numeric(as.character(raw[,4]))),]
raw = raw[order(as.numeric(as.character(raw[,3]))),]
raw = raw[order(as.numeric(as.character(raw[,2]))),]
raw = raw[order(raw[,1], decreasing=TRUE),]

collapsed = cbind(raw[1:24, 2:5], "coverage"=raw[25:48, 5])

xtab = xtable(collapsed, align=c("r", "|", rep("c", 3), "|", rep("c", 2), "|"))
print(xtab, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity,
      hline.after = c(-1, -1, 0, 0, 4, 8, 12, 12, 16, 20, 24, 24),
      file = "simulation_results_avg.tex")


