rm(list = ls())

library(ranger)

load("~/git/split-relabel/FamilySize/R_files/angev80_recode_run1_line525.RData")

FEATURES=data.frame(twoa.agem,
			  	twoa.agefstm,
			  	twoa.educm,
			  	as.numeric(twoa.blackm),
			  	as.numeric(twoa.hispm),
			  	as.numeric(twoa.othracem),
			  	twoa.incomed=round(twoa.incomed))
names(FEATURES)=1:ncol(FEATURES)

#labor income: twoa.incomem, worked for pay: twoa.workedm
DF.all=data.frame(
			  X=FEATURES,
			  Y=as.numeric(twoa.workedm),
			  W=twoa.kidcount,
			  I=as.numeric(twoa.samesex))

# remove ~15% of data with missing father's income
# roughly 4% of fathers have zero income after removing missing
#
# only consider married women
is.ok = !is.na(twoa.incomed) & (twoa.marital==0)
DF=DF.all[is.ok,]

#mat = model.matrix(~ . + 0, DF)
#write.table(mat, file="~/git/split-relabel/experiments/familysize/input.data", quote=FALSE, row.names=FALSE)

forest.iv = ranger(Y ~ ., DF, instrumental = TRUE, num.trees = 2000, min.node.size = 1000, mtry = 4, write.forest = TRUE, status.variable.name="W", instrument.variable.name="I", replace=FALSE, sample.fraction=0.1)

X.test = data.frame(
	median(twoa.agem, na.rm=TRUE),
	median(twoa.agefstm, na.rm=TRUE),
	quantile(twoa.educm, seq(0, 1, by = 0.05)),
	0, 0, 1,
	median(twoa.incomed, na.rm=TRUE))
names(X.test)=1:ncol(X.test)

preds.iv = predict(forest.iv, data.frame(X=X.test, W=-1, I=-1))$predictions
