# This script was run with grf version 1.2.0.
rm(list = ls())
library(grf)
data = read.csv('angev80_recode_run1_line525.csv.xz')

FEATURES=data.frame(data$twoa.agem,
			  	data$twoa.agefstm,
			  	data$twoa.educm,
			  	as.numeric(data$twoa.blackm),
			  	as.numeric(data$twoa.hispm),
			  	as.numeric(data$twoa.othracem),
			  	twoa.incomed=round(data$twoa.incomed))
names(FEATURES)=1:ncol(FEATURES)

#labor income: twoa.incomem, worked for pay: twoa.workedm
DF.all=data.frame(
			  X=FEATURES,
			  Y=1 - as.numeric(data$twoa.workedm), # The outcome is whether the mother did not work
			  W=as.numeric(data$twoa.kidcount > 2),
			  I=as.numeric(data$twoa.samesex))

# remove ~15% of data with missing father's income
# roughly 4% of fathers have zero income after removing missing
#
# only consider married women
is.ok = !is.na(data$twoa.incomed) & (data$twoa.marital==0)
DF=DF.all[is.ok,]

agefstm.vals= c(18, 20, 22, 24)
incomed.vals = quantile(data$twoa.incomed, na.rm=TRUE, seq(0.025, 0.975, by = 0.05))

dummy = expand.grid(AGEFSTM=agefstm.vals, INCOMED=incomed.vals)
X.test = data.frame(
	median(data$twoa.agem, na.rm=TRUE),
	dummy[,1],
	12,
	0, 0, 1,
	dummy[,2])
names(X.test)=1:ncol(X.test)

DF$Y.hat = predict(lm(Y ~ ., data = DF[,1:8]))

W.lr = glm(W ~ ., data = DF[,c(1:7, 9)], family = binomial())
DF$W.hat = predict(W.lr, type="response")

print(Sys.time())

forest.iv = instrumental_forest(
	DF[,1:ncol(FEATURES)], DF$Y, DF$W, DF$I,
	Y.hat = DF$Y.hat, W.hat = DF$W.hat, Z.hat = 0.5,
	min.node.size = 600,
	num.trees = 100000,
	ci.group.size = 125,
	sample.fraction = 0.2)
preds.iv = predict(forest.iv, X.test, estimate.variance = TRUE)
tau.hat = preds.iv$predictions
var.hat = preds.iv$variance.estimates

output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
write.table(output, file="familysize.out")

print(Sys.time())
