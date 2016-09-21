rm(list = ls())

library(ranger)
library(RColorBrewer)

setwd("~/git/split-relabel/experiments/familysize/")

load("../../FamilySize/R_files/angev80_recode_run1_line525.RData")

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
			  Y=twoa.incomem,
			  W=twoa.kidcount,
			  I=as.numeric(twoa.samesex))

# remove ~15% of data with missing father's income
# roughly 4% of fathers have zero income after removing missing
#
# only consider married women
is.ok = !is.na(twoa.incomed) & (twoa.marital==0)
DF=DF.all[is.ok,]

forest.iv = ranger(Y ~ ., DF, instrumental = TRUE, num.trees = 5000, min.node.size = 1000, mtry = 2, write.forest = TRUE, status.variable.name="W", instrument.variable.name="I", replace=FALSE, sample.fraction=0.01)

save(forest.iv, file="forest_income.raw")

agefstm.vals= c(18, 20, 22, 24)
incomed.vals = quantile(twoa.incomed, na.rm=TRUE, seq(0.025, 0.975, by = 0.05))

dummy = expand.grid(AGEFSTM=agefstm.vals, INCOMED=incomed.vals)

X.test = data.frame(
	median(twoa.agem, na.rm=TRUE),
	dummy[,1],
	12,
	0, 0, 1,
	dummy[,2])
names(X.test)=1:ncol(X.test)

preds.iv = predict(forest.iv, data.frame(X=X.test, W=-1, I=-1))$predictions

output = data.frame(dummy, TAU=preds.iv)
write.table(output, file="forest_income.out")

preds.mat.age = matrix(preds.iv, 4, 20)

cols=brewer.pal(4, "Dark2")
idx.toplot = 1:4

pdf("mother_income_vs_mother_age_at_birth_and_father_income.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(incomed.vals[-1]/1000), ylim=range(-preds.mat.age[idx.toplot,]), ylab="CATE [$1/year]", xlab="Father's Income [$1,000/year]")
for(iter in 1:4) {
	lines(incomed.vals[-1]/1000, -preds.mat.age[idx.toplot[iter],-1], lwd = 2, col = cols[iter])
}
legend("topright", sapply(c(18, 20, 22, 24), function(xx) paste(xx, "years")), lwd=2, col=cols, cex=1.5)
abline(h=0, lwd=1, lty = 2)
par=pardef
dev.off()

DF2 = DF[,1:8]
forest.Y = ranger(Y ~ ., DF2, num.trees = 500, mtry = 2, write.forest = TRUE, replace=FALSE, sample.fraction=0.05, min.node.size=100)

preds.Y = predict(forest.Y, data.frame(X=X.test))$predictions
preds.mat.age.Y = matrix(preds.Y[1:(4*20)], nrow=4, ncol=20)


pdf("mother_income_baseline_vs_mother_age_at_birth_and_father_income.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
plot(NA, NA, xlim=range(incomed.vals[-1]/1000), ylim=range(preds.mat.age.Y[idx.toplot,]/1000), ylab="Mother's Baseline Income [$1,000/year]", xlab="Father's Income [$1,000/year]")
for(iter in 1:4) {
	lines(incomed.vals[-1]/1000, preds.mat.age.Y[idx.toplot[iter],-1]/1000, lwd = 2, col = cols[iter])
}
legend("topright", sapply(c(18, 20, 22, 24), function(xx) paste(xx, "years")), lwd=2, col=cols, cex=1.5)
abline(h=0, lwd=1, lty = 2)
par=pardef
dev.off()


run.iv = function(idx) {
	slope.2 = coef(lm(Y ~ I, data = DF[idx,]))[2]
	slope.1 = coef(lm(W ~ I, data = DF[idx,]))[2]
	slope.2/slope.1
}

DF3 = DF[,-9]
forest.YI = ranger(Y ~ ., DF3, num.trees = 500, mtry = 2, write.forest = TRUE, replace=FALSE, sample.fraction=0.01, min.node.size=1000, causal=TRUE, status.variable.name="I")

preds.YI = predict(forest.YI, data.frame(X=X.test, I=-1))$predictions
preds.mat.age.YI = matrix(preds.YI[1:(8*20)], nrow=8, ncol=20)

plot(NA, NA, xlim=range(incomed.vals[-1]/1000), ylim=range(preds.mat.age.YI[idx.toplot,]), ylab="Baseline Prob. Working", xlab="Father's Income [$1k/year]")
for(iter in 1:4) {
	lines(incomed.vals[-1]/1000, preds.mat.age.YI[idx.toplot[iter],-1], lwd = 2, col = cols[iter])
}
legend("topright", sapply(c(18, 20, 22, 24), function(xx) paste(xx, "years")), lwd=2, col=cols, cex=1.5)
abline(h=0, lwd=1, lty = 2)