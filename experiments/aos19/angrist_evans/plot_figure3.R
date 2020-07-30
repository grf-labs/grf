rm(list = ls())
df = read.table("familysize.out")

df$INCOMED = df$INCOMED / 1000
df$ub = df$TAU + sqrt(df$VAR) * qnorm(0.975)
df$lb = df$TAU - sqrt(df$VAR) * qnorm(0.975)

png("plot_figure3.png", width = 900, height = 400)
par(mfrow = c(1, 2))
p.df = df[df$AGEFSTM == 18, ]
plot(p.df$INCOMED, p.df$TAU,
     ylim = c(min(p.df$lb), max(p.df$ub)),
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = " Mother 18 years old at first birth")
lines(p.df$INCOMED, p.df$lb, lty = 2)
lines(p.df$INCOMED, p.df$ub, lty = 2)
abline(h=0, lty = 3)

p.df = df[df$AGEFSTM == 22, ]
plot(p.df$INCOMED, p.df$TAU,
     ylim = c(min(p.df$lb), max(p.df$ub)),
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = "Mother 22 years old at first birth")
lines(p.df$INCOMED, p.df$lb, lty = 2)
lines(p.df$INCOMED, p.df$ub, lty = 2)
abline(h=0, lty = 3)
dev.off()
