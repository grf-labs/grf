D = read.table("harder.dat", header = TRUE)
OUT = read.table("ranger_out.confusion")
plot(NA, NA, xlim = range(D[,1]), ylim = range(OUT), xlab = "x1", ylab = "quantiles")

ord = order(D[,1])

lines(D[ord,1], 1 + 9 * (D[ord, 1] > 0), lwd = 2, col = 2)
lines(D[ord,1], 0 + 0 * (D[ord, 1] > 0), lwd = 2, col = 2)
lines(D[ord,1], -1 - 9 * (D[ord, 1] > 0), lwd = 2, col = 2)

points(D[ord,1], OUT[ord,1])
points(D[ord,1], OUT[ord,2])
points(D[ord,1], OUT[ord,3])

