
trfData <- read.csv2("trf_all_hits.csv", sep = ",", dec = ".", header = FALSE);

#head (trfData)

pdf("cent_seek_graph.pdf")

par(mfrow = c(2,1))
plot(trfData$V5~trfData$V6, xlim = c(0,600), ylim = c(0, .4*max(trfData$V5)), cex = .75, col = rgb(0,0,1,.3), pch = 16, xlab = "Repeat Length", ylab = "Repeat Frequency")
title(main = "Repeats by Length and Frequency (Long X-axis)")
plot(trfData$V5~trfData$V6, xlim = c(0,200), ylim = c(0, .4*max(trfData$V5)), cex = .75, col = rgb(0,0,1,.3), pch = 16, xlab = "Repeat Length", ylab = "Repeat Frequency")
title(main = "Repeats by Length and Frequency (Short X-axis)")

dev.off()