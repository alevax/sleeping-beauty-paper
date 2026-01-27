args <- commandArgs(trailingOnly = TRUE)
libdata <- read.table(args[2], header=TRUE, sep=",", row.names="Library")
pdf(args[4])
boxplot (libdata[, 2:8], main="Library Sequences", ylab="Raw Count", xlab="Sequence Type", names=c(1:7))
legend("topright", c("1:Barcoded sequences", "2:Sequences with mutagen", "3:Mappable sequences", "4:Unique, mappable sequences", "5:Sequences mapped to reference", "6:Unique sequences mapped to reference", "7:Non redundant insertions"), cex=0.7)

regiondata <- read.table(args[3], header=TRUE, sep=",", row.names="Library")
sortedregiondata <- regiondata[order(regiondata[,1], decreasing=TRUE),]
sortedregionmatrix <- data.matrix(sortedregiondata)
barplot(height=t(sortedregionmatrix[c(1:10,(nrow(sortedregionmatrix)-9):nrow(sortedregionmatrix)),1:2]), ylab="Raw Count", beside=TRUE, col=c("blue", "red"), main="Mappable Sequences by Library (Top 10 vs. Bottom 10)", las=2, cex.names=0.7)
legend("topright", c("Mappable sequences", "Sequences with a mapping"), fill=c("blue", "red"), cex=0.7)
abline(v=30.5, lty=2)
axis(1, lab=F, at=(c(0:19)*3)+2)

matrixcalcs <- cbind(sortedregionmatrix, sortedregionmatrix[,2]/sortedregionmatrix[,1])
sortedregionmatrix <- matrixcalcs[order(matrixcalcs[,7], decreasing=TRUE), ]
barplot(height=t(sortedregionmatrix[c(1:10,(nrow(sortedregionmatrix)-9):nrow(sortedregionmatrix)),1:2]), ylab="Raw Count", beside=TRUE, main="Successful Mappings by Library (Top 10 vs. Bottom 10)", las=2, cex.names=0.7, col=c("blue", "red"))
legend("topright", c("Mappable sequences", "Sequences with a mapping"), fill=c("blue", "red"), cex=0.7)
abline(v=30.5, lty=2)
axis(1, lab=F, at=(c(0:19)*3)+2)

matrixcalcs <- cbind(matrixcalcs, matrixcalcs[,6]/matrixcalcs[,1])
sortedregionmatrix <- matrixcalcs[order(matrixcalcs[,8], decreasing=TRUE), ]
barplot(height=t(sortedregionmatrix[c(1:10,(nrow(sortedregionmatrix)-9):nrow(sortedregionmatrix)),c(1,6)]), ylab="Raw Count", beside=TRUE, main="Redundancy of Mapped Regions by Library (Top 10 vs. Bottom 10)", las=2, cex.names=0.7, col=c("blue", "red"))
legend("topright", c("Mappable sequences", "Distinct regions mapped"), fill=c("blue", "red"), cex=0.7)
abline(v=30.5, lty=2)
axis(1, lab=F, at=(c(0:19)*3)+2)

dev.off()


