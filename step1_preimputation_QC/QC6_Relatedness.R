print("R script started")

pdf("6.1_relatedness.pdf")
relatedness = read.table("5_pihat_min0.125.genome", header=T)
par(pch=16, cex=1)
  with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n", main= "Relateness"))
rel = unique(relatedness$RT)
rel = rel[order(rel)]
for (i in seq_along(rel)) {
  with(subset(relatedness,RT==rel[i]) , points(Z0,Z1,col=(1+i)))
}
legend(1,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=(1+seq_along(rel)))

pdf("6.2_zoom_relatedness.pdf")
relatedness_zoom = read.table("5_zoom_pihat.genome", header=T)
par(pch=16, cex=1)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n", main= "Relateness zoom"))
rel_zoom = unique(relatedness_zoom$RT)
rel_zoom = rel_zoom[order(rel_zoom)]
for (i in seq_along(rel_zoom)) {
  with(subset(relatedness_zoom,RT==rel_zoom[i]) , points(Z0,Z1,col=(1+i)))
}
legend(0.02,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness_zoom$RT)), pch=16, col=(1+seq_along(rel_zoom)))

pdf("6.3_hist_relatedness.pdf")
relatedness = read.table("5_pihat_min0.125.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
invisible(dev.off())

print("R script finished")

