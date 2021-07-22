## Want to simulate hybridization through X generations, assuming random mating, no immigration, and no selection

source("simulate_func.R")

## New for sheep paper

########################################
## 1) Simulate different starting proportions - 0.1, 0.25, 0.5
## 20,50,100 ind total, 10 gen

## Result: very skewed starting proportions look similar to empirical results at ~ gen 5, but later generations turn to all hybrids rather than a mix of parentals and backcrosses. Starting pop size doesn't seem to matter much at these numbers.

start0.1nind20 <- simulate.hyb(nind.start=20, prop.sp1=0.1)
start0.25nind20 <- simulate.hyb(nind.start=20, prop.sp1=0.25)
start0.5nind20 <- simulate.hyb(nind.start=20, prop.sp1=0.5)

start0.1nind50 <- simulate.hyb(nind.start=50, prop.sp1=0.1)
start0.25nind50 <- simulate.hyb(nind.start=50, prop.sp1=0.25)
start0.5nind50 <- simulate.hyb(nind.start=50, prop.sp1=0.5)

start0.1nind100 <- simulate.hyb(nind.start=100, prop.sp1=0.1)
start0.25nind100 <- simulate.hyb(nind.start=100, prop.sp1=0.25)
start0.5nind100 <- simulate.hyb(nind.start=100, prop.sp1=0.5)

## 1b) Same thing but with population growth, just 50 ind. Skewed starting proportions still looks most like empirical results
start0.1grow1.25 <- simulate.hyb(nind.start=50, prop.sp1=0.1, growth.rate=1.25)
start0.25grow1.25 <- simulate.hyb(nind.start=50, prop.sp1=0.25, growth.rate=1.25)
start0.5grow1.25 <- simulate.hyb(nind.start=50, prop.sp1=0.5, growth.rate=1.25)


## 2) Simulate replicates of small sample size introductions - 10 reps, 20,50,100 ind, equal proportions from both sources, 5 or 10 gen
## Results: variation across replicates even under simplest conditions (no selection, growth, immigration, etc.). More variability across replicates with fewer individuals, fewer generations.

repsim(nrep=10, nind.start=20, n.generation=10)
repsim(nrep=10, nind.start=20, n.generation=5)

repsim(nrep=10, nind.start=50, n.generation=10)
repsim(nrep=10, nind.start=50, n.generation=5)

repsim(nrep=10, nind.start=100, n.generation=10)
repsim(nrep=10, nind.start=100, n.generation=5)

## Large sample size for contrast too

repsim(nrep=10, nind.start=1000, n.generation=10)
repsim(nrep=10, nind.start=1000, n.generation=5)

## 2b) Add population growth

repsim(nrep=10, nind.start=20, n.generation=10, growth.rate=1.25)
repsim(nrep=10, nind.start=20, n.generation=5, growth.rate=1.25)

repsim(nrep=10, nind.start=50, n.generation=10, growth.rate=1.25)
repsim(nrep=10, nind.start=50, n.generation=5, growth.rate=1.25)

repsim(nrep=10, nind.start=100, n.generation=10, growth.rate=1.25)
repsim(nrep=10, nind.start=100, n.generation=5, growth.rate=1.25)

## 2c) Add some selection. Results: selection adds to stochasticity across replicates. This is true even with a lot of individuals.
repsim(nrep=10, nind.start=20, n.generation=10, sel=0.2)
repsim(nrep=10, nind.start=20, n.generation=5, sel=0.2)

repsim(nrep=10, nind.start=50, n.generation=10, sel=0.2)
repsim(nrep=10, nind.start=50, n.generation=5, sel=0.2)

repsim(nrep=10, nind.start=100, n.generation=10, sel=0.2)
repsim(nrep=10, nind.start=100, n.generation=5, sel=0.2)


## 3) Explore different strengths of selection against hybrids. Proceeding with 50 ind as starting. Results: selection has to be pretty strong in favor of parentals to produce lots of backcrosses/parentals (similar to empirical results).

nind50sel1 <- simulate.hyb(nind.start=50, sel=1)
nind50sel0.75 <- simulate.hyb(nind.start=50, sel=0.75)
nind50sel0.5 <- simulate.hyb(nind.start=50, sel=0.5)
nind50sel0.25 <- simulate.hyb(nind.start=50, sel=0.25)
nind50sel0.1 <- simulate.hyb(nind.start=50, sel=0.1)

pdf("nind50_10gen_selvaries.pdf", width=10, height=3)

par(mfrow=c(1,5))

plot(nind50sel1[[1]][11,],nind50sel1[[2]][11,], type="n", xlab="", ylab="", main="No selection", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel1[[1]][11,], nind50sel1[[2]][11,], col="gray45", cex=1.5)

plot(nind50sel0.75[[1]][11,],nind50sel0.75[[2]][11,], type="n", xlab="", ylab="", main="Sel=0.75", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.75[[1]][11,], nind50sel0.75[[2]][11,], col="gray45", cex=1.5)

plot(nind50sel0.5[[1]][11,],nind50sel0.5[[2]][11,], type="n", xlab="", ylab="", main="Sel=0.5", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.5[[1]][11,], nind50sel0.5[[2]][11,], col="gray45", cex=1.5)

plot(nind50sel0.25[[1]][11,],nind50sel0.25[[2]][11,], type="n", xlab="", ylab="", main="Sel=0.25", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.25[[1]][11,], nind50sel0.25[[2]][11,], col="gray45", cex=1.5)

plot(nind50sel0.1[[1]][11,],nind50sel0.1[[2]][11,], type="n", xlab="", ylab="", main="Sel=0.1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.1[[1]][11,], nind50sel0.1[[2]][11,], col="gray45", cex=1.5)

## Used mtext() to allow a multi-panel plot if desired.
mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)


dev.off()


## Plot at 5 gen
pdf("nind50_5gen_selvaries.pdf", width=10, height=3)

par(mfrow=c(1,5))

plot(nind50sel1[[1]][6,],nind50sel1[[2]][6,], type="n", xlab="", ylab="", main="No selection", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel1[[1]][6,], nind50sel1[[2]][6,], col="gray45", cex=1.5)

plot(nind50sel0.75[[1]][6,],nind50sel0.75[[2]][6,], type="n", xlab="", ylab="", main="Sel=0.75", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.75[[1]][6,], nind50sel0.75[[2]][6,], col="gray45", cex=1.5)

plot(nind50sel0.5[[1]][6,],nind50sel0.5[[2]][6,], type="n", xlab="", ylab="", main="Sel=0.5", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.5[[1]][6,], nind50sel0.5[[2]][6,], col="gray45", cex=1.5)

plot(nind50sel0.25[[1]][6,],nind50sel0.25[[2]][6,], type="n", xlab="", ylab="", main="Sel=0.25", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.25[[1]][6,], nind50sel0.25[[2]][6,], col="gray45", cex=1.5)

plot(nind50sel0.1[[1]][6,],nind50sel0.1[[2]][6,], type="n", xlab="", ylab="", main="Sel=0.1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(nind50sel0.1[[1]][6,], nind50sel0.1[[2]][6,], col="gray45", cex=1.5)

## Used mtext() to allow a multi-panel plot if desired.
mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()

####### Below is simple example plotting code left as a demonstration of  how to pull out and plot specific generations of q and Q

# note that output matrix is n.generation +1 rows, so here 6 rows to include the starting individuals- hence "start0.5[[1]][6,]"
start0.5 <- simulate.hyb(0.5,5)

pdf("sim_5gen.pdf")

plot(start0.5[[1]][6,], start0.5[[2]][6,], type="n", xlab="", ylab="", main="5 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][6,], start0.5[[2]][6,], col="gray45", cex=1.5)

## Used mtext() to allow a multi-panel plot if desired.
mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()


