## Want to simulate hybridization through X generations, assuming random mating, no immigration, and no selection

source("simulate_func.R")

## Simulations for 1000 ind, 11 generations, starting proportions of 0.5, 0.75, 0.9

start0.5 <- simulate.hyb(1000,0.5,11)
start0.25 <- simulate.hyb(1000,0.25,11)
start0.05 <- simulate.hyb(1000,0.05,11)
start0.1 <- simulate.hyb(1000,0.1,12)

## make this fig including colors

##pdf("sim_5gen_1000ind.pdf", width=12,height=4)
pdf("sim_5gen_10gen_1000ind.pdf", width=8,height=6)

par(mfrow=c(2,3))

plot(start0.5[[1]][6,], start0.5[[2]][6,], type="n", xlab="", ylab="", main="5 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][6,], start0.5[[2]][6,], col="gray45", cex=1.5)

plot(start0.25[[1]][6,], start0.25[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][6,], start0.25[[2]][6,], col="gray45", cex=1.5)

plot(start0.1[[1]][6,], start0.1[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][6,], start0.1[[2]][6,], col="gray45", cex=1.5)

plot(start0.5[[1]][11,], start0.5[[2]][11,], type="n", xlab="", ylab="", main="10 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][11,], start0.5[[2]][11,], col="gray45", cex=1.5)

plot(start0.25[[1]][11,], start0.25[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][11,], start0.25[[2]][11,], col="gray45", cex=1.5)

plot(start0.1[[1]][11,], start0.1[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][11,], start0.1[[2]][11,], col="gray45", cex=1.5)

mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()


## Add some migrants

source("simulate_func.R")

## Simulations for 1000 ind, 11 generations, starting proportions of 0.5, 0.75, 0.9

start0.5 <- simulate.hyb(1000,0.5,12,imm.sp1=0.25)
start0.25 <- simulate.hyb(1000,0.25,12,imm.sp1=0.1,imm.sp2=0.1)
start0.05 <- simulate.hyb(1000,0.05,11)
start0.1 <- simulate.hyb(1000,0.1,12)

## make this fig including colors

##pdf("sim_5gen_1000ind.pdf", width=12,height=4)
pdf("sim_5gen_10gen_1000ind.pdf", width=8,height=6)

par(mfrow=c(2,3))

plot(start0.5[[1]][6,], start0.5[[2]][6,], type="n", xlab="", ylab="", main="5 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][6,], start0.5[[2]][6,], col="gray45", cex=1.5)

plot(start0.25[[1]][6,], start0.25[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][6,], start0.25[[2]][6,], col="gray45", cex=1.5)

plot(start0.1[[1]][6,], start0.1[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][6,], start0.1[[2]][6,], col="gray45", cex=1.5)

plot(start0.5[[1]][11,], start0.5[[2]][11,], type="n", xlab="", ylab="", main="10 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][11,], start0.5[[2]][11,], col="gray45", cex=1.5)

plot(start0.25[[1]][11,], start0.25[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][11,], start0.25[[2]][11,], col="gray45", cex=1.5)

plot(start0.1[[1]][11,], start0.1[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][11,], start0.1[[2]][11,], col="gray45", cex=1.5)

mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()


#colnames(qQ) <- c("ind", "trib", "q", "Q", "classification")

## Rewrite this to classify each individual



for(i in 1:length(qQ[,1])){
    if((qQ$q[i] >= 0.9) && (qQ$Q[i] <= 0.25)){
      qQ$classification[i] <- "ysc"
    }else if((qQ$q[i] <= 0.1) && (qQ$Q[i] <= 0.25)){
      qQ$classification[i] <- "rbt"
    }else if((qQ$q[i] > 0.4) && (qQ$q[i] < 0.6) && (qQ$Q[i] > 0.8)){
      qQ$classification[i] <- "f1"
    }else if((qQ$q[i]) > 0.4 && (qQ$q[i] < 0.6) && (qQ$Q[i] > 0.4) && (qQ$Q[i] < 0.6)){
      qQ$classification[i] <- "f2"
    }else if((qQ$q[i]) >= 0.15 && (qQ$q[i] <= 0.35) && (qQ$Q[i] > 0.4) && (qQ$Q[i] < 0.6)){
      qQ$classification[i] <- "bc1.rbt"
    }else if((qQ$q[i]) >= 0.65 && (qQ$q[i] <= 0.85) && (qQ$Q[i] > 0.4) && (qQ$Q[i] < 0.6)){
      qQ$classification[i] <- "bc1.ysc"
    }else if(((qQ$Q[i]-(2*qQ$q[i])) > -0.1) && ((qQ$Q[i]-(2*qQ$q[i])) < 0.1)){
      qQ$classification[i] <- "bc.rbt"
    }else if(((qQ$Q[i]+(2*qQ$q[i])) > 1.9) && ((qQ$Q[i]+(2*qQ$q[i])) < 2.1)){
      qQ$classification[i] <- "bc.ysc"
    }else{
      qQ$classification[i] <- "other"
    }

}


##########################3
    ## Histogram of q for each generation
    ## Replace with triangle plots
    ## hist(q.allgen[i,], breaks=seq(0,1,0.05),xlab="Prop. Sp. 2 ancestry", main=paste("Gen.", i-1, sep=" "), axes=F, col="gray90")
    ## axis(2)
    ## axis(1, at=c(0,0.5,1))

    ## Q.allgen[i,] <- rep(0, nind.start)
    ## ## Calculate Q
    ## for(j in 1:nind.start){

    ##     if(q.allgen[i,j] == 0 | q.allgen[i,j] == 1){
    ##         Q.allgen[i,j] <- 0 ## Parental species

    ##     }else if(q.allgen[i,j] == 0.5 & (parent1[j] == 1 | parent1[j] == 0) & (parent2[j] == 1 | parent2[j] == 0)){

    ##         Q.allgen[i,j] <- 1

    ##     }else{

    ##         Q.allgen[i,j] <- 5

    ##     }


    ## }

##}

dev.off()

## pdf(paste("hist_",n.generation,"gen_",prop.sp1,"sp1_",nind.start,"ind.pdf",sep=""), width=6.5, height=9)


##par(mfrow=c(ceiling(n.generation/4),4))
##hist(q.allgen[1,], breaks=seq(0,1,0.05), xlab="Prop. Sp. 2 ancestry", main="Starting pop.", axes=F, col="gray90")

##axis(2)
##axis(1, at=c(0,0.5,1))

## Histogram of q for each generation
    ## Replace with triangle plots
    ## hist(q.allgen[i,], breaks=seq(0,1,0.05),xlab="Prop. Sp. 2 ancestry", main=paste("Gen.", i-1, sep=" "), axes=F, col="gray90")
    ## axis(2)
    ## axis(1, at=c(0,0.5,1))

    ## Q.allgen[i,] <- rep(0, nind.start)
