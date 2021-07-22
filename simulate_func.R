## The functions in this script allow forward-in-time, individual-based simulations of hybridization/admixture outcomes. Ancestry of hybrids is tracked using parameters identical to those produced by "entropy" (Shastry et al. 2021, Gompert et al. 2014), namely q (proportion of ancestry) and Q (interspecific ancestry). See README for more info and parameter definitions.

## Set parameters manually for testing - holdover from development, delete this block eventually
## nind.start <- 10
## prop.sp1 <- 0.5
## n.generation <- 10
## makeplot <- TRUE
## printoutput <- TRUE
## imm.sp1 <- 0.1
## imm.sp2 <- 0
## growth.rate <- 1.1
## sel <- 0.2
## repID <- "rep1"

simulate.hyb <- function(nind.start = 100, prop.sp1=0.5, n.generation=10, makeplot=TRUE, printoutput=TRUE, imm.sp1=0, imm.sp2=0, growth.rate=1, sel=1, repID="rep1"){

    ## Nind in each generation - set up empty vector, initialized at nind.start
    nind <- numeric(length=(n.generation+1))
    nind[1] <- nind.start

    ## Number of immigrants in each generation, as a proportion of nind.start. Default is 0.
    ## Currently a constant rate of migrants, but wouldn't be hard to change
    nind.imm <- round(nind.start * imm.sp1) + round(nind.start * imm.sp2)
    nind.imm.sp1 <- round(nind.start * imm.sp1)
    nind.imm.sp2 <- round(nind.start * imm.sp2)

    ## Allow growth over time if growth.rate > 1
    ## multiply the nind in the previous generation by the growth.rate, round so as not to get fractional individuals

    for(i in 2:(n.generation+1)){
        nind[i] <- round(nind[i-1]*growth.rate, digits=0) + nind.imm
    }

    ## Set up matrix for keeping track of q, or proportion of ancestry. Note that in this simulation, q is calculated from Q values for each individual
    ## Matrices are n.generation+1 rows and nind[length(nind)] rows (the max number of individuals; need to do this to keep matrix rectangular, will have NA in earlier generations in a model with population growth
    qfromQ.allgen <-matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    qfromQ.allgen[1,1:nind.start] <- c(rep(0, round(prop.sp1*nind.start)), rep(1, (nind.start-round(prop.sp1*nind.start))))

    ## Set up empty matrices for keeping track of Q in all 4 categories
    Q11.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    Q11.allgen[1,1:nind.start] <- c(rep(1,round(prop.sp1*nind.start, digits=0)), rep(0,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Q12.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    Q12.allgen[1,1:nind.start] <- c(rep(0, nind.start))

    Q21.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    Q21.allgen[1,1:nind.start] <- c(rep(0, nind.start))

    Q22.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    Q22.allgen[1,1:nind.start] <- c(rep(0,round(prop.sp1*nind.start, digits=0)), rep(1,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Qinter.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind[length(nind)])
    Qinter.allgen[1,1:nind.start] <- rep(0, nind.start)

    if(makeplot == TRUE){
        pdf(paste("qQ_",n.generation,"gen_",prop.sp1,"sp1_",nind.start,"ind_", nind.imm.sp1, "immsp1_", nind.imm.sp2, "immsp2_", growth.rate, "growth.rate_",sel,"sel_",repID,".pdf",sep=""))
        par(mfrow=c(ceiling(n.generation/4),4))
    }

    for(i in 2:(n.generation+1)){

    ## sample a vector of integers 1:nind, then use that vector to index q.allgen and Q.allgen so we can get matching q and Q
        ## sampling WITH REPLACEMENT, i.e. the fact that an individual has already produced one offspring does not make it ineligible to produce another


        if(sel==1){
            ## Option 1, all individuals from prev. generation are equally likely to produce offspring
            randvec1 <- sample(1:(nind[i-1]+nind.imm), nind[i], replace=T)
            randvec2 <- sample(1:(nind[i-1]+nind.imm), nind[i], replace=T)
        }else{
            ## Option2, weight parental individuals as higher fitness
            ## Note, could adapt to weight parental species similar, but for now they are the same and have higher fitness than hybrids

            ## Set hybrid weight at "sel", a probability weight of being selected as a parent in the next gen
            selectvec <- rep(sel, (nind[i-1]+nind.imm))
            ## replace elements of selection vector with 1 for parental species individuals, indicating greater probability of parental species individuals from prev generation producing offspring in this generation
            selectvec[qfromQ.allgen[i-1,1:nind[i-1]]==0 | qfromQ.allgen[i-1,1:nind[i-1]]==1] <- 1

            randvec1 <- sample(1:(nind[i-1]+nind.imm), nind[i], replace=T, prob=selectvec)
            randvec2 <- sample(1:(nind[i-1]+nind.imm), nind[i], replace=T, prob=selectvec)
        }

        ## sample from prev. generation for parents of new offspring, assume random mating if sel=1, otherwise weighted probability of sampling as above
        ## use random vectors to index both q and Q
        ## Updated Nov. 5, 2018 to include immigration

        parent1.Q11 <- c(Q11.allgen[i-1,1:nind[i-1]], rep(1,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q12 <- c(Q12.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q21 <- c(Q21.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q22 <- c(Q22.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(1,nind.imm.sp2))[randvec1]

        parent2.Q11 <- c(Q11.allgen[i-1,1:nind[i-1]], rep(1,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q12 <- c(Q12.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q21 <- c(Q21.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q22 <- c(Q22.allgen[i-1,1:nind[i-1]], rep(0,nind.imm.sp1), rep(1,nind.imm.sp2))[randvec2]

        Q11.allgen[i,1:nind[i]] <- (parent1.Q11+parent1.Q12)*(parent2.Q11+parent2.Q12)
        Q12.allgen[i,1:nind[i]] <- (((parent1.Q11+parent1.Q12)*(parent2.Q21+parent2.Q22))+((parent1.Q21+parent1.Q22)*(parent2.Q11+parent2.Q12)))/2 ## Need to do this to retain symmetry of Q12,Q21
        Q21.allgen[i,1:nind[i]] <- Q12.allgen[i,1:nind[i]]
        Q22.allgen[i,1:nind[i]] <- (parent1.Q21+parent1.Q22)*(parent2.Q21+parent2.Q22)

        qfromQ.allgen[i,1:nind[i]] <- Q11.allgen[i,1:nind[i]]+((Q12.allgen[i,1:nind[i]]+Q21.allgen[i,1:nind[i]])/2)

        Qinter.allgen[i,1:nind[i]] <- Q12.allgen[i,1:nind[i]]+Q21.allgen[i,1:nind[i]]

        if(makeplot == TRUE){
            plot(qfromQ.allgen[i,1:nind[i]], Qinter.allgen[i,1:nind[i]], main=paste("Gen.", i-1, "nind =", nind[i], sep=" "), xlab="q", ylab="Q", xlim=c(0,1), ylim=c(0,1), type="n", axes=F)
            axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
            axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
            arrows(0,0,0.5,1, length=0, col="gray90")
            arrows(0.5,1,1,0, length=0, col="gray90")
            points(qfromQ.allgen[i,1:nind[i]], Qinter.allgen[i,1:nind[i]], col="gray45")
        }
    }
    if(makeplot == TRUE){
        dev.off()

    }

    if(printoutput == TRUE){
        ## Output of function is list object with rectangular matrices containing q as [[1]] and Q as [[2]]
        return(list(qfromQ.allgen, Qinter.allgen))
    }
}


repsim <- function(nrep=3, nind.start = 100, prop.sp1=0.5, n.generation=10, makeplot=FALSE, printoutput=TRUE, imm.sp1=0, imm.sp2=0, growth.rate=1, sel=1){

    ## Plot just the final generation for each rep

    pdf(paste("repsim",nrep,"_",n.generation,"gen_",prop.sp1,"sp1_",nind.start,"ind_", imm.sp1, "imm.sp1_", imm.sp2, "immsp2_", growth.rate, "growth.rate_",sel,"sel.pdf",sep=""), width=8, height=15)
    par(mfrow=c(5,2), mar=c(2,2,1,1))

    for(j in 1:nrep){
        tmpsim <- simulate.hyb(nind.start=nind.start, prop.sp1=prop.sp1, n.generation=n.generation, imm.sp1=imm.sp1, imm.sp2=imm.sp2, growth.rate=growth.rate, sel=sel, repID=paste("rep",j,sep=""), makeplot=makeplot, printoutput=printoutput)

        plot(tmpsim[[1]][n.generation+1,], tmpsim[[2]][n.generation+1,], type="n", xlab="", ylab="", main=paste("Rep", j, ", ", n.generation," gen.",sep=""), axes=F, xlim=c(0,1), ylim=c(0,1))
        axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
        axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
        arrows(0,0,0.5,1, length=0, col="gray90")
        arrows(0.5,1,1,0, length=0, col="gray90")
        points(tmpsim[[1]][n.generation+1,], tmpsim[[2]][n.generation+1,], col="gray45", cex=1.5)

    }

    dev.off()

    #return(tmpsim)


}
