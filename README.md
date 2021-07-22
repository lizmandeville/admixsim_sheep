# admixsim_sheep

This is simulation code developed by Liz Mandeville for a collaborative paper on admixture between translocated populations (Jahner et al. in preparation). Code is modified from Mandeville et al. 2019, Molecular Ecology

The simulations are individual-based, forward-in-time, and not genetically explicit. Ancestry of hybrids is tracked using parameters identical to those produced by "entropy" (Shastry et al. 2021, Gompert et al. 2014), namely q (proportion of ancestry) and Q (interspecific ancestry). 

Changes made for sheep project specifically are:
1) Add in population growth, through the new parameter "growth.rate"
2) Add option to specify differential fitness for ancestry classes. Currently this is achieved through the "sel" parameter, which allows the user to weight probability of an individual being a parent in the next generation, conditional on ancestry of that individual. We assume here that parental species have higher fitness, and that hybrids have lower fitness, quantified by magnitude of "sel".
    Currently the two parental classes have identical fitness, but it would be fairly easy to modify to have either a) one parental species with higher fitness, or 2) hybrids with higher fitness than parentals. 
3) Addition of the "repsim( )" function to easily run and plot replicate runs of the same simulation conditions. This function is simply a wrapper for the main "simulate.hyb( )" function.

Definition of parameters for the simulate.hyb() function:
nind.start = number of individuals at the start of the simulation. (Default=100)

prop.sp1 = proportion of the total number nind.start that come from species 1 (or source
population 1). Note that prop.sp2 is defined within the function as 1-prop.sp1. (Default=0.5)

n.generation = number of generations to run the simulation (Default=10)

makeplot = boolean, whether to output a plot in PDF form (Default=TRUE)
printoutput = boolean, whether to print the output of the function to screen (Default=TRUE)

imm.sp1 = immigrants from species/source 1, as a proportion of nind.start. Note that this will be a constant number of immigrants per generation even with a model that includes population growth. (Default=0)
imm.sp2 = same as imm.sp1, but for species/source 2. (Default=0)

growth.rate = rate of population growth. Values >1 will cause the population to grow by nind[i-1]*growth.rate per generation. (Default=1, no growth)

sel = specification of the fitness reduction in hybrids. This is used as a weight in the prob parameter of sample() within the function to reduce the probability of an individual of hybrid ancestry being selected as a parent in the next generation. (Default=1, bounded 0-1)

repID = a convenenience parameter to allow output files to be labeled differently. Can use any string here. (Default="rep1")

