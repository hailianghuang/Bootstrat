#!/usr/local/bin/Rscript
##Load in regular SKAT
library("SKAT")
##functions for running skat with PIP
source("./SKAT-src/Main.R")
source("./SKAT-src/Null_Model.R")
source("./SKAT-src/KMTest_Logistic_VarMatching.R")
source("./SKAT-src/KMTest_Optimal_VarMatching.R")
source("./SKAT-src/Function.R")
source("./SKAT-src/skat.functions.for.mod.perm.R")

#set to TRUE if running interactively
INTERACTIVE <- FALSE

args <- commandArgs(TRUE)
pheno.ID <- args[1]
ID.shuffle <- args[2]
geno_file <- args[3]

if(INTERACTIVE){
	pheno.ID <- "pheno.ID.txt"
	ID.shuffle <- "ID.shuffle.txt"
	geno_file <- "wgas1.geneA.raw"
}

##linker between the bootstrapped samples and the observed data
numSAM.truSAM <- read.table(pheno.ID)
names(numSAM.truSAM) <- c("FID", "IID", "PHENOTYPE", "ID")
numSAM.truSAM <- numSAM.truSAM[, -3]

##call in set bootstrapped samples from plink
shuffle <- read.table(ID.shuffle)
shuffle <- t(shuffle)
resample.n <- dim(shuffle)[2] ##number of resampling

##call in pheno

##matrix of genotypes want to aggregate
pheno.geno <- read.table(geno_file,header=T,as.is=T)
pheno.geno$PHENOTYPE[with(pheno.geno, !is.na(PHENOTYPE) & PHENOTYPE!=1 & PHENOTYPE!=2) ] <- NA
pheno.geno$PHENOTYPE <- pheno.geno$PHENOTYPE -1

###check stopped here
##merge data with the bootstrapped samples and put in correct order
dat <- merge(numSAM.truSAM, pheno.geno, by=c("FID","IID"))
dat <- dat[order(dat$ID),]

perm <- cbind(numSAM.truSAM,  apply(shuffle, 2 ,function(x, dat){dat$PHENOTYPE[x] }, dat))

##genotype matrix
gene.geno <- as.matrix(dat[,-c(1:7)])
pheno.shuffled <- perm[, -c(1:3)]

##weights to use
my.maf <- apply(gene.geno,2,mean,na.rm=T)
my.weights <- 1/sqrt(my.maf*(1-my.maf))
my.weights <- my.weights[which(!(my.weights == Inf))]

##SKAT
set.seed(1)
obj.noadj <- SKAT_Null_Model(PHENOTYPE~1,data=dat,out_type="D",n.Resampling=resample.n)
skat.noadj <- SKAT(gene.geno,obj.noadj,weights=my.weights,kernel="linear.weighted")
print(paste("SKAT P:", format(skat.noadj$p.value, digits=3)))
print(paste("SKAT-resampling P:", format(Get_Resampling_Pvalue(skat.noadj)$p.value, digits=3)))

##Bootstrat
obj.boot <- SKAT_Null_Model.boot(PHENOTYPE~1,data=dat,out_type="D",n.Resampling=resample.n, perm.phen=pheno.shuffled)
skat.boot <- SKAT(gene.geno, obj.boot, weights=my.weights,kernel="linear.weighted")
print(paste("SKAT-BOOTSTRAT P:", format(Get_Resampling_Pvalue(skat.boot)$p.value, digits=3) ))


