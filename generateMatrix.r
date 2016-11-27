#!/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/bin/Rscript

#Replace this with the user's environment variable
#usage: ./generateMatrix.r dat.eval dat.pca.evec covar kappa
#dat.eval: eigen values 
#dat.pca.evec: engein vectors
#both files from eigenstrat
#covar: covariates to use, example: 1,2,3,5-7
#kappa: optional parameters to tune the randomization, recommend to use 2 or 3. If leaves blank, I will generate matrices for kappa=c(Inf, -4:8, -Inf), you may scan through the files and identify the most optimal one based on genomic inflation factor. The Inf and -Inf are for the controls (Inf=no shuffling, -Inf=ordinary shuffling)

#read the parameters
args <- commandArgs(TRUE)
eval <- scan(args[1])
dat <- read.table(args[2], header=F)
element <- strsplit(strsplit(args[3], ",")[[1]], "-") 

covar <- c()
for (i in  element){
	if(length(i)==1){
		covar <- c(covar, as.numeric(i) )
	}else if(length(i)==2){
		covar <- c(covar, c(as.numeric(i[1]) :as.numeric(i[2]) ))
	}else{
		print("covar specification error")
		stop();
	}
}

print(paste("Using covar #:", paste(covar, collapse=",")))

kappa_list <- as.numeric(args[4])
if(is.na(kappa_list)){
	kappa_list <-  c(Inf, -4:8, -Inf)
}

dat <- dat[,-c(1,2)]
subset <- dat[,covar]
weight<- sqrt(eval[covar]/sum(eval[covar]))
for (i in c(1:length(covar))){
        subset[,i] <- scale(subset[,i])*weight[i]
}
dist <- as.matrix(dist(subset, diag=T))

N <- dim(dist)[1]
for (kappa in kappa_list){
	print(kappa)
	lambda <- 2^kappa
	prob <- exp(-lambda*dist)
	diag(prob) <- 1
	for (i in c(1:N)){
    	prob[,i] <- cumsum(prob[,i])
	}
	prob_formatted <- format(prob, digits=5)
	write.table(prob_formatted, file=paste(sep="", "prob_", kappa, ".txt"), col.names=F, row.names=F, quote=F, sep="\t")
}
