SKAT_Null_Model = function(formula, data=NULL, out_type="C", n.Resampling=0, type.Resampling="bootstrap", Adjustment=TRUE){
	
	SKAT_MAIN_Check_OutType(out_type)
	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)

	# Check whether n < 1000 and out_type="D", apply the adjustment 
	# if No_Adjustment = FALSE
	if(n< 3000 && out_type=="D" && Adjustment){
		MSG<-sprintf("The sample size is %d, which is < 2000, so the small sample adjustment is applied!\n",n )
		cat(MSG)
		n.Resampling.kurtosis=10000
		#if(n > 1000){
		#	n.Resampling.kurtosis = floor(10000 - (n-1000) * 5)	
		#} 
		#if(n.Resampling.kurtosis < 5000){
		#	n.Resampling.kurtosis = 5000
		#}

		
		re<-SKAT_Null_Model_MomentAdjust(formula, data, n.Resampling, type.Resampling="bootstrap", is_kurtosis_adj=TRUE, n.Resampling.kurtosis=n.Resampling.kurtosis)
		return(re)
	}

	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	if(out_type=="C"){
		re<-Get_SKAT_Residuals.linear(formula, data, n.Resampling, type.Resampling, id_include )
	} else {
		re<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
	}

	class(re)<-"SKAT_NULL_Model"
	return(re)
	
}

Get_SKAT_Residuals.linear = function(formula, data, n.Resampling, type.Resampling, id_include ){

	
 	mod = lm(formula, data=data)
	X1<-model.matrix(formula,data=data)

  	s2 = summary(mod)$sigma**2
  	res = mod$resid
	n1<-length(res)
	res.out<-NULL
	
	if(n.Resampling > 0){

		if(type.Resampling=="permutation"){
			res.out<-res %x% t(rep(1,n.Resampling))
			res.out<-apply(res.out,2,sample)
		} else if(type.Resampling=="bootstrap"){
			res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=sqrt(s2)),ncol=n.Resampling)
			X1_inv<-solve(t(X1) %*% X1)
			res.out<- res.out - (X1 %*% X1_inv) %*% (t(X1) %*% res.out)
		} else if(type.Resampling=="perturbation"){
			res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=1),ncol=n.Resampling)
			res.out<-res.out * res
			stop("Error: Perturbation is no more provided!")
		} else {
			stop("Error: Wrong resampling method!")
		}
	}

  	return(list(res=res, X1=X1,res.out=res.out,out_type="C", 
	n.Resampling=n.Resampling, type.Resampling=type.Resampling,
	id_include=id_include, s2=s2))
}

Get_Resampling_Pvalue<-function (obj, ...){


	ml<-match.call()
	ml.name<-names(ml)
	IDX1<-which(ml.name == "p.value")
	if(length(IDX1) > 0){
		re<-Get_Resampling_Pvalue.numeric(...)
	} else {
		if(class(obj) == "SKAT_OUT"){
			re<-Get_Resampling_Pvalue.SKAT_OUT(obj, ...)
		} else{
			re<-Get_Resampling_Pvalue.numeric(obj, ...)
		}
	}

	return(re)
}

Get_Resampling_Pvalue.SKAT_OUT<-function(obj){

	if(is.null(obj$p.value.resampling)){
		stop("No resampling was applied!")
	}

	n<-length(obj$p.value.resampling)
	n1<-length(which(obj$p.value >= obj$p.value.resampling))
	pval1<-(n1+1)/(n+1)
	
	re<-list(p.value=pval1, is_smaller=FALSE)
	if(n1==0){
		re$is_smaller=TRUE
	}
	
	return(re)
}

Get_Resampling_Pvalue.numeric<-function(p.value,p.value.resampling){

	if(is.null(p.value.resampling)){
		stop("No resampling was applied!")
	}

	n<-length(p.value.resampling)
	n1<-length(which(p.value >= p.value.resampling))
	pval1<-(n1+1)/(n+1)

	re<-list(p.value=pval1, is_smaller=FALSE)
	if(n1==0){
		re$is_smaller=TRUE
	}

	return(re)
}

SKAT_MAIN_Check_OutType<-function(out_type){
 	
	if(out_type != "C" && out_type != "D"){
		stop("Invalid out_type!. Please use either \"C\" for the continous outcome or \"D\" for the dichotomous outcome.")
	}
}

SKAT_Null_Model.boot <- function(formula, data=NULL, out_type="C", n.Resampling=0, type.Resampling="bootstrap", perm.phen, Adjustment=TRUE){
	
	SKAT_MAIN_Check_OutType(out_type)
	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)

	# Check whether n < 1000 and out_type="D", apply the adjustment 
	# if No_Adjustment = FALSE
	if(n< 3000 && out_type=="D" && Adjustment){
		MSG<-sprintf("The sample size is %d, which is < 2000, so the small sample adjustment is applied!\n",n )
		cat(MSG)
		n.Resampling.kurtosis=10000
		#if(n > 1000){
		#	n.Resampling.kurtosis = floor(10000 - (n-1000) * 5)	
		#} 
		#if(n.Resampling.kurtosis < 5000){
		#	n.Resampling.kurtosis = 5000
		#}
		
		re<-SKAT_Null_Model_MomentAdjust.boot(formula, data, n.Resampling, type.Resampling="bootstrap",perm.phen,  is_kurtosis_adj=TRUE, n.Resampling.kurtosis=n.Resampling.kurtosis)
		return(re)
	}

	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	if(out_type=="C"){
		re<-Get_SKAT_Residuals.linear.boot(formula, data, n.Resampling, type.Resampling, id_include, perm.phen)
	} else {
		re<-Get_SKAT_Residuals.logistic.boot (formula, data, n.Resampling, type.Resampling, id_include, perm.phen)
	}

	class(re)<-"SKAT_NULL_Model"
	return(re)
	
}

SKAT_Null_Model_MomentAdjust.boot <- function(formula, data=NULL, n.Resampling=0, type.Resampling="bootstrap", perm.phen, is_kurtosis_adj=TRUE, n.Resampling.kurtosis=10000){
	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-as.numeric(rownames(obj1))

	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	re1<-Get_SKAT_Residuals.logistic.boot (formula, data, n.Resampling, type.Resampling, id_include, perm.phen)
	re2<-NULL

	if(is_kurtosis_adj == TRUE){
		re2<-Get_SKAT_Residuals.logistic.boot(formula, data, n.Resampling, type.Resampling, id_include,  perm.phen)
	}

	re<-list(re1=re1, re2=re2, is_kurtosis_adj= is_kurtosis_adj, type = "binary")

	class(re)<-"SKAT_NULL_Model_ADJ"
	return(re)
	
}

Get_SKAT_Residuals.logistic.boot <- function(formula, data, n.Resampling, type.Resampling,id_include, perm.phen){

 	mod = lm(formula, data)
	X1<-model.matrix(formula,data=data)
	
	glmfit= glm(formula, data=data, family = "binomial")
 	betas = glmfit$coef
  	mu    = glmfit$fitted.values
  	eta   = glmfit$linear.predictors
	n.case = sum(glmfit$y)

	pi_1 = mu*(1-mu)
  	res = glmfit$y- mu
	n1<-length(res)
	res.out<-NULL
	
	mu1<-mu/sum(mu)	# all prob
	res.out<-perm.phen[,1:n.Resampling]

	res.out<-res.out - mu

  	return(list(res=res, X1=X1,res.out=res.out,out_type="D", 
	n.Resampling=n.Resampling, type.Resampling=type.Resampling,
	id_include=id_include, mu=mu,pi_1=pi_1))

}

