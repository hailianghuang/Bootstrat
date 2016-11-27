//////////////////////////////////////////////////////////////////////
//                                                                  //
//  GWiS core module for PLINK (c) 2011 Hailiang Huang, Dan Arking  //
//                                      and Joel S. Bader           //
//                                                                  //
// This file is distributed under the GNU General Public            //
// License, Version 2.  Please see the file COPYING for more        //
// details                                                          //
//                                                                  //
//////////////////////////////////////////////////////////////////////


#include <cmath>
#include <vector>

#include "GWiS_core.h"

SNP_SUMMARY::SNP_SUMMARY(){

}


double & GWiS_math::getLDlower(vector<vector<double> > & LD, int i, int j){
  
  //get the lower triangle of the LD matrix, including the diagonals.  The lower triangle and the diagonals are updated in the regression.

  if(i>j)  return LD[i][j];
  else return LD[j][i];
}


bool GWiS_math::UpdateCovP(int SNP_selected, vector <SNP_SUMMARY> & SNP, vector<vector<double> > & LD){

  //update the covariance between phenotype and the genotype

  vector <SNP_SUMMARY>::iterator snp = SNP.begin();
  double cov=0;
  
  for (;snp != SNP.end(); advance(snp, 1)){

    //skip the selected SNP 
    if(snp->gene_id == SNP_selected)
      continue;

    //skip SNPs that have 0 covariance between the phenotype and genotype (possibly from SNPs already in the model, or SNPs highly correlated to those)
    if(fabs(snp->pheno_cov)<EPS)  continue;
    
    else
      snp->pheno_cov -= getLDlower(LD, snp->gene_id,SNP_selected )*SNP[SNP_selected].pheno_cov/LD[SNP_selected][SNP_selected];
  }
  SNP[SNP_selected].pheno_cov = 0;
}

bool GWiS_math::UpdateCovMat(int SNP_selected, vector<vector<double> > & LD){
  
  //update the LD matrix given the selected SNP
  
  int i, j, nSNP;
  double cov1, cov2;
  
  nSNP=LD.size();
  for (i=0; i<nSNP; i++)
    for (j=0; j<=i; j++){
      if(i == SNP_selected || j == SNP_selected)
	continue;
      getLDlower(LD, i, j)-= getLDlower(LD, i, SNP_selected)*getLDlower(LD, j, SNP_selected)/getLDlower(LD, SNP_selected, SNP_selected);
    }
  
  for (i=0; i<nSNP; i++)
    getLDlower(LD, i, SNP_selected)=0;
}

double GWiS_math::getSSM(SNP_SUMMARY & snp, double var){
  
  //calculate the sum of square of the model
  
  //SNPs that have very small variances are likely to be highly correlated with the SNPs already in the model. 
  if(fabs(var) <EPS)
    return 0;
  else   
    return snp.pheno_cov * snp.pheno_cov / var;
}


int GWiS_math::getBestSNP(vector <SNP_SUMMARY> & SNP, vector<vector<double> > & LD, vector <int> & SNP_selected){
  
  //find the SNP that has the best SSM

  vector <SNP_SUMMARY>::iterator snp = SNP.begin();
  double bestSSM = -1;
  double ssm=0;
  int snp_id=-1;

  for (;snp != SNP.end(); advance(snp, 1)){

    //skip the SNPs that have been added to the model already
    for (int i=0; i< SNP_selected.size(); i++){
      if(snp->gene_id == SNP_selected[i] )
	continue;
    }
    
    //skip the SNPs that have high multiple correlation with SNPs in the model (multiple R2 > 1-VIF_R2), or var=0
    if(snp->var==0 || LD[snp->gene_id][snp->gene_id]/snp->var < VIF_R2 || fabs(LD[snp->gene_id][snp->gene_id]) < EPS)
      continue;

    ssm = getSSM(*snp, LD[snp->gene_id][snp->gene_id]);
    
    if(ssm > bestSSM){
      bestSSM=ssm;
      snp_id = snp->gene_id;
    }
  }
  return snp_id;    
}

double GWiS_math::SSM2BIC(double SSM, double T, int k, int nsample, double RSS){

  //calculate the change in logProb (BIC score) from k-1 to k

  double increment;
  
  increment = log(T-k+1) - log(k);
  increment += log(nsample)/2;
  
  //RSS is the RSS for k-1
  increment += log(RSS - SSM)/2 * nsample;
  increment -= log(RSS)/2 *nsample;
  
  /*
  //Historically had P0=Pr(K=0)
  if(k==1)
    increment +=log(T+1) -log((1-P0)/(P0+ (1-P0)/(1+T)));
  */

  return -increment;

}


double GWiS_math::findBestModel(double TSS, int nSNP, double T,  int nSample, vector <int> &SNP_selected, vector <SNP_SUMMARY> &Z, vector<vector<double> > &LD ){

  //perform Bayes regulated forward regression to find the best SNPs that fit the model

  int k_model=0;
  int bestSNP=-1;
  double delta_BIC=-1;
  double bestSSM =-1;
  
  double BIC = 0;
  
  do{
    bestSNP = GWiS_math::getBestSNP(Z, LD, SNP_selected);
    if(bestSNP<0)
      break;
    bestSSM = GWiS_math::getSSM(Z[bestSNP], LD[bestSNP][bestSNP]);
    delta_BIC = GWiS_math::SSM2BIC(bestSSM, T, ++k_model, nSample, TSS);
    BIC+= delta_BIC;
    if(delta_BIC<0){
      k_model--;
      BIC-= delta_BIC;
      break;
    }
    SNP_selected.push_back(bestSNP);
    TSS-= bestSSM;
    GWiS_math::UpdateCovP(bestSNP, Z, LD);
    GWiS_math::UpdateCovMat(bestSNP, LD);
  }while(k_model < MAX_K && k_model < T);
  
  return BIC;
}

double GWiS_math::calT(int nSNP, vector<vector<double> > &LD ){

  //calculate the effective number of tests

  int k_model=0;
  int bestSNP=0;
  double T=0;
  double bestVar;
    
  //all SNPs have initial weight of 1.
  vector<double> var;
  for (int i =0; i<nSNP; i++)
    var.push_back(1);
  
  vector<double> norm_original;
  for (int i =0; i<nSNP; i++)
    norm_original.push_back(LD[i][i]);
  
  do{

    bestVar = -1;

    //get the weight of the best SNP
    double w =  var[bestSNP];
    T += w;

    for (int i =0; i<nSNP; i++){
      
      //check whether a SNP has no significant weight remaining, or has high multiple correlation with selected SNPs.
      if(fabs(var[i]) <EPS || LD[i][i]/norm_original[i] < VIF_R2){
	var[i]=0;
	continue;
      }
      	
      //adjust the weights of the SNPs 
      double r =getLDlower(LD, i,bestSNP)/sqrt(LD[bestSNP][bestSNP])/sqrt(LD[i][i]);
      var[i] -= r*r*w;
      if(var[i]<0) var[i]=0;
    }

    GWiS_math::UpdateCovMat(bestSNP, LD);
    
    //find out the best SNP for the next iteration
    for(int i=0; i<nSNP; i++){
      if(var[i] > bestVar){
	bestVar=var[i];
	bestSNP=i;
      }
    }
    k_model++;
    
  }while(k_model < nSNP && bestVar > EPS);

  return T;
}
