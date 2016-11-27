//////////////////////////////////////////////////////////////////
//                                                              //
//  GWiS module for PLINK (c) 2011 Hailiang Huang, Dan Arking   //
//                                 and Joel S. Bader            //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "plink.h"
#include "perm.h"
#include "GWiS.h"
#include "GWiS_core.h"
#include "helper.h"

using namespace std;

//extern ofstream LOG;
GWiS::GWiS(Plink & pref){

  //setup the output file
  ofstream ASC;
  int nPheno;
  string f= par::output_file_name + ".GWiS";
  ASC.open(f.c_str(),ios::out);
  ASC.precision(2);

  int nSets = pref.snpset.size();

  //print the headers
  pref.printLOG("Start to run GWiS on "+ int2str(nSets) +" sets\n");
  pref.printLOG("\n");
  pref.printLOG("GWiS v3 is a gene-based test of association developped by the Bader group and the Arking group in the Johns Hopkins University\n");
  pref.printLOG("Reference: Huang H, Chanda P, Alonso A, Bader JS, Arking DE (2011) Gene-Based Tests of Association. PLoS Genetics\n");
  pref.printLOG("\n");

  //calculate the mean and the variance of the phenotype
  double pheno_mean, pheno_var;
  pheno_mean= getPhenoMean(pref.sample, &nPheno);
  pheno_var= getPhenoVar(pref.sample, pheno_mean, &nPheno);
  pref.printLOG(int2str(nPheno)+" individuals, "+"pheno_mean="+ dbl2str(pheno_mean) +", pheno_var=" +  dbl2str(pheno_var) +"\n");

  //prepare for the permutations. 
  pref.pperm->setPermClusters(pref);
  //set the number of tests to be 1 for confidence interval.
  pref.pperm->GWiSTests(1);

  if(par::GWiS_SampleSize == 0)
    pref.printLOG("Sample size will be determined by the SNP pair having the highest missing rate\n");
  if(par::GWiS_SampleSize == 1)
    pref.printLOG("Sample size will be determined by the SNP having the highest missing rate\n");
  if(par::GWiS_SampleSize == 2)
    pref.printLOG("Sample size will be the number of individuals with available phenotypes\n");

  
  ASC << setw(12) << "SET"<< " "
      << setw(6) << "nSNP"<< " "
      << setw(8) << "T"<< " "
      << setw(6) << "K"<< " "
      << setw(8) << "BIC"<< " "
      << setw(12) << "nPerformed"<< " "
      << setw(12) << "nBetter"<< " "
      << setw(12) << "EMP"<< " "
      << setw(50) << "SNP"<< "\n";

  for (int i=0; i<nSets;i++){

    int nSNP = pref.snpset[i].size();
    if(nSNP==0){
      pref.printLOG("Skip running GWiS on empty set "+ int2str(i+1) + " out of " + int2str(nSets) + " total sets\n");
      continue;
    }else{
      pref.printLOG("Running GWiS on set "+ int2str(i+1) + " out of " + int2str(nSets) + " total sets\n");
    }
    //calculate the summary statistics that are needed by GWiS, including LD, genotype variance and covariance between genotype and phenotype
    vector<vector<double> > LD(nSNP, vector<double>(nSNP));
    vector <SNP_SUMMARY> Z(nSNP);
    vector <int> SNP_selected;
    string SNP_selected_name ="";
    pref.printLOG("Done memory allocation\n");
    int nGeno=nPheno;
    pref.printLOG(int2str(nSNP)+ " SNPs\n");

    for (int j=0; j<nSNP; j++){
      int nSize=0;
      pref.printLOG(int2str(j) +"\n");
      CSNP* s=pref.SNP[pref.snpset[i][j]];
      pref.printLOG("SNP acquaired\n");
      Z[j].mean=getGenoMean(s, pref.sample, &nSize);
      pref.printLOG(int2str(j) +"\n");
      if (par::GWiS_SampleSize<2 && nSize < nGeno)
	nGeno=nSize;
      Z[j].var=getGenoVar(s, Z[j].mean, pref.sample, &nSize);
      pref.printLOG(dbl2str(Z[j].var) +"\n");
      Z[j].pheno_cov=getSNPPhenoCov(s,pref.sample , Z[j].mean, pheno_mean, &nSize);
      pref.printLOG(dbl2str(Z[j].pheno_cov) +"\n");
      Z[j].gene_id = j;
    }

    for (int j=0; j<nSNP; j++){
      int nSize=0;
      CSNP* s1=pref.SNP[pref.snpset[i][j]];
      for (int k=0; k<nSNP; k++){
	CSNP* s2=pref.SNP[pref.snpset[i][k]];
	double cov, mean_s1, mean_s2;
	mean_s1 =  getGenoMean(s1, pref.sample, &nSize);
	mean_s2 =  getGenoMean(s2, pref.sample, &nSize);
	cov =  getSNPCov(s1,s2, mean_s1, mean_s2, pref.sample, &nSize);
	LD[j][k] = cov;
	if (par::GWiS_SampleSize ==0 &&  nSize < nGeno)
	  nGeno=nSize;
      }
    }

    //calculate the effective number of tests
    double T=GWiS_math::calT(nSNP, LD);
    pref.printLOG("Number of tests (T) is " + dbl2str(T) + ", total number of SNPs is " + int2str(nSNP) + ", sample size is " + int2str(nGeno) + "\n");

    //the lower triangle, including the diagonals of LD are destroyed after each use.  Restore the LD matrix
    for (int j=0; j<nSNP; j++)
      LD[j][j]=Z[j].var;
    for (int j=0; j<nSNP; j++)
      for (int k=j+1; k<nSNP; k++)
	LD[k][j] = LD[j][k];

    //calculate the original BIC
    double BIC_original = GWiS_math::findBestModel(pheno_var, nSNP, T, nGeno, SNP_selected, Z, LD);
    pref.printLOG(pref.setname[i] + " has k="+ int2str(SNP_selected.size()) + ", BIC="+ dbl2str(BIC_original) + "\n");

    ASC << setw(12) << pref.setname[i]<< " "
	<< setw(6) << int2str(nSNP)<< " "
	<< fixed << setw(8) << T<< " "
	<< setw(6) << SNP_selected.size()<< " "
	<< fixed << setw(8) << BIC_original << " ";
    
    //no need to do permutation if K=0
    if(SNP_selected.size()==0 ){
      ASC<< setw(12) << '0' << " "
	 << setw(12) << '0'<< " "
	 << setw(12) << "-" << " "
	 << setw(50) << "-" << "\n";
      continue;
    }

    for(int j=0; j<SNP_selected.size(); j++)
      SNP_selected_name +=  pref.locus[pref.snpset[i][SNP_selected[j]]]->name +":" ;
    SNP_selected_name.erase(SNP_selected_name.size()-1);

    //skip permutation if --perm is not set
    if(!par::permute ){
      ASC<< setw(12) << '0' << " "
	 << setw(12) << '0'<< " "
	 << setw(12) << "-" << " "
	 << setw(50) << SNP_selected_name << "\n";
      continue;
    }
	   

    bool finished=false;
    while(!finished){
      pref.pperm->permuteInCluster();
      
      //calculate the covariance between phenotype and the genotype
      for (int j=0; j<nSNP; j++){
	int nSize=0;
	CSNP* s=pref.SNP[pref.snpset[i][j]];
	Z[j].pheno_cov=getSNPPhenoCov(s,pref.sample , Z[j].mean, pheno_mean, &nSize ) ;
	LD[j][j]=Z[j].var;
      }
      //Restore the LD matrix
      for (int j=0; j<nSNP; j++)
	for (int k=j+1; k<nSNP; k++)
	  LD[k][j] = LD[j][k];
      
      //calculate the new BIC using shuffled trait
      SNP_selected.clear();
      double BIC = GWiS_math::findBestModel(pheno_var, nSNP, T, nGeno, SNP_selected, Z, LD);
      
      //update the permutation counters
      finished=pref.pperm->updateGWiS(BIC, BIC_original);
    }

    //print info and move to the next set
    pref.printLOG(int2str( pref.pperm->current_reps()) + " performed, " +  int2str(pref.pperm->current_success()) + " success\n");
    ASC<< setw(12) <<  pref.pperm->current_reps() << " "
       << setw(12) <<  pref.pperm->current_success()<< " "
       <<scientific << setw(12) << ((double)pref.pperm->current_success() + 1)  / (pref.pperm->current_reps()+1)  << " "
       << setw(50) <<  SNP_selected_name << "\n";

    pref.pperm->nextSet();

  }
  ASC.close();
  pref.printLOG("Writing GWiS results to [ " + f + " ] \n");

}



double GWiS::getPhenoVar(vector<Individual*> sample, double mean, int* nSize){
  
  //calculate the variance of phenotype
    
  vector<Individual*>::iterator gperson = sample.begin();
  double sum=0;
  *nSize=0;

  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if(!pperson->missing){
      (*nSize)++;
      if(par::qt )
	sum+=pperson->phenotype*pperson->phenotype;
      else if(par::bt)
	if(pperson->aff)
	  sum++;
    }
    gperson++;
  }

  return sum/(*nSize) - mean*mean;

}

double GWiS::getPhenoMean(vector<Individual*> sample, int*nSize){

  //calculate the mean of phenotype

  vector<Individual*>::iterator gperson = sample.begin();
  double sum=0;
  *nSize=0;

  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if(!pperson->missing){
      (*nSize)++;
      if(par::qt )
	sum+=pperson->phenotype;
      else if(par::bt)
	if(pperson->aff)
	  sum++;
    }
    gperson++;
  }

  return sum/(*nSize);
}

double GWiS::getGenoVar(CSNP * s, double mean, vector<Individual*> sample, int* nSize ){

  //calculate the variance of genotype for a particular SNP s
  
  vector<bool>::iterator i1 = s->one.begin();
  vector<bool>::iterator i2 = s->two.begin();
  vector<Individual*>::iterator gperson = sample.begin();

  double sum=0;
  *nSize=0;

  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if (!pperson->missing){
      if(*i2)
	if(*i1){
	  sum+=4;
	  (*nSize)++;
	}
	else{
	  sum+=1;
	  (*nSize)++;
	}
      else
	if(!*i1){
	  //sum+=0;
	  (*nSize)++;
	}
    }
    gperson++;
    i1++;
    i2++;
  }
    
  return sum/(*nSize) - mean*mean;
}

double GWiS::getGenoMean(CSNP *  s, vector<Individual*> sample, int *nSize){

  //calculate the mean of genotype for a particular SNP s
  printf("start\n");
  vector<bool>::iterator i1 = s->one.begin();
  vector<bool>::iterator i2 = s->two.begin();
  vector<Individual*>::iterator gperson = sample.begin();
  
  double sum=0;
  *nSize=0;
  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if (!pperson->missing){
      printf("%d, %g\n", *nSize, sum);
      if(*i2)
	if(*i1){
	  sum+=2;
	  (*nSize)++;
	}else{
	  sum+=1;
	  (*nSize)++;
	}
      else
	if(!*i1){
	  //sum+=0;
	  (*nSize)++;
	}
    }
    gperson++;
    i1++;
    i2++;
  }
  
  return sum/(*nSize);
}

double GWiS::getSNPCov(CSNP* s1 , CSNP* s2, double mean_s1, double mean_s2,  vector<Individual*> sample, int* nSize){

  //calculate the covariance between SNPs s1 and s2

  vector<bool>::iterator i11 = s1->one.begin();
  vector<bool>::iterator i12 = s1->two.begin();
  vector<bool>::iterator i21 = s2->one.begin();
  vector<bool>::iterator i22 = s2->two.begin();
  vector<Individual*>::iterator gperson = sample.begin();

  double sum=0;
  *nSize=0;

  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if (!pperson->missing){
      //s1=11 and s2=11 -> s1s2=4
      if(*i11 && *i12 && *i21 && *i22){
	sum+=4;
	(*nSize)++;
      }
      //s1=10 or s2=10 -> missing
      else if((*i11 && !*i12 )||(*i21 && !*i22)){
	//missing
	//error("No missing genotype is allowed for GWiS");
      }
      //s1=00 or s2=00 -> s1s2=0
      else if((!*i11 && !*i12 )||(!*i21 && !*i22)){
	sum+=0;
	(*nSize)++;
      }
      //s1=01 and s2=01 -> s1s2=1
      else if((!*i11 && *i12 )&&(!*i21 && *i22)){
	sum+=1;
	(*nSize)++;
      }
      //(s1=01 and s2=11) or (s1=11 and s2=01) -> s1s2=2
      else {
	sum+=2;
	(*nSize)++;
      }
    }
    i11++;
    i12++;
    i21++;
    i22++;
    gperson++;
  }
  return sum/(*nSize) - mean_s1 *mean_s2;  
}

double GWiS::getSNPPhenoCov(CSNP * s, vector<Individual*> sample, double mean_s, double mean_pheno, int*nSize){
  
  //calculate the covariance between SNP s and phenotype

  vector<Individual*>::iterator gperson = sample.begin();
  vector<bool>::iterator i1 = s->one.begin();
  vector<bool>::iterator i2 = s->two.begin();

  double sum=0;
  *nSize=0;
  double pheno=0;

  while ( gperson !=  sample.end() ){
    Individual * pperson = (*gperson)->pperson;
    if (!pperson->missing){
      if(par::qt )
	pheno=pperson->phenotype;
      else if(par::bt){
	if(pperson->aff)
	  pheno=1;
	else
	  pheno=0;
      }
      if(!*i1 && !*i2){
	//sum+=0;
	(*nSize)++;
      }
      else if(*i1 && *i2){
	sum+=2*pheno;
	(*nSize)++;
      }
      else if(!*i1 && *i2){
	sum+=pheno;
	(*nSize)++;
      }
    }
    i1++;
    i2++;
    gperson++;
  }
  return sum/(*nSize) - mean_s*mean_pheno;
}
