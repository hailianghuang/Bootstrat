
//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>

#include "linear.h"
#include "helper.h"
#include "plink.h"
#include "options.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"


// Usage of Model

// Fit model                       LinearModel lm;
// Give pointer to PLINK           lm.setPlink(this);
// Set missing data                lm.setMissing();
// Set dependent (adds intercept)  lm.setDependent(Y);
// Addive effects, labels          lm.addAdditiveSNP(CSNP*); 
//                                 lm.label.push_back("ADD");
//	                           lm.addDominanceSNP(CSNP*);
// Covariates?                     lm.addCovariate(int);
//              	           lm.label.push_back("COV"+int2str(c+1));
// Interactions?       	           lm.addInteraction(int,int);
// Build design matrix             lm.buildDesignMatrix();
// Prune out any missing           lm.pruneY();
// Fit logistic model              lm.fitLM();


vector_t Plink::linearAssoc(bool print_results, Perm & perm)
{

  vector<double> results(nl_all);

  ofstream ASC;
  if (print_results)  
    {
      string f = par::output_file_name + ".assoc.linear";
      //if (par::twoDFmodel) f += ".genotypic";
      printLOG("Writing linear model association results to [ " + f + " ] \n");
      ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " " 
	    << setw(par::pp_maxsnp) << "SNP" << " " 
	    << setw(10) << "TEST" << " "
	    << setw(8) << "NMISS" << " " 
	    << setw(10) << "BETA" << " " 
	  //	    << setw(10) << "SE" << " " 
	    << setw(12) << "CHISQ" << " " 
	    << setw(12) << "P" << " " 
	    << "\n";
	ASC.precision(4);
    }
 

  ////////////////////////////
  // Iterate over each locus
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // Skip possibly
      if (par::adaptive_perm && !perm.test[l])
	{
	  l++;
	  s++;
	  continue;
	}
      
      // Haploid?
      if ( par::chr_sex[locus[l]->chr] ||
	   par::chr_haploid[locus[l]->chr] )
	{
	  l++;
	  s++;
	  continue;
	}
      
      // Fit model (sets dependent to pheno)
      LinearModel lm(this);
      
      // Set missing data
      lm.setMissing();
      
      // Main effect of SNP
      lm.addAdditiveSNP(*s); 
      lm.label.push_back("ADD");

      // Dominance?
      if (par::twoDFmodel)
	{
	  lm.addDominanceSNP(*s);
	  lm.label.push_back("DOM");
	}
      // Covariates?
      if (par::clist)
	{
	  for (int c=0; c<par::clist_number; c++)
	    {
	      lm.addCovariate(c);
	      lm.label.push_back("COV"+int2str(c+1));
	    }
	}
      
      // 0 intercept
      // 1 A
      // 2 D
      // 3 cov1
      // 4 cov2

      
      // Basic SNP x covariate interaction? 
      if (par::simple_interaction && par::clist)
	{
	  if ( !par::twoDFmodel )
	    for (int c=0; c<par::clist_number; c++)
	      {
		lm.addInteraction(1,c+2);
		lm.label.push_back("ADDxCOV"+int2str(c+1));	  
	      }
	  else 
	    for (int c=0; c<par::clist_number; c++)
	      {
		lm.addInteraction(1,c+3);
		lm.label.push_back("ADDxCOV"+int2str(c+1));	  
		
		lm.addInteraction(2,c+3);
		lm.label.push_back("DOMxCOV"+int2str(c+1));	  
	      }
	}

      // Build design matrix
      lm.buildDesignMatrix();

      // Prune out any remaining missing individuals
      lm.pruneY();
      
      // Fit logistic model (Newton-Raphson)
      lm.fitLM();

      // Check for multi-collinearity
      lm.validParameters();
      
      // Obtain estimates and statistic
      if (print_results)
	lm.displayResults(ASC,locus[l]);
      
      if (par::twoDFmodel)
 	{
 	  vector_t h;
 	  h.resize(1,0);
 	  matrix_t H;
 	  sizeMatrix(H,2,3);

	  // Joint test of dominant and recessive
 	  H[0][0] = 0; 	  H[0][1] = 1;   H[0][2] = 0; 
	  H[1][0] = 0; 	  H[1][1] = 0;   H[1][2] = 1; 

	  double chisq = lm.isValid() ? lm.linearHypothesis(H,h) : 0;
	  double pvalue = chiprobP(chisq,2);

	  ASC << setw(4) << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 
	      << setw(10) << "GENO_2DF" << " "
	      << setw(8) << lm.Ysize() << " " 
	      << setw(10) << "NA" << " ";

	  if (lm.isValid() && realnum(chisq) )
	    ASC << setw(12) << chisq << " " 
		<< setw(12) << pvalue << "\n"; 
	  else
	    ASC << setw(12) << "NA" << " " 
		<< setw(12) << "NA" << "\n"; 
	    
 	}
      
      // Store statistic
      results[l] = lm.getStatistic();
     
      // Next SNP
      s++;
      l++;

    }
  
  if (print_results)
    ASC.close();
  
  return results;

}

