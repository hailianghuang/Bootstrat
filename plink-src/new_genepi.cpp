

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
#include <fstream>
#include <algorithm>
#include <functional>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "sets.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"

void calcGENEPIMeanVariance(vector<CSNP*> &,
			    int, int,
			    bool,
			    vector<double> &,
			    vector<vector<double> > &,
			    vector<Individual*>&,
			    vector<int> &,
			    vector<int> &);

vector<double> CCA(bool perm, 
		   vector<vector<int> > & blperm,
		   vector<vector<int> > & blperm_control,
		   Set & S,
		   Plink & P);


matrix_t t(matrix_t a)
{
  matrix_t r = a;
  return a;
}

void Plink::driverSCREEPI()
{

  ///////////////////////////////
  // Gene-based epistasis
  
  printLOG("\n\n*** Warning: method under development, not for general use *** \n\n");
  

  //////////////////////////////////////////
  // Requires that sets have been speciefied

  if (par::set_test) readSet();
  else error("Need to specify genes with --set {filename} when using --genepi\n");

    
  //////////////////
  // SET statistics

  Set S(snpset);


  //////////////////////////////////////////////
  // Prune SET (0-sized sets, MAF==0 SNPs, etc) 

  S.pruneSets(*this);
  
  int ns = snpset.size();
  if (ns < 2)
    error("Need to specify at least two fully valid sets\n");


  int ncase = 0;
  int ncontrol = 0;

  /////////////////////////////////////////////////////////
  // Prune based on VIF, jointly for cases and controls

  printLOG("\nConsidering cases and controls: ");
  setFlags(false);
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( ! (*person)->missing )
	{
	  (*person)->flag = true;
	  ncase++;
	}
      person++;
    }

  string original_outfile = par::output_file_name;
  par::output_file_name += ".all";
  S.pruneMC(*this,false,par::vif_threshold);
  par::output_file_name = original_outfile;


  // Write finalized set
  ofstream SET1, SET2;
  string f = par::output_file_name + ".all.set.in";
  
  printLOG("Writing combined pruned-in set file to [ " + f + " ]\n");
  SET1.open(f.c_str(),ios::out);

  f = par::output_file_name + ".all.set.out";
  printLOG("Writing combined pruned-out set file to [ " + f + " ]\n");
  SET2.open(f.c_str(),ios::out);
  
  for (int s=0; s<snpset.size(); s++)
    {
      
      int nss = snpset[s].size();
      
      SET1 << setname[s] << "\n";
      SET2 << setname[s] << "\n";
      
      for (int j=0; j<nss; j++)
	{
	  if (S.cur[s][j])
	    SET1 << locus[snpset[s][j]]->name << "\n";
	  else
	    SET2 << locus[snpset[s][j]]->name << "\n";
	}
      
      SET1 << "END\n\n";
      SET2 << "END\n\n";
    }
  
  SET1.close();
  SET2.close();
  

  // Prune empty sets once more:

  S.pruneSets(*this);
  
  ns = snpset.size();
  if (ns < 2)
    error("Need to specify at least two fully valid sets\n");



  ////////////////////////////////
  // Set up permutation structure

  // Specialized (i.e. cannot use Perm class) as this 
  // requires a block-locus permutation

  // First block is fixed
  
  vector<vector<int> > blperm(ns);
  vector<vector<int> > blperm_control(ns);

  for (int i=0; i<ns; i++)
    {
      // A slot for each individual per locus
      for (int j=0; j<n; j++)
	if ( sample[j]->aff && ! sample[j]->missing )
	  blperm[i].push_back(j);
      
      for (int j=0; j<n; j++)
	if ( (!sample[j]->aff) && ! sample[j]->missing )
	  blperm_control[i].push_back(j);
      
    }




  /////////////////
  // Original data

  vector<double> original = CCA(false, 
				blperm, 
				blperm_control,
				S,*this);
  


  if (!par::permute) 
   return;
  
}


// MY MAJOR CHANGES START HERE


///////////////////////////////////////////////////////////
// First CCA function: use for case-control logit analysis

vector<double> CCA_logit(bool perm, 
		         vector<vector<int> > & blperm,
		         vector<vector<int> > & blperm_control,
		         Set & S,
		         Plink & P)
  
{


  // For each pair of sets (genes), this function will do the following:
  
  // 
  // Step 1.  Construct a p+q by p+q SNP covariance matrix
  // (sigma). Cases and controls together. Partition into
  // sigma11,sigma12,sigma21,sigma22.

  // Step 2.  Calculate the p x p matrix M = inv(sqrt(sig11)) %*%
  // sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11)).

  // Step 3.  Determine the p eigenvalues of M. The sqrt(eigen(M)) = p
  // canonical correlations

  // Step 4a. Calculate the p x p eigenvectors of M (e). These are
  // required to compute the coefficients used to build the p
  // canonical variates a[k] for gene1 (see below)

  // Step 4b. Calculate the q x q eigenvectors of inv(sqrt(sig22)) %*%
  // sig21 %*% inv(sig11) %*% sig12 %*% inv(sqrt(sig22)) (f). These
  // are required to compute the coefficients used to build the p
  // canonical variates a[k] for gene2 (see below). The first p of the
  // f eigenvectors se can be calculated as: f[k] = (1/eigen[k]) *
  // inv(sqrt(sig22)) %*% sig21 %*% inv(sqrt(sig11)) %*% e[k]

  // Step 5.  Calculate the gene1 (a[k], pxp) and gene2 (b[k], pxq)
  // coefficients used to create the canonical variates associated
  // with the p canonical correlations a[k] = t(e[k]) %*%
  // inv(sqrt(sig11)) b[k] = t(f[k]) %*% inv(sqrt(sig22)), where k =
  // 1,2,...,p

  // Step 6.  Compute the gene1 (V[1]) and gene2 (W[1]) canonical
  // variates associated with the highest canonical correlation (ie
  // k=1) NOTE: the original variables of data need to be standardised
  // first! Otherwise, the resulting correlation between variate.gene1
  // and variate.gene1 != estimated cancor.  V[1] = data %*% a[1] W[1]
  // = data %*% b[1]

  // Step 7.  Fit a logistic regression model with V and W as main
  // effects + an interaction term : y ~ V + W + V*W




  ///////////////
  // Output file

  ofstream EPI;

  if (!perm)
  {  
   string f = par::output_file_name+".genepi";
   P.printLOG("\nWriting gene-based epistasis tests to [ " + f + " ]\n");
   EPI.open(f.c_str(), ios::out);
   EPI.precision(4);
  }
  

  //////////////////////////////////
  // Canonical correlation analysis

  vector<double> res;

  int ns = P.snpset.size();
  
  // Consider each pair of genes
  
  for (int s1=0; s1 < ns-1; s1++)
      for (int s2 = s1+1; s2 < ns; s2++)
      {
	  

 	  ////////////////////////////////////////////////////////
	  // Step 0. Do we want to recode the SNPs into additive coding?			// *********



 	  ////////////////////////////////////////////////////////
	  // Step 1. Construct covariance matrix (cases and controls together)
          //    And partition covariance matrix:
	  //    S_11  S_21
	  //    S_12  S_22
	  
	  int n1=0, n2=0;
	  
	  vector<vector<double> > sigma(0);
	  vector<double> mean(0);
	  vector<CSNP*> pSNP(0);
	  
	  /////////////////////////////
	  // List of SNPs for both loci
	  
	  for (int l=0; l<P.snpset[s1].size(); l++)
	      if ( S.cur[s1][l] )
	      {
		  pSNP.push_back( P.SNP[ P.snpset[s1][l] ] );
		  n1++;
	      }
	  for (int l=0; l<P.snpset[s2].size(); l++)
	      if ( S.cur[s2][l] )
	      {
		  pSNP.push_back( P.SNP[ P.snpset[s2][l] ] );
		  n2++;
	      }
	  

          // NOTE: we need to make sure that n1 < n2. Migth cause problems below if this is not the case.	// *********
	  int n12 = n1 + n2;
	  int ne = n1 < n2 ? n1 : n2;		// ne = min(p,q)

  
	  ///////////////////////////////////
	  // Construct covariance matrix (cases and controls together)
	  
	  P.setFlags(false);
	  vector<Individual*>::iterator person = P.sample.begin();
	  while ( person != P.sample.end() )
	  {
	      // if ( (*person)->aff && !(*person)->missing )                       // MARF: this was choosing cases only, right?
              if ( !(*person)->missing )                      
		  (*person)->flag = true;
	      person++;
	  }
	  
	  calcGENEPIMeanVariance(pSNP, 
				 n1,n2,
				 false,
				 mean, 
				 sigma, 
				 P.sample , 
				 blperm[s1], 
				 blperm[s2] );
	  
	  
	  ///////////////////////////
	  // Partition covariance matrix
	  
	  vector<vector<double> > I11;
	  vector<vector<double> > I12;
	  vector<vector<double> > I21;
	  vector<vector<double> > I22;
	  vector<vector<double> > inv_sqrt_I22;
	  
	  sizeMatrix( I11, n1, n1);
	  sizeMatrix( I12, n1, n2);
	  sizeMatrix( I21, n2, n1);
	  sizeMatrix( I22, n2, n2);
	  sizeMatrix( inv_sqrt_I22, n2, n2);
	  
	  for (int i=0; i<n1; i++)
	      for (int j=0; j<n1; j++)
		  I11[i][j] = sigma[i][j];
	  
	  for (int i=0; i<n1; i++)
	      for (int j=0; j<n2; j++)
		  I12[i][j] = sigma[i][n1+j];
	  
	  for (int i=0; i<n2; i++)
	      for (int j=0; j<n1; j++)
		  I21[i][j] = sigma[n1+i][j];
	  
	  for (int i=0; i<n2; i++)
	      for (int j=0; j<n2; j++)
		  I22[i][j] = sigma[n1+i][n1+j];
	  
	  
 	  ////////////////////////////////////////////////////////
	  // Step 2. Calculate the p x p matrix M1 = inv(sqrt(sig11))
	  // %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))
       	  
	  I11 = msqrt(I11);
	  I11 = svd_inverse(I11);
	  inv_sqrt_I22 = msqrt(I22);			// For Step 4b
	  inv_sqrt_I22 = svd_inverse(inv_sqrt_I22);
	  I22 = svd_inverse(I22);
	  
	  matrix_t tmp;
	  matrix_t M1;
	  
	  multMatrix(I11, I12, tmp);
	  multMatrix(tmp, I22, M1);
	  multMatrix(M1, I21, tmp);
	  multMatrix(tmp, I11, M1);


 	  ////////////////////////////////////////////////////////
	  // Step 3. Determine the p eigenvalues of M1. The
	  // sqrt(eigen(M)) = p canonical correlations
	  
	  vector_t gene1_eigvalues = eigenvalues(M1);
	  
	  // Display eigen values
	  cout << "E (1) values\n";
	  display(gene1_eigvalues);
	  cout << "\n";

	  matrix_t MT;
	  sizeMatrix(MT,3,3);

	  MT[0][0] =  1.11067269;  MT[0][1] =-0.1104839; MT[0][2] = 0.07441985;
	  MT[1][0] =  -0.11048390; MT[1][1] = 0.9777328; MT[1][2] = -0.14673707;
	  MT[2][0] =  0.07441985;  MT[2][1] =-0.1467371; MT[2][2] = 0.88281471;

	  Eigen gene1_eigen = eigenvectors(M1);
	  
	  // Display eigen values - check with above
	  cout << "E values\n";
	  display(gene1_eigen.d);
	  cout << "\n";
	  
	  // Display and store eigen vectors
	  cout << "E vectors\n";
	  display(gene1_eigen.z);
	  cout << "\n\n";

	  exit(0);



 	  ////////////////////////////////////////////////////////
          // Step 4a. Calculate the p x p eigenvectors of M (e). These
          // are required to compute the coefficients used to build
          // the p canonical variates a[k] for gene1 (see below)


	  Eigen gene1_eigen = eigenvectors(M1);
	  
	  // Display eigen values - check with above
	  cout << "E values\n";
	  display(gene1_eigen.d);
	  cout << "\n";
	  
	  // Display and store eigen vectors
	  cout << "E vectors\n";
	  display(gene1_eigen.z);
	  cout << "\n\n";
          matrix_t gene1_eigvectors = gene1_eigen.z;		// *********
	  
  
	  sort(gene1_eigvalues.begin(),gene1_eigvalues.end(),greater<double>());
	  
	  for (int i=0; i<ne; i++)
	  {
	      if ( gene1_eigvalues[i] < 0 ) 
		  gene1_eigvalues[i] = 0;
	      else if ( gene1_eigvalues[i] >= 1 )
		  gene1_eigvalues[i] = 0.999;
	  }

	  // Display canonical correlations, sorted
	  cout << "All canonical correlations \n";
	  for (int i=0; i<ne; i++)
	    cout << "Cancor " << i << " = " << sqrt(gene1_eigenvalues[i]) << "\n";
	  cout << "\n";
	  

 	  ////////////////////////////////////////////////////////
          // Step 4b. Calculate the q x q eigenvectors of M2
          // (f). These are required to compute the coefficients used
          // to build the p canonical variates b[k] for gene2 (see
          // below). The first p are given by: 

	  // f[k] = (1/eigen[k]) * inv_sqrt_I22 %*% I21 %*% inv_sqrt_sig11 %*% e[k] 

	  // for (k in 1:p) { e.vectors.gene2[,k] =
          // (1/sqrt(e.values[k])) * inv.sqrt.sig22 %*% sig21 %*%
          // inv.sqrt.sig11 %*% e.vectors.gene1[,k] }
         
          matrix_t gene2_eigvectors;		// *********
	  sizeMatrix( gene2_eigvectors, n1, n1);		

	  matrix_t tmp2;	  
	  multMatrix(inv_sqrt_I22, I21, tmp);
	  multMatrix(tmp, I11, tmp2);

	  for (int i=0; i<n1; i++)				
          {
	    matrix_t tmp3(ne,1);
	    for (int j=0;j<ne;j++)
	      tmp3[j][0] = gene1_eigvectors[j][i];
	    multMatrix(tmp2, tmp3, tmp);

	    for (int j=0;j<ne;j++) 
	      gene2_eigvectors[j][i] = (1/sqrt(gene1_eigvalues[i])) * tmp[j][0];
          }						       



 	  ////////////////////////////////////////////////////////
          // Step 5.  Calculate the gene1 (pxp) and gene2 (pxq)
          // coefficients used to create the canonical variates
          // associated

          //          with the p canonical correlations
          //          coeff.gene1 = t(e.vectors.gene1) %*% inv.sqrt.sig11 
          //          coeff.gene2 = t(e.vectors.gene2) %*% inv.sqrt.sig22 


          matrix_t gene1_coeff;
          matrix_t gene1_coeff;
	  
          multMatrix(t(gene1_eigvectors), inv_sqrt_I11, gene1_coeff);	// *********
          multMatrix(t(gene2_eigvectors), inv_sqrt_I22, gene2_coeff);   // *********

	  

 	  ////////////////////////////////////////////////////////
          // Step 6.  Compute the gene1 and gene2 canonical variates
          // associated with the highest canonical correlation (if <
          // 1) NOTE: the original variables of data need to be
          // standardised first! Otherwise, the resulting correlation
          // between variate.gene1 and variate.gene1 != estimated
          // cancor.
 
          // Choose the canonical correlation to use and store the
          // corresponding coefficients to create the can variates
          // Starts with the last can correlation and updates with the
          // previous if greater and < max(cancor)

          double max_cc = 0.9;
          int use_ccn = ne-1;						// *********
          double use_cc = sqrt(gene1_eigvalues[use_ccn]);
          vector<double> use_coeff;					// Make 1 single vector of length n12 (n1+n2), with coef for gene1, then gene2
          for (int i=0; i<n1; i++)
            use_coeff[i] = gene1_coeff[use_ccn][i];
          for (int i=0; i<n2; i++)
            use_coeff[n1+1+i] = gene1_coeff[use_ccn][i];

          for (int i=ne-2; i>0; i--)
          {
            if ( sqrt(gene1_eigvalues[i]) >= use_cc && sqrt(gene1_eigvalues[i]) < max_cc )
            {
              use_ccn = i;
              use_cc = sqrt(gene1_eigvalues[i]);
              for (int i=0; i<n1; i++)
                use_coeff[i] = gene1_coeff[use_ccn][i];
              for (int i=0; i<n2; i++)
                use_coeff[n1+1+i] = gene1_coeff[use_ccn][i];
            }
          }

          // Compute the gene1 and gene2 canonical variates after standardizing, as:
          // variates[x] = SUM ( (Xi - mean(Xi))/sqrt(var(Xi)) * coef(Xi[use_ccn,]) ), where i=1,...,p (gene1) or q (gene2)
          // pSNP_z = apply(pSNP,2,standardize)
          // variate.gene1 = pSNP_z[,1:n1] %*% gene1_coeff[1,]  
          // variate.gene2 = pSNP_z[,(n2+1):n12] %*% gene2_coeff[1,]  
          // cor(variate.gene1,variate.gene2)                 // Check		// *********
          // Use Mean and Variances calculated by calcGENEPIMeanVariance() to standardise
          // Can we do this?							// *********


	  vector<vector<double> > sigma(0);
	  vector<double> mean(0);

          calcGENEPIMeanVariance(pSNP, 
				 n1,n2,
				 false,
				 mean,
				 sigma,
				 P.sample , 
				 blperm[s1], 
				 blperm[s2] );
	  

	  // We now have mean, variance


          ///////////////////////////
          // Iterate over individuals
     
          for (int i=0; i<n; i++)
          {	      
	    
	    Individual * person = sample[i];

            ///////////////////////////
            // Iterate over all SNPs 

            double variate_gene1 = 0;
            double variate_gene2 = 0;
            int cur_snp = 0;

            for (int j=0; j<n12; j++)
            {
      
              cur_snp++;

              CSNP * ps = pSNP[i][j];
        
              bool a1 = ps->one;
              bool a2 = ps->two;
	  
              if ( a1 )
	      {		      
	        if ( a2 )    // 11 homozygote (genotye value = +1)
	        {
                  // SNP in gene 1?
                  if (cur_snp <= n1)
                  {
                    variate_gene1 = variate_gene1 + (1 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
                  else
                  {
                    variate_gene2 = variate_gene2 + (1 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
		} 
              }
	      else 
	      {
	        if ( ! a2  ) 		// 00 homozygote (genotye value = -1)
		{
                  // SNP in gene 1?
                  if (cur_snp <= n1)
                  {
                    variate_gene1 = variate_gene1 + (-1 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
                  else
                  {
                    variate_gene2 = variate_gene2 + (-1 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
	        }
		else
                {
                  // SNP in gene 1?
                  if (cur_snp <= n1)
                  {
                    variate_gene1 = variate_gene1 + (0 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
                  else
                  {
                    variate_gene2 = variate_gene2 + (0 - means[j])/sqrt(variances[j]) * use_coeff[j] ;
                  }
                }
              }
		  
	    } // Next individual	  

            
	    // Can i store the two new variables in clist[]?
	    
	    // TODO -- check for existing covariates, make space allocation
	    
            person->clist[0] = variate_gene1;
            person->clist[1] = variate_gene2;
      
          } // Next SNP 


          // If all working, at this stage we should have for each
          // individual in person->clist the values for the 2
          // variables (variate_gene1 and 2)

          // associated with max can correlation
  

          // These can now be used to fit a logistic regression using GLM
          // y = variate_gene1 + variate_gene2 +  variate_gene1 * variate_gene2



	  
	  //////////////////////////////////////////////////////////
	  // A new GLM
	  
	  Model * lm;
	  
	  
	  //////////////////////////////////////////////////////////
	  // Linear or logistic?
	  
	  if (par::bt)
	  {
	      LogisticModel * m = new LogisticModel(&P);
	      lm = m;
	  }
	  else
	  {
	      LinearModel * m = new LinearModel(&P);
	      lm = m;
	  }
	  

	  //////////////////////////////////////////////////////////
	  // Set missing data
	  
	  lm->setMissing();
	  
	  
	  
	  //////////////////////////////////////////////////////////
	  // Main effect of GENES
	  
	  lm->addCovariate(0); 
	  lm->label.push_back("GENE1");
	  
	  lm->addCovariate(1); 
	  lm->label.push_back("GENE2");

	  
	  // addInteraction() takes parameter numbers
	  // 0 intercept
	  // 1 X
	  // 2 Y
	  //   {interaction}
	  
	  
	  //////////////////////////
	  // Gene-based epistasis 
	  
	  lm->addInteraction(1,2);
	  lm->label.push_back("EPI");	  
	  
      
	  //////////////////////////////
	  // Build design matrix
	  
	  lm->buildDesignMatrix();
      
      
	  //////////////////////////////////////////////////
	  // Fit linear or logistic model (Newton-Raphson)

	  lm->fitLM();
      

	  ////////////////////////////////////////
	  // Check for multi-collinearity

	  lm->validParameters();
      

	  ////////////////////////////////////////
	  // Obtain estimates and statistic

//	  lm->displayResults(EPI,locus[l]);
// Make our own output
      
	  ////////////////////////////////////////
	  // Store statistic (1 df chisq)
	  // ( based on value of testParameter )

//	  results[l] = lm->getStatistic();
	  

	  /////////////////////////////////
	  // Clear up linear model
	  
	  delete lm;
     


	  
// 	  EPI << setw(8) << P.setname[s1] << " " 
// 	      << setw(4) << n1 << " "
// 	      << setw(8) << P.setname[s2] << " "
// 	      << setw(4) << n2 << " "
// 	      << setw(4) << i << " "
// 	      << setw(8) << sqrt(eig[i]) << " ";
// 	  if (!maxchisq) EPI << setw(8) << chisq << " ";
// 	  else EPI << setw(8) << "500+" << " ";
// 	  EPI << setw(8) << df << " "	    
// 	      << setw(8) << pvalue << " "
// 	      << "\n";
	  
//       }
// 	  EPI <<"\n";
	  
      }
  

  EPI.close();
  
  return res;
  
}  // END of CCA_logit







///////////////////////////////////////////////////////////
// SECOND CCA function: use for case-only analysis

vector<double> CCA_caseonly(bool perm, 
		            vector<vector<int> > & blperm,
		            vector<vector<int> > & blperm_control,
		            Set & S,
		            Plink & P)
  
{

  // For each pair of sets (genes), this function will do the following:
  // 
  // Step 1.  Construct a p x q SNP covariance matrix (sigma). Cases only. Partition into sigma11,sigma12,sigma21,sigma22.
  // Step 2.  Calculate the p x p matrix M = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11)).
  // Step 3.  Determine the p eigenvalues of M. The sqrt(eigen(M)) = p canonical correlations
  // Step 4.  Compute the significance of the largest canonical correlation using the c.d.f. approximation of Pillai (1964)


  ///////////////
  // Output file

  ofstream EPI;

  if (!perm)
  {  
   string f = par::output_file_name+".genepi";
   P.printLOG("\nWriting gene-based epistasis tests to [ " + f + " ]\n");
   EPI.open(f.c_str(), ios::out);
   EPI.precision(4);
  }
  

  //////////////////////////////////
  // Canonical correlation analysis

  vector<double> res;

  int ns = P.snpset.size();
  
  // Consider each pair of genes
  
  for (int s1=0; s1 < ns-1; s1++)
      for (int s2 = s1+1; s2 < ns; s2++)
      {
	  

 	  ////////////////////////////////////////////////////////
	  // Step 0. Do we want to recode the SNPs into additive coding?			// *********



 	  ////////////////////////////////////////////////////////
	  // Step 1. Construct covariance matrix (cases and controls together)
          //    And partition covariance matrix:
	  //    S_11  S_21
	  //    S_12  S_22
	  
	  int n1=0, n2=0;
	  
	  vector<vector<double> > sigma(0);
	  vector<double> mean(0);
	  vector<CSNP*> pSNP(0);
	  
	  /////////////////////////////
	  // List of SNPs for both loci
	  
	  for (int l=0; l<P.snpset[s1].size(); l++)
	      if ( S.cur[s1][l] )
	      {
		  pSNP.push_back( P.SNP[ P.snpset[s1][l] ] );
		  n1++;
	      }
	  for (int l=0; l<P.snpset[s2].size(); l++)
	      if ( S.cur[s2][l] )
	      {
		  pSNP.push_back( P.SNP[ P.snpset[s2][l] ] );
		  n2++;
	      }
	  

          // NOTE: we need to make sure that n1 < n2. Migth cause problems below if this is not the case.	// *********
	  int n12 = n1 + n2;
	  int ne = n1 < n2 ? n1 : n2;		// ne = min(p,q)

  
	  ///////////////////////////////////
	  // Construct covariance matrix (cases and controls together)
	  
	  P.setFlags(false);
	  vector<Individual*>::iterator person = P.sample.begin();
	  while ( person != P.sample.end() )
	  {
	      if ( (*person)->aff && !(*person)->missing )                       // MARF: this now chooses cases only
		  (*person)->flag = true;
	      person++;
	  }
	  
	  calcGENEPIMeanVariance(pSNP, 
				 n1,n2,
				 false,
				 mean, 
				 sigma, 
				 P.sample , 
				 blperm[s1], 
				 blperm[s2] );
	  
	  
	  ///////////////////////////
	  // Partition covariance matrix
	  
	  vector<vector<double> > I11;
	  vector<vector<double> > I12;
	  vector<vector<double> > I21;
	  vector<vector<double> > I22;
	  
	  sizeMatrix( I11, n1, n1);
	  sizeMatrix( I12, n1, n2);
	  sizeMatrix( I21, n2, n1);
	  sizeMatrix( I22, n2, n2);
	  
	  for (int i=0; i<n1; i++)
	      for (int j=0; j<n1; j++)
		  I11[i][j] = sigma[i][j];
	  
	  for (int i=0; i<n1; i++)
	      for (int j=0; j<n2; j++)
		  I12[i][j] = sigma[i][n1+j];
	  
	  for (int i=0; i<n2; i++)
	      for (int j=0; j<n1; j++)
		  I21[i][j] = sigma[n1+i][j];
	  
	  for (int i=0; i<n2; i++)
	      for (int j=0; j<n2; j++)
		  I22[i][j] = sigma[n1+i][n1+j];
	  
	  
 	  ////////////////////////////////////////////////////////
	  // Step 2. Calculate the p x p matrix M1 = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))
       	  
	  I11 = msqrt(I11);
	  I11 = svd_inverse(I11);
	  I22 = svd_inverse(I22);
	  
	  matrix_t tmp;
	  matrix_t M1;
	  
	  multMatrix(I11, I12, tmp);
	  multMatrix(tmp, I22, M1);
	  multMatrix(M1, I21, tmp);
	  multMatrix(tmp, I11, M1);



 	  ////////////////////////////////////////////////////////
	  // Step 3. Compute the p eigenvalues of M1. The sqrt(eigen(M)) = p canonical correlations
	  
	  vector_t gene1_eigvalues = eigenvalues(M1);
	  
	  // Display eigen values
	  cout << "E (1) values\n";
	  display(gene1_eigvalues);
	  cout << "\n";

	  sort(gene1_eigvalues.begin(),gene1_eigvalues.end(),greater<double>());
	  
	  for (int i=0; i<ne; i++)
	  {
	      if ( gene1_eigvalues[i] < 0 ) 
		  gene1_eigvalues[i] = 0;
	      else if ( gene1_eigvalues[i] >= 1 )
		  gene1_eigvalues[i] = 0.999;
	  }

	  // Display canonical correlations, sorted
	  cout << "All canonical correlations \n";
	  display(sqrt(gene1_eigvalues));			// *********
	  cout << "\n";
	  



 	  ////////////////////////////////////////////////////////
          // Step 4.  Compute the significance of the largest canonical correlation using the c.d.f. approximation of Pillai (1964)

          int p=n1;
          int s=n1;
          int q=n2;
          int N=ncase;
          if (n1 > n2)
          {
            p=n2; q=n1; s=n2;
          }
          double m = 0.5 * (q-p-1);
          double n = 0.5 * (N-p-q-2);
          double x = gene1_eigvalues[1] * gene1_eigvalues[1];
          double Pvalue = pLargRoot(s,m,n,x);


      }
  

  EPI.close();
  
  return res;

}    // END of CCA_logit






void calcGENEPIMeanVariance(vector<CSNP*> & pSNP,
			    int n1,
			    int n2,
			    bool perm, 
			    vector<double> & mean,
			    vector<vector<double> > & variance,
			    vector<Individual*> & sample,
			    vector<int> & gp1,
			    vector<int> & gp2 )
  
{

  // Calculate mean and variance for n1+n2 x n1+n2 matrix

  // Individual order in n1 , n2 deteremined by g1, g2
  // (i.e. block-based permutation)
  
  // Under permutations, mean and variances won't change
  // Store means only for now

  int n = gp1.size();
  
  int nss = pSNP.size();
  
  // Original calculation?
  if (!perm)
    mean.resize(nss,0);
  
  vector<int> cnt(nss,0);
  variance.resize(nss);
  
  for (int j=0; j<nss; j++)
    variance[j].resize(nss,0);
  
      
  /////////
  // Means

  /////////////////////////////////
  // Consider each SNP in this set
  
  for (int j=0; j<nss; j++)
    {
      
      CSNP * ps = pSNP[j];
      
      ///////////////////////////
      // Iterate over individuals
      
      for (int i=0; i<n; i++)
	{	      
	  // Only need to look at one perm set
	  bool a1 = ps->one[gp1[i]];
	  bool a2 = ps->two[gp2[i]];
	  
	  if ( a1 )
	    {		      
	      if ( a2 ) // 11 homozygote
		{
		  mean[j]++;
		  cnt[j]++;
		}
	    }
	  else 
	    {
	      cnt[j]++;
	      if ( ! a2  ) // 00 homozygote
		mean[j]--;
	      
	    }
	  
	} // Next individual	  
      
    } // Next SNP in set
  
  for (int j=0; j<nss; j++)
    mean[j] /= (double)cnt[j];
  
  
  
  /////////////////////////////////////
  // Iterate over pairs of SNPs in SETs
  
  // First SNP 
  for (int j1=0; j1<nss; j1++)
    {
      CSNP * ps1 = pSNP[j1];

      // Second SNP
      for (int j2=0; j2<nss; j2++)
	{
	  CSNP * ps2 = pSNP[j2];

	  // Iterate over individuals
	  
	  for (int i=0; i<n; i++)
	    {
	      
	      bool a1, a2;
	      if (j1<n1)
		{
		  a1 = ps1->one[gp1[i]];
		  a2 = ps1->two[gp1[i]];
		}
	      else
		{
		  a1 = ps1->one[gp2[i]];
		  a2 = ps1->two[gp2[i]];
		}
	      
	      bool b1, b2;
	      if (j1<n1)
		{
		  b1 = ps2->one[gp1[i]];
		  b2 = ps2->two[gp1[i]];
		}
	      else
		{
		  b1 = ps2->one[gp2[i]];
		  b2 = ps2->two[gp2[i]];
		}
	      
	      
	      // Mean substitution
	      double v1=mean[j1], v2=mean[j2];
	      
	      // First SNP
	      if ( a1 )
		{		      
		  if ( a2 )    // 11 homozygote
		    {
		      v1 = 1;
		    }
		}
	      else 
		{
		  if ( ! a2  ) // 00 homozygote
		    {
		      v1 = -1;
		    }
		  else
		    v1 = 0;       // 01 heterozygote
		}
	      
	      
	      // Second SNP
	      if ( b1 )
		{		      
		  if ( b2 )    // 11 homozygote
		    {
		      v2 = 1;
		    }
		}
	      else 
		{
		  if ( ! b2  ) // 00 homozygote
		    {
		      v2 = -1;
		    }
		  else
		    v2 = 0;       // 01 heterozygote
		}
	      
	      
	      // Contribution to covariance term
	      variance[j1][j2] += ( v1 - mean[j1] ) * ( v2 - mean[j2] );
	      
		} // Next individual
	    } // Second SNP
	} // First SNP
        
      // Make symmetric covariance matrix
      for (int i=0; i<nss; i++)
	for (int j=i; j<nss; j++)
	  {
	    variance[i][j] /= (double)(n);
	    variance[j][i] = variance[i][j];
	  }
      
      return;
      
}
  
