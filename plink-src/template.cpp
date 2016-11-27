
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
#include <fstream>
#include <iomanip>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"

vector_t Plink::myFunction1(bool print_results, Perm & perm)
{

  /////////////////////////////////////////////////////////////////////  
  // Many parameters passed in via static, global "par" class
  // These controlled by Shaun (i.e. request for me to add)
  // Defined in options.*
  // Set from command line in parse.cpp


  ////////////////////////////
  // Several inbuilt functions
  // Probably requires consultation with me to use these
  
  // e.g.
  //    Haplotype phasing (see phase.cpp)
  //    Matrix operations, e.g. inversion (see helper.cpp)
  //    SVD / linear and logistic model fitting (see glm.cpp)
  //    Statistical functions, distributions (see helper.cpp)
  //    Permutation (see perm.cpp)


  /////////////////////////////////////////////////////////////////////
  // SNP-major versus individual-major coding
  
  if ( par::SNP_major )
    SNP2Ind();   // helper function (helper.h)
  

  /////////////////////////////////////////////////////////////////////
  // Use "printLOG(string)" to write messages
  
  printLOG("This is my function # " + int2str(1) 
	   + " " + dbl2str(1.0) + "\n");


  /////////////////////////////////////////////////////////////////////
  // Write results to file
  
  string f = par::output_file_name + ".myfunc1";
  ofstream OUT(f.c_str(),ios::out);
  printLOG("Writing results to [ " + f + " ]\n");
  

  /////////////////////////////////////////////////////////////////////
  // "sample" is a vector of pointers to "Individual"'s
  
  iIndividual person = sample.begin();
      
  while ( person != sample.end() )
    {
      
      if ( (*person)->founder ) 
	OUT << "Founder\t";
      
      if ( (*person)->aff )
	OUT << "Affected\t";
      
      if ( (*person)->sex )
	OUT << "Male\t";
      
      if ( (*person)->missing )
	OUT << "Missing\t";
      
      OUT << "Full phenotype = " << (*person)->phenotype << "\n";
      
      person++;
    }
  
  OUT << "\n\n";



  /////////////////////////////////////////////////////////////////////
  // Individual major iteration over genotypes, using array indexes

  for (int i=0; i<n; i++)
    {
      
      Individual * person = sample[i];
      
      OUT << person->fid << " " 
	  << person->iid << " ";

      if ( person->aff )
	OUT << "AFF   ";
      else
	OUT << "UNAFF ";
      	          
      	  
      //////////////////////////////////////////
      // Consider each SNP for this person
      
      for (int l=0; l<nl_all; l++)
	{
	  
	  //////////////////////////////////////////
	  // Skip possibly

// 	  if (par::adaptive_perm && !perm.test[l])
// 	    continue;
	  
	  
	  //////////////////////////////////////////
	  // X-chromosome, haploid?
	  
	  if ( par::chr_sex[locus[l]->chr] ||
	       par::chr_haploid[locus[l]->chr] )
	    continue;

	  
 	  OUT << locus[l]->chr << " "
 	      << locus[l]->name << " : ";
	  
	  // Also bp member of Locus with physical position

	  //////////////////////////////////////////
	  // Get and parse genotypes

	  bool one = sample[i]->one[l];
	  bool two = sample[i]->two[l];
	  

	  if ( one )
	    {
	      if ( two )
		OUT << locus[l]->allele2 << "/" 
		    << locus[l]->allele2 << "  ";
	      else
		OUT << "0/0  ";
	    }
	  else
	    {
	      if ( two )
		OUT << locus[l]->allele1 << "/" 
		    << locus[l]->allele2 << "  ";
	      else
		OUT << locus[l]->allele1 << "/" 
		    << locus[l]->allele1 << "  ";
	    }


	} // Next SNP
      

      OUT << "\n";

    } // Next Individual





  /////////////////////////////////////////////////////////////////////
  // SNP-major mode genotype, using iterators

  // Would require:
  // if ( !par::SNP_major )
  //   Ind2SNP();   
  // i.e. these will just be empty in individual-major mode

  // Iterator over SNPs

  iSNP s = SNP.begin();
  
  while ( s != SNP.end() )
    {
      
      // Iterate over individuals

      iAllele i1 = (*s)->one.begin();
      iAllele i2 = (*s)->two.begin();
      iIndividual gperson = sample.begin();                             

      while ( gperson != sample.end() )
	{
	  
	  // "Phenotype" individual (typically self, but allows for permutation)

	  Individual * pperson = (*gperson)->pperson;


	  // Alleles for this SNP

	  bool one = *i1;
	  bool two = *i2;
	  
	  // { ... do something here....}

	  
	  // Next individual

	  gperson++;
	  i1++;
	  i2++;

	}
      
      // Summarize counts for SNP, calculate statistic, report output, etc

      double chisq = 22;
      int df = 1;
      OUT << "pvalue = " << chiprobP(chisq,df) << "\n";
      
    }




  


  /////////////////////////////////////////////////////////////////////
  // Consider families

  // Requires that family structure have been previously 
  // set (ask me)

  OUT << "\n\nFamily Structure\n";
  
  for (int f = 0 ; f < family.size() ; f++)
    {
      
      Family * fam = family[f];
     
      // pat, mat and kid[] are pointers to Individuals
      
      if ( fam->singleton )
	{
	  OUT << "SINGLETON(S)\t" 
	       << fam->kid[0]->fid << " : ";
	  for (int k=0; k < fam->kid.size() ;k++)
	    OUT << fam->kid[k]->iid << " ";
	  OUT << "\n";
	}
      else if ( fam->sibship )
	{
	  OUT << "SIBSHIP  \t" << fam->kid[0]->fid << " : ";
	  for ( int k=0; k<fam->kid.size(); k++)
	    OUT << fam->kid[k]->iid << " ";
	  OUT << "\n";
	}
      else if ( fam->parents )
	{
	  OUT << "W/ PARENTS\t" << fam->pat->fid << " : ";
	  OUT << fam->pat->iid << " x " << fam->mat->iid << " -> ";
	  for ( int k=0; k<fam->kid.size(); k++)
	    OUT << fam->kid[k]->iid << " ";
	  OUT << "\n";	     
	}
      else
	OUT << "UNDEFINED\t" 
	     << fam->pat->fid << " "
	     << fam->pat->iid << "\n";

    }
  


  /////////////////////////////////////////////////////////////////////
  // Close output stream

  OUT.close();

  
  /////////////////////////////////////////////////////////////////////
  // Typically, return a vector of statistics

  vector_t results(nl_all);
  
 
  return results;

}

