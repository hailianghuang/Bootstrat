

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
#include <sstream>
#include <cmath>
#include <vector>
#include <map>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"


void HaploPhase::haplotypicMH(map<int,int> & tests, int nt)
{

  vector<double> caseN(nt);
  vector<double> controlN(nt);
  
  // Consider each individual
  for (int i=0; i<P.n; i++)
    {

      Individual * person = P.sample[i];

      // Case?
      if ( ! person->missing )
	{
	  
	  if (person->aff )
	    {
	      
	      for (int z = 0 ; z < hap1[i].size(); z++)
		{
				  
		  map<int,int>::iterator i1 = tests.find(hap1[i][z]);
		  map<int,int>::iterator i2 = tests.find(hap2[i][z]);

		  if ( i1 != tests.end() )
		    {
		      if (!ambig[i]) 
			caseN[i1->second]++;
		      else
			caseN[i1->second] += pp[i][z];
		    }
		  
		  if ( i2 != tests.end() )
		    {
		      if (!ambig[i]) 
			caseN[i2->second]++;
		      else
			caseN[i2->second] += pp[i][z];
		    }


		}  
	    }
	  // Or control?
	  else 
	    {
	      for (int z = 0 ; z < hap1[i].size(); z++)
		{
		  
		  map<int,int>::iterator i1 = tests.find(hap1[i][z]);
		  map<int,int>::iterator i2 = tests.find(hap2[i][z]);
		  
		  if ( i1 != tests.end() )
		    {
		      if (!ambig[i]) 
			controlN[i1->second]++;
		      else
			controlN[i1->second] += pp[i][z];
		    }
		  
		  if ( i2 != tests.end() )
		    {
		      if (!ambig[i]) 
			controlN[i2->second]++;
		      else
			controlN[i2->second] += pp[i][z];
		    }

		  
		  
		}  	  
	    }
	}
      
    } // next individual
  

  HTEST << setw(10) << hname << " ";
  
  if (nt==2)
    {

      // Find test haplotype
      int hh=0;
      map<int,int>::iterator i1 = tests.begin();
      while ( i1 != tests.end() )
	{
	  if ( i1->second == 0 )
	    hh = i1->first;
	  i1++;
	}

      HTEST << setw(12) << haplotypeName(hh) << " ";
      
      HTEST << setw(10) << caseN[0] / ( caseN[0] + caseN[1] ) << " "
	    << setw(10) << controlN[0] / ( controlN[0] + controlN[1] ) << " ";
      
    }
  else
    {
      HTEST << setw(12) << "OMNIBUS" << " "
	    << setw(10) << "NA" << " "
	    << setw(10) << "NA" << " ";
    }

  vector<double> rowT(nt);
  double caseT = 0;
  double controlT = 0;

  for (int h=0; h<nt; h++)
    {
      rowT[h] = caseN[h] + controlN[h];
      caseT += caseN[h];
      controlT += controlN[h];       
    }

  double chi2 = 0;
  for (int h=0; h<nt; h++)
    {
      double exp = ( rowT[h] * caseT ) / (caseT + controlT);
      chi2 += ( ( caseN[h] - exp ) * ( caseN[h] - exp ) ) / exp ; 
      
      exp = ( rowT[h] * controlT ) / (caseT + controlT);
      chi2 += ( ( controlN[h] - exp ) * ( controlN[h] - exp ) ) / exp ;  
    }

  
  if ( realnum(chi2) )
    {
      HTEST << setw(10) << chi2 << " "
	    << setw(4) << nt-1 << " "
	    << setw(10) << chiprobP(chi2,nt-1) << " ";
    }
  else
    {
      HTEST << setw(10) << "NA" << " "
	    << setw(4) << "NA" << " "
	    << setw(10) << "NA" << " ";
    }
  
  for (int snps=0; snps<ns-1; snps++)
    HTEST << P.locus[S[snps]]->name << "|";
    
  HTEST << P.locus[S[ns-1]]->name << "\n";


}



///////////////////////////////////
// Multimarker test with weighting

void HaploPhase::haplotypicWeightedCC()
{

  vector_t weights;

  for (int i=0; i<nh; i++)
    {
      map<string,double>::iterator whap = new_pred_weighted_allele[current].find( haplotypeName(i) );

      if ( whap != new_pred_weighted_allele[current].end() )
	{
	  weights.push_back( whap->second );
	}
      else
	{
	  weights.push_back( 0 );
	}
    }
  



   vector<double> caseN(2);
   vector<double> controlN(2);
  
   // Consider each individual
   for (int i=0; i<P.n; i++)
     {

       Individual * person = P.sample[i];

       // Case?
       if ( ! person->missing )
	 {
	   
	   if (person->aff )
	     {
	       
	       for (int z = 0 ; z < hap1[i].size(); z++)
		 {
		   
		   int h1 = hap1[i][z];
		   int h2 = hap2[i][z];
		   
		   if (!ambig[i]) 
		     {
		       caseN[0] += weights[h1];
		       caseN[1] += 1-weights[h1];
		     }
		   else
		     {
		       caseN[0] += weights[h1] * pp[i][z];
		       caseN[1] += (1-weights[h1]) * pp[i][z];
		     }

		   if (!ambig[i]) 
		     {
		       caseN[0] += weights[h2];
		       caseN[1] += 1-weights[h2];
		     }
		   else
		     {
		       caseN[0] += weights[h2] * pp[i][z];
		       caseN[1] += (1-weights[h2]) * pp[i][z];
		     }

		 }
	     }
	   // Or control?
	   else 
	     {
	       for (int z = 0 ; z < hap1[i].size(); z++)
		 {
		   
		   int h1 = hap1[i][z];
		   int h2 = hap2[i][z];
		   
		   if (!ambig[i]) 
		     {
		       controlN[0] += weights[h1];
		       controlN[1] += 1-weights[h1];
		     }
		   else
		     {
		       controlN[0] += weights[h1] * pp[i][z];
		       controlN[1] += (1-weights[h1]) * pp[i][z];
		     }

		   if (!ambig[i]) 
		     {
		       controlN[0] += weights[h2];
		       controlN[1] += 1-weights[h2];
		     }
		   else
		     {
		       controlN[0] += weights[h2] * pp[i][z];
		       controlN[1] += (1-weights[h2]) * pp[i][z];
		     }


		 }
	     }
	 }

     } // next individual

   

   HTEST << setw(10) << hname << " ";
   
   // set hh to integer
   
   HTEST << setw(12) << new_map[current]->allele1 << " ";
   
   HTEST << setw(10) << caseN[0] / ( caseN[0] + caseN[1] ) << " "
	 << setw(10) << controlN[0] / ( controlN[0] + controlN[1] ) << " ";
   
  
   vector<double> rowT(2);
   double caseT = 0;
   double controlT = 0;
  
   for (int h=0; h<2; h++)
     {
       rowT[h] = caseN[h] + controlN[h];
       caseT += caseN[h];
       controlT += controlN[h];       
     }

   double chi2 = 0;
   for (int h=0; h<2; h++)
     {
       double exp = ( rowT[h] * caseT ) / (caseT + controlT);
       chi2 += ( ( caseN[h] - exp ) * ( caseN[h] - exp ) ) / exp ; 
     
       exp = ( rowT[h] * controlT ) / (caseT + controlT);
       chi2 += ( ( controlN[h] - exp ) * ( controlN[h] - exp ) ) / exp ;  
     }

  
   if ( realnum(chi2) )
     {
       HTEST << setw(10) << chi2 << " "
 	    << setw(4) << 1 << " "
 	    << setw(10) << chiprobP(chi2,1) << " ";
     }
 else
     {
       HTEST << setw(10) << "NA" << " "
 	    << setw(4) << "NA" << " "
 	    << setw(10) << "NA" << " ";
     }
  
   for (int snps=0; snps<ns-1; snps++)
     HTEST << P.locus[S[snps]]->name << "|";
    
   HTEST << P.locus[S[ns-1]]->name << "\n";

}

