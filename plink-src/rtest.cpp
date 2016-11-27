

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2010 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <iomanip>

#include "rtest.h"
#include "plink.h"
#include "options.h"
#include "helper.h"
#include "perm.h"
#include "sets.h"
#include "model.h"
#include "stats.h"


// TODO: when using pooled data, freqs are not populated yet

void Plink::rTestWrapper()
{


  if ( par::rtest_type == RTEST_VARIABLE_THRESHOLD )
    printLOG("Performing (weighted) VT rare-variant test\n");
  else if ( par::rtest_type == RTEST_FREQWEIGHT )
    printLOG("Performing frequency-weighted rare-variant test\n");
  else if ( par::rtest_type == RTEST_VANILLA )
    printLOG("Performing rare-variant burden test, " 
	     + dbl2str( par::rtest_freq_threshold ) + " threshold\n");


  if ( ! par::permute ) 
    error("You need to specify --mperm {N} or --perm when using --vt-test");

  bool useSets = par::read_set || par::make_set;


  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Do we need to use full GLM? If covariates, 2-sided tests            //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  if ( par::rtest_hyp == RTEST_2SIDED || par::clist ) 
    par::rtest_glm = true;
  
  // Add an extra slot for the rare-variant count in the set
  
  if ( par::rtest_glm )
    {
      par::clist = true;
      ++par::clist_number;            
      clistname.push_back("RVCNT");
      for (int i=0; i<n; i++)
	sample[i]->clist.push_back(0);      
    }
  
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Temporary measure: flip phenotype for basic test                    //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  if ( ! par::rtest_glm ) 
    {
      for (int i=0; i<n; i++)
	sample[i]->phenotype *= -1;
    }

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Using weights?  Assumes a 0..1 scale, where 1 is high weight        //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  vector_t wgt(nl_all,0); 

  ifstream IN1;

  map<string,int> mlocus;
  for (int l=0; l<nl_all; l++) 
    mlocus.insert(make_pair( locus[l]->name, l ));
    
  if ( par::rtest_weights )
    {
      checkFileExists( par::rtest_weight_file );
        printLOG("Reading weights from [ " + par::rtest_weight_file + " ]\n");
	
	IN1.open( par::rtest_weight_file.c_str() , ios::in );
	int wc= 0;
	while ( ! IN1.eof() )
	  {
	    string snp;
	    string w_str;
	    double w;
	    IN1 >> snp >> w_str;
	    map<string,int>::iterator i = mlocus.find( snp );
	    if ( i == mlocus.end() )
	      continue;
	    if ( from_string<double>( w , w_str , std::dec ) )
	      {  
		// Assume 0..1 coding for weight
		wgt[ i->second ] = w < 0 ? 0 : w > 1 ? 1 : w;
	      }
	    ++wc;
	  }
	printLOG("Read weights for " + int2str( wc ) + " variants\n");  
	IN1.close();
    } 
  
  

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Construct list of rare individual/genotype entries to cycle through //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  vector<int> counts(nl_all);

  if ( ! par::pool_input )
    {

      // Using all variants, or a set?

      if ( useSets ) 
	{
	  
	  rdata.resize( pS->snpset.size() );
	  
	  for (int s=0; s < pS->snpset.size(); s++)
	    {
	      for (int j=0;j<pS->snpset[s].size();j++)
		{	      
		  for (int i=0; i<n; i++)
		    {
		      int l = snpset[s][j];
		      
		      bool s1 = SNP[ l ]->one[i];
		      bool s2 = SNP[ l ]->two[i];
		      
		      // Store only minor alleles
		      
		      if ( !s1 )
			{
			  if ( !s2 )
			    {
			      rdata[s].push_back( RVData(i,l,2) );
			      counts[l] += 2;
			    }
			  else 
			    {
			      rdata[s].push_back( RVData(i,l,1) );
			      ++counts[l];
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  // A single set
	  rdata.resize(1);
	  
	  for (int l=0; l<nl_all; l++)
	    for (int i=0; i<n; i++)
	      {
		bool s1 = SNP[l]->one[i];
		bool s2 = SNP[l]->two[i];
		
		// Store only minor alleles
		
		if ( !s1 )
		  {
		    if ( !s2 )
		      {
			rdata[0].push_back( RVData(i,l,2) );
			counts[l] += 2;
		      }
		    else 
		      {
			rdata[0].push_back( RVData(i,l,1) );
			++counts[l];	
		      }
		    
		  }
	      }
	}
    }



  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Populate locus frequencies                                          //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  
  if ( par::pool_input )
    {
      
      vector<vector<RVData> >::iterator si = rdata.begin();
      vector_t tc(nl_all);
      
      set<int2> counted;
      
      for (int l=0;l<nl_all;l++)
	locus[l]->freq = 0;
      
      int scnt=0;
      while ( si != rdata.end() )
	{
	  
	  vector<RVData>::iterator i = si->begin();      
	  while ( i != si->end() )
	    {	  
	      // Only count each variant once
	      if ( counted.find( int2( i->l , i->n ) ) == counted.end() )
		{
		  locus[ i->l ]->freq += i->c;
		  tc[ i->l ] += i->tc;
		  counted.insert(int2(i->l,i->n));
		}
	      ++i;
	    }	  
	  ++si;
	  ++scnt;
	}
      

//       for (int l=0;l<nl_all;l++)
// 	{
// 	  cout << "MAF for locus " 
// 	       << l << " is " 
// 	       << locus[l]->freq << " " 
// 	       << tc[l] << "\n";

// 	  locus[l]->freq /= tc[l];

// 	  cout << "MAF for locus " << l 
// 	       << " is " << locus[l]->freq 
// 	       << " " << tc[l] << "\n\n";
// 	}
    }


  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Set up permutations                                                 //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  int nTests = rdata.size();
  
  pperm->setTests( nTests );
  pperm->setPermClusters(*this);
  pperm->originalOrder();
  


  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Frequency-based weight for each SNP                                 //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  vector_t fwgt;

  if ( par::rtest_type == RTEST_FREQWEIGHT )
    {

      // CHECK: why this form of freq (1+c)/(2+2n)? 
      //       double f2 = (1.0+tot)/(2.0+2.0*n);
      //       double fweight = 1 / sqrt( f2 * (1-f2) );

      fwgt.resize(nl_all,1);
      for (int l=0; l<nl_all; l++)
	{

	  // NOTE: currently ignore missing phenotype/genotype status here

	  double f = ( 1.0 + counts[l] ) / ( 2.0 + 2.0 * n );
	  
	  fwgt[l] =  f > 0 && f < 1 ? 
	    1.0 / sqrt( f * (1-f) ) : 
	    0 ;
	}
    }

  

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Mean phenotype value                                                //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  double pmean = 0; int n1 = 0;
  for (int i=0; i<n; i++) 
    {
      if ( ! sample[i]->missing )
	{
	  pmean += sample[i]->phenotype;
	  ++n1;
	}
    }
  pmean /= (double)n1;


    
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Summary statistics                                                  //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  ofstream SUMM;
  ofstream SUMM2;
    
  if ( par::rtest_type == RTEST_VARIABLE_THRESHOLD )
    {
      
      printLOG("Writing VT-test summary to [ " 
	       + par::output_file_name + ".vt ]\n");
      
      printLOG("Writing VT-test variant summary to [ " 
	       + par::output_file_name + ".vt.var ]\n");
      
      SUMM.open( (par::output_file_name+".vt").c_str() , ios::out );
      SUMM2.open( (par::output_file_name+".vt.var").c_str() , ios::out );
    
      
      SUMM << setw(20) << "SET" << " "
	   << setw(8) << "NSNP" << " "
	   << setw(6) << "TC" << " "
	   << setw(10) << "TF" << " "
	   << setw(8) << "NSNP2" << " "; 
      if  ( par::bt ) 
	SUMM << setw(8) << "CNTA" << " "
	     << setw(8) << "CNTU" << " ";
      else
	{
	  SUMM << setw(8) << "CNT1" << " "
	       << setw(8) << "CNT0" << " "
	       << setw(8) << "MEAN1" << " "
	       << setw(8) << "MEAN0" << " ";
	}
      SUMM << "\n";
      
      SUMM2 << setw(20) << "SET" << " "
	    << setw(par::pp_maxsnp) << "SNP" << " ";
      
      if ( par::rtest_weights)
	SUMM2 << setw(8) << "WGT" << " ";
      
      if ( par::rtest_type == RTEST_FREQWEIGHT )
	SUMM2 << setw(8) << "FWGT" << " ";
      
      SUMM2 << setw(8) << "CNT" << " "
	    << setw(12) << "F" << " ";
      
      if ( attribs )
	SUMM2 << "ATTRIB";
      SUMM2 << "\n";
      
      SUMM.precision(4);
      SUMM2.precision(4);
    }

  
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Analyse original data                                               //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////


  vector_t original = doRTests(true, pmean, wgt, fwgt, rdata, counts, SUMM, SUMM2);


  SUMM.close();
  SUMM2.close();
  


  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Permutations                                                        //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////


  bool finished = false;

  while ( ! finished )
    {
      pperm->permuteInCluster();
      vector_t pr = doRTests(false, pmean, wgt, fwgt, rdata, counts, SUMM, SUMM2);
      finished = pperm->update(pr, original);
    }
  


  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Output results                                                      //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  string label = par::rtest_type == RTEST_VARIABLE_THRESHOLD ? ".vt.mperm" : 
    par::rtest_type == RTEST_FREQWEIGHT ? ".fw.mperm" : ".rv.mperm";


  printLOG("Writing permutation results to [ " 
	   + par::output_file_name 
	   + label + " ]\n");

  ofstream POUT( ( par::output_file_name + label ).c_str() , ios::out ); 
  
  int ns = useSets ? setname.size() : 1;

  // For now, skip EMP2

  POUT << setw(12) << "SET" << " "       
       << setw(12) << "EMP1" << "\n";
    //       << setw(12) << "EMP2" << "\n";

  POUT.precision(6);

  for (int s = 0 ; s < ns ; s++) 
    {
      if ( useSets ) 
	POUT << setw(12) << setname[s] << " ";
      else
	POUT << setw(12) << "ALL" << " ";	
      POUT << setw(12) << pperm->pvalue(s) << "\n";

      //POUT << setw(12) << pperm->max_pvalue(s) << "\n";
    }
  
  POUT.close();

}






//////////////////////////////////////////////////////
//                                                  //
// Main function to perform RV-tests                //
//                                                  //
//////////////////////////////////////////////////////

vector_t Plink::doRTests(bool print, 
			 double pmean, 
			 vector_t & wgt, 
			 vector_t & fwgt, 
			 vector<vector<RVData> > & rdata,
			 vector<int> & allele_counts , 
			 ofstream & SUMM,
			 ofstream & SUMM2 )
{  
  
  bool useSets = par::read_set || par::make_set;

  // One statistic per set tested

  int ns = rdata.size();

  vector_t results;

  // Current thresholds on tot # of alleles per (rare) variant

  const int MINCOUNT = 1;
  const int MAXCOUNT = 1000; 



  //////////////////////////////////////////////////////
  //                                                  //
  // Iterate over sets                                //
  //                                                  //
  //////////////////////////////////////////////////////

  int scnt = 0;
  vector<vector<RVData> >::iterator si = rdata.begin();
  while ( si != rdata.end() )
    {
      
      double score = 0;
      int bestk = 0;
 

      //////////////////////////////////////////////////////
      //                                                  //
      // Use fast, non-GLM tests                          //
      //                                                  //
      //////////////////////////////////////////////////////
      
      if ( ! par::rtest_glm )
	{
	  
	  // Weighted variable threshold test
	  // key = sample count 
	  // value = component of statistic for vt-test

	  map<int,double> sumx;
	  map<int,double> sum1;
	  map<int,double> sum2;
	  
	  
	  //////////////////////////////////////////////////////
	  //                                                  //
	  // Iterate over variants in set                     //
	  //                                                  //
	  //////////////////////////////////////////////////////
	  
	  vector<RVData>::iterator i = si->begin();
	  
	  while ( i != si->end() )
	    {
	      // Phenotype
	      double p = sample[ i->n ]->pperson->phenotype;

	      // Allele count for this individual (1/2)
	      int c = i->c;
	      
	      // Total sample allele count/freq
	      int tot = allele_counts[ i->l ];
	      double f = locus[ i->l ]->freq;

	      // Weight for this variant 
	      double pweight = wgt[ i->l ];
	      

	      // CHECK: Adjust polyphen score for common alleles

	      if ( par::rtest_ppweights )
		if ( pweight < 1.0 && f >= 0.01 ) pweight = 0.5;
	      
	      
	      // 
	      // Calculate test statistic as appropriate
	      //

	      if ( par::rtest_type == RTEST_VANILLA )
		{
		  if ( f < par::rtest_freq_threshold )
		    score += par::rtest_weights ? p * c * pweight : p * c ;
		}
	      
	      
	      //
	      // Frequency-weighted test
	      //

	      else if ( par::rtest_type == RTEST_FREQWEIGHT )
		{		  
		  score += par::rtest_weights ? 
		    p * c * fwgt[ i->l ] * pweight : 
		    p * c * fwgt[ i->l ] ;
		}
	  

	      //
	      // Variable threshold test 
	      //

	      else if ( par::rtest_type == RTEST_VARIABLE_THRESHOLD )
		{
		  		  		  
		  if ( tot < MAXCOUNT ) 
		    {
		      if ( par::rtest_weights )
			{
			  sumx[ tot ] += p * c * pweight;
			  sum1[ tot ] += c * pweight;
			  sum2[ tot ] += c * c * pweight * pweight; 
			}
		      else
			{
			  sumx[ tot ] += p * c;
			  sum1[ tot ] += c;
			  sum2[ tot ] += c * c;
			}
		    }
		}
	      

	      //
	      // Next rare variant
	      //

	      ++i;
	    }
	
      
	  //////////////////////////////////////////////////////
	  //                                                  //
	  // For VT-test, iterate over thresholds             //
	  //                                                  //
	  //////////////////////////////////////////////////////

	  if ( par::rtest_type == RTEST_VARIABLE_THRESHOLD )
	    {
	      
	      double tsumx = 0, tsum1 = 0, tsum2 = 0;

	      map<int,double>::iterator i = sumx.begin();

	      // Consider each possible threshold
	      
	      while ( i != sumx.end() )
		{

		  int th = i->first;
		  
		  tsumx += i->second;
		  tsum1 += sum1[th];
		  tsum2 += sum2[th];
		  
		  double tscore = ( tsumx - tsum1 * pmean ) / sqrt( tsum2 );
		  
		  if ( tscore > score ) 
		    {
		      score = tscore;
		      bestk = th; // and track threshold
		    }


		  //
		  // Next possible count-threshold
		  //

		  ++i;
		}

	    }
	}
      

    
      
      //////////////////////////////////////////////////////
      //                                                  //
      // Use alternative GLM tests, for covariates, etc   //
      //                                                  //
      //////////////////////////////////////////////////////
      
      if ( par::rtest_glm )
	{	  
	  
	  int start = par::rtest_type == 1 ? 1 : MAXCOUNT;
	  int finish = MAXCOUNT;
	  
	  score = 0;

	  for (int k = start; k <= finish; k++)
	    {
	      
	      vector_t sc(n,0);
	      
	      vector<RVData>::iterator i = si->begin();
	      while ( i != si->end() )
		{	      
		  Individual * person = sample[ i->n ]->pperson;
		  int c = i->c;
		  
		  double f = locus[ i->l ]->freq;	      
		  int cnt = allele_counts[ i->l ];
		  double weight = wgt[ i->l ];
		  
		  if(cnt >= par::rtest_polycount_threshold )
		    weight = 1.0;
		  
		  if ( cnt <= k )
		    sc[i->n] += i->c * weight;
		  
		  ++i;
		}

	      
	      //	      
	      // Place (weighted) count as last covariate
	      //

	      for (int i=0; i<n; i++)
		sample[i]->clist[ par::clist_number - 1 ] = sc[i];

	      
	      //
	      // Regression of phenotype on count + covariates
	      //

	      glmAssoc( false , *pperm );


	      //
	      // Get result (retain X^2 statistic)
	      //

	      model->testParameter = par::clist_number;
	      vector_t b = model->getCoefs();
	      double chisq = model->getStatistic();
	      double pvalue = chiprobP(chisq,1);
	      double beta = b[ par::clist_number ];
	      double thisScore = model->isValid() ? 
		model->getStatistic() :
		0 ;

	      
	      //
	      // Is this a good score?
	      //

	      if ( thisScore > score ) 
		{
		  score = thisScore;
		  bestk = k;
		}


	      //
	      // Clean up
	      //

	      delete model;
	      
	    }

	}



      //////////////////////////////////////////////////////
      //                                                  //
      // Output results, if looking at original data      //
      //                                                  //
      //////////////////////////////////////////////////////
      
      if ( print && par::rtest_type == RTEST_VARIABLE_THRESHOLD ) 
	{
	  
	  int nsnp = useSets ? snpset[scnt].size() : nl_all ;
	  
	  string label = useSets ? setname[scnt] : "ALL";
	  int nsnp2 = 0;
	      
	  double qtmeana = 0, qtmeanu = 0;
	  int cnta = 0 , cntu = 0;
	  
	  for (int l0=0; l0<nsnp; l0++)
	    {
	      
	      int l = useSets ? snpset[scnt][l0] : l0 ;
	      
	      bool display = false;

	      if ( par::rtest_type == RTEST_VANILLA 
		   && locus[l]->freq < par::rtest_freq_threshold )
		display = true;
	      else if ( par::rtest_type == RTEST_FREQWEIGHT )
		display = true;
	      else if ( par::rtest_type == RTEST_VARIABLE_THRESHOLD 
			&& allele_counts[l] <= bestk )
		display = true;

	      if ( display )
		{
		  
		  SUMM2 << setw(20) << label << " "  
			<< setw(par::pp_maxsnp) << locus[l]->name << " ";

		  if ( par::rtest_weights ) 
		    SUMM2 << setw(8) << wgt[l] << " ";

		  if ( par::rtest_type == RTEST_FREQWEIGHT ) 
		    SUMM2 << setw(8) << fwgt[l] << " ";

		  SUMM2 << setw(8) << allele_counts[l] << " "
			<< setw(12) << locus[l]->freq << " ";
		  
		  ++nsnp2;
		  
		  if ( attribs )
		    {
		      map<Locus*,string>::iterator i = 
			attrib.find( locus[l] );
		      if ( i != attrib.end() )
			SUMM2 << i->second;		
		    }
		  SUMM2 << "\n";
		}
	    }
	  
	  vector<RVData>::iterator j = si->begin();
	  set<Individual*> hitinds;
	  
	  while ( j != si->end() )
	    {	 
	      	      
	      // Only look at below-threshold SNPs
	      
	      if ( allele_counts[ j->l ] > bestk )
		{
		  ++j;
		  continue;
		}
	      
	      if ( par::bt )
		{
		  if ( !sample[ j->n ]->missing )
		    {
		      if ( sample[ j->n ]->aff )
			cnta += j->c;
		      else
			cntu += j->c;
		    }
		}
	      else
		{
		  // Individual has 1+ copy 
		  
		  hitinds.insert( sample[ j->n ] );	
		}
	      
	      ++j;
	    }
	  
	  if ( ! par::bt )
	    for (int i=0; i<n; i++)
	      {
		if ( !sample[i]->missing )
		  {
		    if ( hitinds.find( sample[i] ) != hitinds.end() )
		      {
			qtmeana += sample[i]->phenotype;
			++cnta;
		      }
		    else
		      {
			qtmeanu += sample[i]->phenotype;
			++cntu;
		      }		      
		  }
	      }
	  
	  SUMM << setw(20) << label << " "
	       << setw(8) << nsnp << " "
	       << setw(6) << bestk << " "
	       << setw(10) << (double)bestk/(double)(2*n) << " "
	       << setw(8) << nsnp2 << " "; 
	  
	  if  ( par::bt ) 
	    SUMM << setw(8) << cnta << " "
		 << setw(8) << cntu << " ";
	  else
	    SUMM << setw(8) << cnta << " "
		 << setw(8) << cntu << " "
		 << setw(8) << -qtmeana/(double)cnta << " "   // Note: flipping sign of pheno back
		 << setw(8) << -qtmeanu/(double)cntu << " ";  // Note: flipping sign of pheno back
	  
	  SUMM << "\n";	      
	}
      
            

      //////////////////////////////////////////////////////
      //                                                  //
      // Store association statistic for emp. signif.     //
      //                                                  //
      //////////////////////////////////////////////////////

      // If 1-sided/low-score-enrichment test, pass the 
      // negative of the score to the permutation class, 
      // otherwise just pass the score


      if ( par::rtest_hyp == RTEST_1SIDED_LOW )
	results.push_back( -score );	
      else
	results.push_back( score );	


      //////////////////////////////////////////////////////
      //                                                  //
      // Look at the next set of rare variants            //
      //                                                  //
      //////////////////////////////////////////////////////

      ++scnt;
      ++si;
      
    }
  
  
  return results;
}
  
