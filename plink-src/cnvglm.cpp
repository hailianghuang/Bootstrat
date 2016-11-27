

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
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
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "cnv.h"
#include "model.h"

extern Plink * PP;

void Plink::runTestCNVwithGLM(Perm & perm)
{

  // Permutation test for mean difference in QT between people with
  // versus without a CNV.  By default two-sided, unless
  // par::segment_test_force_1sided = T


  ///////////////////////////////////////////////////
  // Set up association model

  bool OLD_assoc_glm_without_main_snp = par::assoc_glm_without_main_snp;
  bool OLD_clist = par::clist;

  par::assoc_glm_without_main_snp = true;
  par::clist = true;
  
  par::clist_number++;
  for (int i=0; i<n; i++)
    sample[i]->clist.resize( par::clist_number );
  
  int testTerm = par::clist_number - 1;
  clistname.resize( par::clist_number );
  clistname[ testTerm ] = "CNV";

  
  //////////////////////////////////////////
  // Test positons = MAP positions (nl_all)

  int nt = nl_all;
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Initialise permutation procedures                              //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  perm.setTests(nt);
  perm.setPermClusters(*this);
  perm.originalOrder();
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Standard positional tests                                      //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  
  if ( par::cnv_indiv_perm )
    error("--cnv-indiv-perm is not implemented --cnv-glm yet");
      
  vector<int> count;

  vector_t original = testCNVwithGLM( true, perm , count );
  
  


  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Run permutations                                               //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  if ( ! par::permute ) 
    return;

  bool finished = false;

  while(!finished)
    {      
      perm.permuteInCluster();
      vector_t pr = testCNVwithGLM(false, perm,count);
      finished = perm.update(pr,original);      
    }
  
  if (!par::silent)
    cout << "\n\n";
  
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Display permuted results                                       //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  string f = par::output_file_name + ".cnv.glm.summary.mperm";
  printLOG("Writing CNV GLM permutation results to [ "+f+" ]\n");

  ofstream FOUT;
  FOUT.open( f.c_str() , ios::out );
  FOUT.precision(4);
  
  FOUT << setw(4) << "CHR" << " "
       << setw(par::pp_maxsnp) << "SNP" << " "
       << setw(12) << "BP" << " "
       << setw(12) << "EMP1" << " "
       << setw(12) << "EMP2" << "\n";
  
  for (int l=0; l<nt; l++)
    {
      FOUT << setw(4) << locus[l]->chr << " "
	   << setw(par::pp_maxsnp) << locus[l]->name << " "
	   << setw(12) << locus[l]->bp << " "
	   << setw(12) << perm.pvalue( l ) << " " 
	   << setw(12) << perm.max_pvalue( l ) << "\n";      
    }
  
  FOUT.close(); 
  
}



vector_t Plink::testCNVwithGLM(bool print, Perm & pperm , vector<int> & counts )
{
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Report to summary file                                         //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  ofstream FOUT;

  if ( print )
    {

      string f = par::output_file_name + ".cnv.glm.summary";
      printLOG("Writing CNV GLM summary to [ "+f+" ]\n");
            
      FOUT.open( f.c_str() , ios::out );
      FOUT.precision(4);
      
      string effectLabel = par::bt ? "OR" : "BETA"; 
      string testLabel = par::seg_test_region ? "GRP" : "POS";
      
      if ( par::seg_test_region )
	{
	  FOUT << setw(4) << "CHR" << " " 
	       << setw(16) << "REGION" << " " 
	       << setw(12) << "BP1" << " "
	       << setw(12) << "BP2" << " "
	       << setw(8) << "NCNV" << " "
	       << setw(12) << effectLabel << " "
		   << setw(12) << "P" << "\n";	  
	}
      else
	{
	  FOUT << setw(4) << "CHR" << " "
	       << setw(par::pp_maxsnp) << testLabel << " "
	       << setw(12) << "BP" << " "
	       << setw(8) << "NCNV" << " "
	       << setw(12) << effectLabel << " "
	       << setw(12) << "P" << "\n";
	}
    }
  
    
  vector_t score;

  int testTerm = par::clist_number - 1;
  cout << "testTerm = " << testTerm << "\n";

  if ( par::seg_test_region )
    {
      
      //////////////////////////////////
      // Testing regions, not positions
      
      int nGenes = geneList.size();
      score.resize( nGenes, 0 );
      counts.resize( nGenes, 0 );

      int l = 0;

      set<Range>::iterator g = geneList.begin();
      
      while ( g != geneList.end() )
	{

	  // Clear covariate slot
	  for (int i=0; i<n; i++)
	    sample[i]->clist[ testTerm ] = 0;
	  
	  Range tr = *g;
	  
	  // Either calculate, of if already done, just lookup the CNVs that
	  // fall in this range (i.e. this function will be called many times
	  // when using permutation)
	  
	  map<Range, set<Segment> >::iterator i = gene2segment.find( tr );
	  if ( i == gene2segment.end() )
	    {
	      set<Segment> s = allSegmentsIntersecting( tr );
	      gene2segment.insert(make_pair( tr, s ) );
	      i = gene2segment.find( tr );
	    }
	  
	  set<Segment>::iterator si = i->second.begin();
	  while ( si != i->second.end() )
	    {
	      if ( si->p1->pperson->missing )
		continue;
	      si->p1->pperson->clist[ testTerm ] += 1;		  
	      counts[l]++;
	      ++si;
	    }
	  
	  
	  // Perform test

	  glmAssoc( false , pperm );
	  
	  bool valid = model->isValid();
	  vector_t b = model->getCoefs();
	  vector_t pvs = model->getPVals();
	  double pv = pvs[ testTerm ];
	  double chisq = model->getStatistic();
	  double z=0;
	  double unreal = 1/z;
	  score[l] = valid ? chisq : unreal ;
	  double beta = par::bt ? exp(b[testTerm+1]) : b[testTerm+1] ;

	  // A one-sided test?	  
 	  if ( par::segment_test_force_1sided && valid )
 	    {	      
 	      if ( b[ testTerm + 1 ] < 0 ) 
 		{
		  score[ l ] = 0;
		  pv = 1;
		}
 	    }
	  
	  if ( print )
	    {
	      FOUT << setw(4) << g->chr << " " 
		   << setw(16) << g->name << " " 
		   << setw(12) << g->start << " "
		   << setw(12) << g->stop << " "
		   << setw(8) << counts[l] << " "
		   << setw(12) << beta << " "
		   << setw(12) << pv << "\n";
	    }
 

	  // Clean up
	  
	  delete model;
	  
	  ++g;
	  ++l;
	}
    }
  else  
    {
	  
      ////////////////////////////////////////
      // Consider each position at a time
      
      score.resize( nl_all, 0 );
      counts.resize( nl_all , 0 );

      for (int l=0; l<nl_all; l++)
	{
	  
	  // Clear covariate slot
	  
	  for (int i=0; i<n; i++)
	    sample[i]->clist[ testTerm ] = 0;
	  
	  // Count actual segments
	  
	  vector<Segment>::iterator s = segment.begin();
	  while ( s != segment.end() )
	    {    		  
	      if ( l >= s->start && l <= s->finish  ) 
		{
		  s->p1->pperson->clist[ testTerm ] += 1;
		  counts[l]++;
		}
	      s++;
	    }	
	  
	  // Allow for optional 'wings'
	      
	  if ( par::seg_test_window )
	    {		  
	      vector<Segment>::iterator s = segment.begin();
	      while ( s != segment.end() )
		{    	      
		  // Shift left from start
		  int l = s->start;
		  Locus * loc1 = locus[s->start];
		  while ( 1 )
		    {
		      --l;
		      if ( l < 0 ) break;
		      Locus * loc2 = locus[l];
		      if ( loc2->chr != loc1->chr ) break;
		      if ( loc1->bp - loc2->bp > par::seg_test_window_bp  ) break;		
		      s->p1->pperson->clist[ testTerm ] += 1;
		      counts[l]++;
		    }
		  
		  // Shift right from start
		  l = s->finish;
		  loc1 = locus[s->finish];
		  while ( 1 )
		    {
		      ++l;
		      if ( l == nl_all ) break;
		      Locus * loc2 = locus[l];
		      if ( loc2->chr != loc1->chr ) break;
		      if ( loc2->bp - loc1->bp > par::seg_test_window_bp  ) break;
		      s->p1->pperson->clist[ testTerm ] += 1;
		      counts[l]++;
		    }
		  
		  // Next segment
		  s++;
		}  
	    }
	  
	  
	  // Perform test
	  
	  glmAssoc( false , pperm );
	  
	  bool valid = model->isValid();
	  vector_t b = model->getCoefs();
	  vector_t pvs = model->getPVals();
	  double pv = pvs[ testTerm ];
	  double chisq = model->getStatistic();
	  double z=0;
	  double unreal = 1/z;
	  score[l] = valid ? chisq : unreal ;

	  double beta = par::bt ? exp(b[testTerm+1]) : b[testTerm+1] ;

// 	  cout << valid << "\t";
// 	  display(b);
// 	  display(pvs);
// 	  cout << "-----\n";

	  // A one-sided test?	  
 	  if ( par::segment_test_force_1sided && valid )
 	    {	      
 	      if ( b[ testTerm + 1 ] < 0 ) 
 		{
		  score[ l ] = 0;
		  pv = 1;
		}
 	    }
	  
	  if ( print )
	    {
	      FOUT << setw(4) << locus[l]->chr << " "
		   << setw(par::pp_maxsnp) << locus[l]->name << " "
		   << setw(12) << locus[l]->bp << " "
		   << setw(8) << counts[l] << " "
		   << setw(12) << beta << " "
		   << setw(12) << pv << "\n";
	    }
   

	  // Clean up
	  
	  delete model;
	  
	} // Next position
      
    }


  if ( print )
    FOUT.close();
      
  return score;
  
}
 


