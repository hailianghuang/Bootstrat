

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2010 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <vector>
#include <iomanip>

#include "plink.h"
#include "phase.h"
#include "options.h"
#include "helper.h"
#include "stats.h"
#include "model.h"
#include "whap.h"
#include <cmath>
#include "linear.h"
#include "logistic.h"

extern Plink * PP;

bool mytest(int t, 
	    bool allele, 
	    vector<int> & snps, 
	    double & lowestAIC, 
	    vector<set<int> > & bestg, 
	    int & spechap,
	    double origfrq,
	    ofstream * FOUT );

string writeSNPs(vector<int> & snps, string d=",")
{
  string str;
  for (int s = 0 ; s < snps.size(); s++)
    {
      str += PP->locus[ snps[s] ]->name;
      if ( s < snps.size()-1 ) str += d;
    }
  return str;
}

bool testHap()
{
  PP->glmAssoc(false,*(PP->pperm));
  return PP->model->isValid();
}

bool testHaplotypesTwoHap()
{
  par::assoc_glm_without_main_snp = true;
  par::proxy_glm = true;
  par::test_hap_GLM = false;
  par::chap_test = false;
  return testHap();  
}

// Alternate wrapper: use WHAP coding
bool testHaplotypesWHAP(vector< set<int> > & grp)
{

  Chap thisCModel(PP, PP->haplo);
  
  ChapModel myModel;
  myModel.group = grp;  

  thisCModel.setModels(myModel,myModel);
  PP->whap = & thisCModel;
  PP->whap->current = & myModel;
  
  par::assoc_glm_without_main_snp = true;
  par::proxy_glm = false;
  par::test_hap_GLM = false;
  par::chap_test = true;
  
  return testHap();  

}

// Alternate wrapper, that allows for >2 haplotypes (i.e. not proxy framing)
bool testHaplotypesFull()
{
  par::assoc_glm_without_main_snp = true;
  par::proxy_glm = false;
  par::test_hap_GLM = true;
  par::chap_test = false;
  return testHap();
}


vector_t coefs()  { return PP->model->getCoefs(); }
double statistic() { return PP->model->getStatistic(); }
double pval(int df) { return chiprobP( PP->model->getStatistic() , df );  }
double pval(double stat, int df) { return chiprobP( stat , df );  }
int nump() { return PP->model->getNP(); }

double coef(int i) 
{
  vector_t c = PP->model->getCoefs();
  return i >= c.size() ? 0 : c[i];
}

double loglik()
{
  // Note: getLnLk() returns -2LL
  LogisticModel * lmodel = (LogisticModel*)PP->model;
  return lmodel->getLnLk();
}

double AIC(double ll, int np)
{
  return 2 * np + ll;
}

double AIC()
{
  return 2 * nump() + loglik();
}

void freeModel()
{
  if ( PP->model )
    delete PP->model;
}



/////////////////////////////////////////////////
// Main wrapper for LF-SCAN

void Plink::lfsearch( string s )
{
  
  int ref = -1;
  for (int l=0; l<nl_all; l++)
    if ( locus[l]->name == s )
      {
	ref = l;
	break;
      }
  
  if ( ref < 0 ) 
    error("Cannot find target SNP [ " + s + " ]\n");


  // Determine which allele is the risk increasing one (i.e. 
  // assume we are looking for rare disease alleles)

  vector<int> dummy;
  dummy.push_back( ref );

  haplo->phaseSNPs( dummy );
    

  // Get direction of effect for target SNP
  // to determine risk allele
  
  haplo->testSet.clear();
  for (int h2=0; h2 < PP->haplo->nh; h2++)
    {
      if ( PP->haplo->hap[h2][0] )
	PP->haplo->testSet.insert(make_pair(h2,0));
      else
	PP->haplo->testSet.insert(make_pair(h2,1));
    }
  

  // Perform regression of disease on haplotypes

  bool okay = testHaplotypesTwoHap();
  
  // Get results
  
  if ( ! par::bt )
    error("This routine currently only set up for disease traits");

  if ( ! okay ) 
    error("Could not perform basic single SNP test... quitting");

  double b      = exp( coef(1) );
  double m1_aic = AIC( loglik(), nump() );
  double pv     = pval(1);


  // Clean up association model

  freeModel();



  // Determine the high-risk allele
  
  bool allele = b < 1 ;
  if ( allele ) b = 1.0/b;  
  double frq = 0;
  string allelename;
  
  for (int h=0; h<haplo->nh; h++)
    {
      if ( allele )
	{
	  if ( !haplo->hap[h][0] ) 
	    {
	      allelename = haplo->haplotypeName(h);
	      frq = haplo->f[h];
	    }
	}
      else
	{
	  if ( haplo->hap[h][0] ) 
	    {
	      allelename = haplo->haplotypeName(h);
	      frq = haplo->f[h];
	    }
	}      
    }
  




  ///////////////////////////////////////////////////////////////
  // Consider all surrounding haplotypes, exhaustive search
    
  // Consider haplotypes up to N SNPs (excluding target)

  int minsnp = 3;
  int maxsnp = 3;

  // For now, we will search *all* SNPs in dataset
  // i.e. necessary to use --snp and --window or 
  //  some such
  
  int search_cnt = nl_all - 1 ;
  
  if ( search_cnt > 1000 ) 
    error("Too large a search space: reduce --window");

  vector<int> mapback;
  for (int z=0; z < nl_all ; z++)
    if ( z != ref )
      mapback.push_back(z);
  

  // Use these to save the best model and 
  // keep track of progress

  vector<set<int> > bestg;
  double myAIC = m1_aic;
  vector<int> bestSet;
  bool found1 = false;
  int ntests = 0;
  int spechap = -1;
  
  for ( int i=minsnp ; i<=maxsnp; i++)
    {
      
      // For an i-SNP of cnt-SNP, consider all the permutations
      
      vector<unsigned long> pos1(i);
      vector<int> d( search_cnt );

      for (int z=0; z<search_cnt; z++) 
	d[z] = z;
      
      vector<vector<int> > collection;

      combinations_recursive(d,i,pos1,0,0,collection);
      
      if ( ! par::silent )
	cerr << "Searching all " << collection.size() << " haplotypes of size " << i << "\n";
      
      ntests += collection.size();

      for (int c1 = 0; c1 < collection.size() ; c1++ )
	{
	  
	  if ( ! par::silent ) 
	    cerr << c1 << " tests performed          \r";

	  vector<int> posit;

	  // Add target SNP first
	  
	  posit.push_back( ref );
	  
	  // And now the tagging SNPs

	  for (int k=0; k<collection[c1].size(); k++)
	    posit.push_back( mapback[ collection[c1][k] ] );


	  // Phase and perform series of haplotype tests

	  if ( mytest( ref, allele, posit , myAIC, bestg, spechap, frq, NULL ) )
	    {
	      bestSet  = posit;
	      found1 = true;
	    }
	  
	}

    } // Next haplotype size
  

  if ( ! par::silent ) 
    cerr << "\n";


  printLOG("Tested " + int2str(ntests) + " haplotypes in total\n");

  if ( ! found1 ) 
    {
      printLOG("LF-scan did not yield any potential haplotypes\n");
      return;
    }

  // Report results for this single SNP

  ofstream FOUT;
  printLOG("Writing LF-scan results to [ " + par::output_file_name + "." + s + ".lfscan ]\n");
  FOUT.open( (par::output_file_name+"." + s + ".lfscan").c_str() , ios::out );
  FOUT.precision(3);  
  FOUT.setf(ios::fixed);
  
  FOUT << "\n";
  FOUT << setw(10) << "TEST" << " "
       << setw(8) << "ALLELE" << " "
       << setw(8) << "FRQ" << " "
       << setw(8) << "OR" << " "
       << setw(12) << "AIC" << " "
       << setw(10) << "P" << "  "
       << "SNPS" << "\n";
  
  FOUT << setw(10) << "TARGET" << " " 
       << setw(8) << allelename << " "
       << setw(8)  << frq << " "
       << setw(8) << b << " "
       << setw(12) << m1_aic << " ";
  FOUT.setf( ios::scientific );
  FOUT << setw(10) << pv << "  ";
  FOUT.unsetf( ios::scientific );
  FOUT << locus[ ref ]->name << "\n"; 

  printLOG("LF-scan found a potential haplotype\n");


  // Repeat analyses and show results for best set

  mytest( ref, allele, bestSet , myAIC, bestg, spechap, frq, &FOUT );
    

  // All done

  FOUT.close();
  return;
}



bool mytest(int t,
	    bool allele,
	    vector<int> & snps,
	    double & myAIC,
	    vector<set<int> > & bestg,
	    int & spechap,
	    double frq,
	    ofstream * FOUT)
{
  
  
  // Apply all tests to the set of SNPs given by 'snps'. We expect the
  // target/reference SNP to be the first entry here
  
  // Do we find a better model, based on lowest AIC?

  bool foundModel = false;


  // Phase these markers w/ EM
  
  PP->haplo->phaseSNPs( snps );
  

  // Count number of common haplotypes
  // and SNPs
  
  int nch = 0;  
  set<int> commonHaplotypes;
  set<int> riskBackground;
  for (int h=0; h < PP->haplo->nh; h++)
    {
      
      if ( allele != PP->haplo->hap[h][0] )
	riskBackground.insert(h);

      if ( PP->haplo->f[h] >= par::min_hf )
	{
	  ++nch;
	  commonHaplotypes.insert(h);
	}
    }
  int ns = PP->haplo->ns;
  

  // Search for a lower freq. haplotype 
  
  map<string,int> hap1;    // Found on high-risk background
  set<string> hap2;        // Found on low-risk background


  // Consider each haplotype, is it a potential candidate?
  
  for (int h=0; h< PP->haplo->nh; h++)
    {
      
      // On background?
      
      if ( allele != PP->haplo->hap[h][0] )
	{
	  if ( PP->haplo->f[h] > par::min_hf ) 
	    {
	      hap1.insert(make_pair( PP->haplo->haplotypeName(h) , h ) );		      
	    }
	}
      else
	{
	  
	  // Hard-coded EPS value for haplotype frequency on non-risk
	  // allele
	  
	  if ( PP->haplo->f[h] > 0.0001 ) 
	    {
	      hap2.insert( PP->haplo->haplotypeName(h).substr(1,ns-1) );
	    }
	}
      
    }
  

  // Get list haplotypes that are "unique" to the high-risk 
  // common allele. These are potential candidates to test
  // for stronger association, i.e. if the common alelle 
  // associaton in fact indexes a lower frequency, stronger
  // effect

  set<int> fnd;        
  
  for (int h=0; h< PP->haplo->nh; h++)
    {
      
      if ( PP->haplo->f[h] < par::min_hf )
	continue;
      
      string subname = PP->haplo->haplotypeName(h).substr(1,ns-1);
      string name = PP->haplo->haplotypeName(h);

      map<string,int>::iterator i0 = hap1.find( name );
      
      if ( i0 != hap1.end() && hap2.find( subname ) == hap2.end() )
	{ 	  
	  fnd.insert(h);
	}

    }


 
  /////////////////////////////////////////
  // Quit, if we did not find any potential 
  // candidates

  if ( fnd.size() == 0 ) 
    return foundModel;


  // Get AIC for primary model, using this exact phasing
  // (i.e. just if misining data in the target SNP, etc, check
  // that we are comparing like with like
  
  PP->haplo->testSet.clear();
  
  for (int h2=0; h2 < PP->haplo->nh; h2++)
    {
      // Note: direction of effect doesn't matter here
      PP->haplo->testSet.insert(make_pair(h2, (int)(PP->haplo->hap[h2][0]) ));
    }
  

  // Test haplotype (i.e. reconstitute single SNP test within this
  // haplotypic framework), e.g.

  //  A-----
  //  A-----
  //    vs.
  //  B-----
  //  B-----
  //  B-----
  

  testHaplotypesTwoHap();
  

  // Only need AIC and -2LL from here

  double m1_aic = AIC();
  double m1_loglik = loglik();
  int    m1_np     = nump();

    
  // Clean up

  freeModel();

  


  ////////////////////////////////////
  // Test alternate haplotypic models
  
  PP->haplo->sets.clear();
  

  // Find a reference haplotype that is the first 
  // non-risk hi-freq allele, i.e. just to make 
  // the ORs in the omnibus model nicer
  
  int reference = 0;
  const int tsnp = 0;

  for (int h=0; h<PP->haplo->nh; h++)
    {
      if ( allele == PP->haplo->hap[h][tsnp] &&
	   fnd.find( h ) == fnd.end() &&
	   commonHaplotypes.find(h) != commonHaplotypes.end() ) 
	{
	  reference = h;
	  break;
	}
    }
  


  //////////////////////////////////
  // Fit omnibus model

  map<int,int> omnimap;
  int cnter = 1; // start at 1, to ignore intercept

  set<int>::iterator i = commonHaplotypes.begin();
  while ( i != commonHaplotypes.end() )
    {
      if ( *i != reference ) 
	{
	  PP->haplo->sets.insert(*i);
	  omnimap.insert(make_pair( *i , cnter++ ) );
	}
      ++i;
    }
  
  
  // Use the full wrapper, as >2 haplotypes being tested
  
  testHaplotypesFull();
  

  // Get -2LL and # of parameters under omnibus model

  double omni_loglik = loglik();
  int    omni_np     = nump();
  double omni_aic    = AIC(omni_loglik,omni_np);
  vector_t omni_beta = coefs();
  

  // Save a best AIC? 
  
  if ( omni_aic < m1_aic && omni_aic < myAIC )
    {
      myAIC = omni_aic;

      spechap = -1; // Omnibus model

      // save omnibus model
      bestg.clear();
      for (int h=0; h<PP->haplo->nh; h++)
	{
	  if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
	    {
	      set<int> t;
	      t.insert(h);
	      bestg.push_back(t);
	    }
	}

      foundModel = true;
    }


  
  // Clean up
  
  freeModel();
  
      

  ///////////////////////////////////////////////////////////////////
  // Omnibus submodel 1: constrain all haplotypes on risk background

  // Specify tests using WHAP coding

  vector<set<int> > g(1);
  
  i = commonHaplotypes.begin();
  
  while ( i != commonHaplotypes.end() )
    {
      if ( riskBackground.find( *i ) != riskBackground.end() )
	g[0].insert( *i );
      else
	{
	  set<int> t;
	  t.insert( *i );
	  g.push_back(t);
	}
      ++i;
    }
  
  
  // Use the full wrapper, as >2 haplotypes being tested
  
  testHaplotypesWHAP(g);

  // Get -2LL and # of parameters under omnibus model

  double omni_s1_loglik = loglik();
  int    omni_s1_np     = nump();
  double omni_s1_aic    = AIC(omni_s1_loglik,omni_s1_np);


  // Best model yet?

  if ( omni_s1_aic < m1_aic && omni_s1_aic < myAIC )
    {
      myAIC = omni_s1_aic;
      bestg = g;
      spechap = -2; // S1 Omnibus model
      foundModel = true;
    }


  // Clean up
  
  freeModel();

  

  ///////////////////////////////////////////////////////////////////////
  // Omnibus submodel 2: constrain all haplotypes on non-risk background
  
  // Specify tests using WHAP coding

  g.clear();
  set<int> tmp;
  g.push_back(tmp);

  i = commonHaplotypes.begin();
  
  while ( i != commonHaplotypes.end() )
    {
      if ( riskBackground.find( *i ) == riskBackground.end() )
	g[0].insert( *i );
      else
	{
	  set<int> t;
	  t.insert( *i );
	  g.push_back(t);
	}
      ++i;
    }
  
  
  // Use the full wrapper, as >2 haplotypes being tested

  testHaplotypesWHAP(g);


  // Get -2LL and # of parameters under omnibus model

  double omni_s2_loglik = loglik();
  int    omni_s2_np     = nump();
  double omni_s2_aic    = AIC(omni_s2_loglik,omni_s2_np);
  
  // Best model yet?

  if ( omni_s2_aic < m1_aic && omni_s2_aic < myAIC )
    {
      myAIC = omni_s2_aic;
      spechap = -3; // S2 Omnibus model  
      bestg = g;
      foundModel = true;
    }
  
  
  // Clean up
  
  freeModel();

  

  ////////////////////////////////////////////////////
  // Now consider all other haplotype-specific models
      
  set<int>::iterator i0 = fnd.begin();
      
  while ( i0 != fnd.end() )
    {

      // Test each haplotype
      	  
      PP->haplo->testSet.clear();
      
      for (int h2=0; h2 < PP->haplo->nh; h2++)
	{
	  if ( *i0 == h2 )
	    PP->haplo->testSet.insert(make_pair(h2,0));
	  else
	    PP->haplo->testSet.insert(make_pair(h2,1));
	}
      

      // Test haplotype
      
      testHaplotypesTwoHap();


      // Get results
      
      double m2_odds = exp(coef(1));
      double m2_pval = pval(1);
      double m2_aic  = AIC();
      
      
      // Clean up
      
      freeModel();
	  

      // Does this potential haplotype better explain the 
      // original association?  Assess based on AIC
	  
      if ( m2_aic < myAIC && m2_aic < m1_aic ) 
	{	  
	  foundModel = true;
	  myAIC = m2_aic;

	  spechap = *i0;

	  bestg.clear();
	  set<int> t;
	  for (int h=0; h<PP->haplo->nh; h++)
	    {
	      if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
		if ( h != *i0 )
		  t.insert(h);	      
	    }
	  bestg.push_back(t);
	  t.clear();
	  t.insert( *i0 );
	  bestg.push_back(t);

	}
          
      // Next root of interest	  
      ++i0;
	  
    }


  /////////////////////////////////////////////////////////
  // If we have a non-null pointer to an ofstream,
  // this must mean we've a) finished, b) have something
  // to report. The best model will be in bestg

  
  if ( FOUT ) 
    {
      
      *FOUT << "\nMarkers (target first): " << writeSNPs( snps ) << "\n\n";

      *FOUT << setw(12) << "HAP" << " " 
	    << setw(8)  << "FRQ" << " " 
	    << setw(8)  << "OR" << " "
	    << setw(6)  << "FLAG" << "\n";

      *FOUT << setw(12) << "------------" << " " 
	    << setw(8)  << "--------" << " " 
	    << setw(8)  << "--------" << " "
	    << setw(6)  << "------" << "\n";


      // Order as High-risk; candidate/non-candidate; low-risk

      vector<int> ord;

      int last1 = -1;
      int last2 = -1;
  
      for (int h=0; h< PP->haplo->nh; h++)
	if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
	  if ( riskBackground.find(h) != riskBackground.end() )
	    if ( fnd.find(h) != fnd.end() )
	      {
		ord.push_back(h);
		last1 = h;
	      }

      for (int h=0; h< PP->haplo->nh; h++)
	if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
	  if ( riskBackground.find(h) != riskBackground.end() )
	    if ( fnd.find(h) == fnd.end() )
	      {
		ord.push_back(h);
		last2 = h;
	      }

      for (int h=0; h< PP->haplo->nh; h++)
	if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
	  if ( riskBackground.find(h) == riskBackground.end() )
	    ord.push_back(h);
	    
     
      for (int h2=0; h2< ord.size(); h2++)
 	{
 	  int h = ord[h2];

	  string name = PP->haplo->haplotypeName(h).substr(1,ns-1);
	  string fullname = PP->haplo->haplotypeName(h);
	  string label = fullname.substr(0,1) + "-" + name;

	  *FOUT << setw(12) << label << " "  
		<< setw(8)  << PP->haplo->f[h] << " ";
	  
	  map<int,int>::iterator i1 = omnimap.find( h );
	  if ( i1 == omnimap.end() )
	    *FOUT << setw(8) << "-ref-" << " ";
	  else
	    *FOUT << setw(8) << exp(omni_beta[ i1->second ] ) << " "; 
	  
	  map<string,int>::iterator i0 = hap1.find( fullname );
	  if ( i0 != hap1.end() && hap2.find( name ) == hap2.end() )
	    *FOUT << setw(6) << "*";
	  else
	    *FOUT << setw(6) << "-";
	  
	  *FOUT << "\n";
	  
	  if ( h == last1 || h == last2 ) 
	    *FOUT << "\n";

	}
      *FOUT << "\n";



      ////////////////////////////////////////////////////
      // Recalculate association under best-fitting model

      testHaplotypesWHAP( bestg );
      

      // Get -2LL and # of parameters under omnibus model
      
      double best_loglik = loglik();
      int    best_np     = nump();
      double best_aic    = AIC();
      vector_t best_beta = coefs();
      
      // Clean up

      freeModel();



      ////////////////////////////////////////////////////
      // Fit completely null model
      
      vector<set<int> > nullg;
      set<int> t;
      for (int h=0; h<PP->haplo->nh; h++)
	{
	  if ( commonHaplotypes.find(h) != commonHaplotypes.end() )
	    t.insert(h);	  
	}
      nullg.push_back(t);

      testHaplotypesWHAP( nullg );

      double null_loglik = loglik();
      int    null_np     = nump();
      double null_aic    = AIC();

      // Clean up

      freeModel();



      ////////////////////////////////////////////////////
      // Report

      *FOUT << setw(12) << "MODEL" << " "
	    << setw(12) << "-2LL" << " "
	    << setw(5) << "NP" << " "
	    << setw(12) << "AIC" << " "
	    << setw(12) << "LRT(1)" << " "
	    << setw(12) << "LRT(2)" << "\n";
	
      *FOUT << setw(12) << "------------" << " "
	    << setw(12) << "------------" << " "
	    << setw(5) << "-----" << " "
	    << setw(12) << "------------" << " "
	    << setw(12) << "------------" << " "
	    << setw(12) << "------------" << "\n";

      *FOUT << setw(12) << "NULL" << " "
	    << setw(12) << null_loglik << " "
	    << setw(5) << null_np << " "
	    << setw(12) << null_aic << " "
	    << setw(12) << "NA" << " "
	    << setw(12) << "NA" << "\n";
      
      *FOUT << setw(12) << "SNP" << " "
	    << setw(12) << m1_loglik << " "
	    << setw(5) << m1_np << " "
	    << setw(12) << m1_aic << " ";
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(null_loglik - m1_loglik,m1_np - null_np) << " ";
      FOUT->unsetf(ios::scientific);

      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval( m1_loglik - omni_loglik , omni_np - m1_np ) << "\n";
      FOUT->unsetf(ios::scientific);


      *FOUT << setw(12) << "OMNI" << " "
	    << setw(12) << omni_loglik << " "
	    << setw(5) << omni_np << " "
	    << setw(12) << omni_aic << " ";
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(null_loglik - omni_loglik,omni_np - null_np) << " ";
      FOUT->unsetf(ios::scientific);      
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << "NA" << "\n";
      FOUT->unsetf(ios::scientific);

      *FOUT << setw(12) << "SUB1" << " "
	    << setw(12) << omni_s1_loglik << " "
	    << setw(5) << omni_s1_np << " "
	    << setw(12) << omni_s1_aic << " ";
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(null_loglik - omni_s1_loglik,omni_s1_np - null_np) << " ";
      FOUT->unsetf(ios::scientific);

      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval( omni_s1_loglik - omni_loglik, omni_np - omni_s1_np) << "\n";
      FOUT->unsetf(ios::scientific);

      *FOUT << setw(12) << "SUB0" << " "
	    << setw(12) << omni_s2_loglik << " "
	    << setw(5) << omni_s2_np << " "
	    << setw(12) << omni_s2_aic << " ";
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(null_loglik - omni_s2_loglik,omni_s2_np - null_np) << " ";
      FOUT->unsetf(ios::scientific);

      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval( omni_s2_loglik - omni_loglik, omni_np - omni_s2_np) << "\n";
      FOUT->unsetf(ios::scientific);


      *FOUT << "\n";

      *FOUT << setw(12) << "BEST" << " "
	    << setw(12) << best_loglik << " "
	    << setw(5) << best_np << " "
	    << setw(12) << best_aic << " ";
      
      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(null_loglik - best_loglik,best_np - null_np) << " ";
      FOUT->unsetf(ios::scientific);

      FOUT->setf(ios::scientific);
      *FOUT << setw(12) << pval(best_loglik - omni_loglik, omni_np - best_np) << "\n";
      FOUT->unsetf(ios::scientific);


      *FOUT << "\n\n";



      *FOUT << setw(8) << "RISKGRP" << " " 
	    << setw(8) << "OR" << " " 
	    << setw(8) << "HAP" << " "
	    << setw(8) << "FRQ" << "\n";

      *FOUT << setw(8) << "--------" << " " 
	    << setw(8) << "--------" << " " 
	    << setw(8) << "--------" << " " 
	    << setw(8) << "--------" << "\n"; 

      for (int hg=0; hg<bestg.size(); hg++)
	{

	  *FOUT << setw(8) << hg+1 << " ";
	  if ( hg == 0 ) 
	    *FOUT << setw(8) << "-ref-" << " ";
	  else
	    *FOUT << setw(8) << exp(best_beta[hg]) << " ";
	  bool done1 = false;

	  set<int>::iterator i = bestg[hg].begin();
	  while ( i != bestg[hg].end() )
	    {
	      
	      // either reference, or hg+1 (i.e. ignore intercept)

	      if ( commonHaplotypes.find( *i ) != commonHaplotypes.end() )
		{
		  string fullname = PP->haplo->haplotypeName( *i );
		  string label = fullname.substr(0,1)+ "-" + fullname.substr(1,ns-1);
		  
		  if ( ! done1 ) 
		    {
		      *FOUT << setw(8) << label << " "
			    << setw(8) << PP->haplo->f[ *i ] << "\n";
		      done1 = true;
		    }
		  else
		    {
		      *FOUT << setw(8) << " " << " "
			    << setw(8) << " " << " "
			    << setw(8) << label << " "
			    << setw(8) << PP->haplo->f[ *i ] << "\n";
		    }
		}

	      ++i;
	    }

	  *FOUT << "\n";
	 
	}
      
      *FOUT << "\nDIFF_AIC " << m1_aic - best_aic << "\n";
      
      if ( spechap >= 0 )
	{
	  *FOUT << "MODEL SPEC_HAP " << PP->haplo->haplotypeName(spechap) << "\n"
		<< "DIFF_FRQ " << frq << " " << PP->haplo->f[spechap] << "\n";
	}
      else if ( spechap == -1 ) 
	*FOUT << "MODEL OMNI\n";
      else if ( spechap == -2 )
	*FOUT << "MODEL SUB1\n";
      else if ( spechap == -3 )
	*FOUT << "MODEL SUB0\n";
      else
	*FOUT << "MODEL UNKNOWN\n";
            
    }

  return foundModel;
  
}
