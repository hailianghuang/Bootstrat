

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __HAPPHASE_H_
#define __HAPPHASE_H__


class HaploPhase
{
 public:

  // Reference to main class
  Plink & P;
  
  int ns;         // Number of SNPs in haplotype
  int nh;         // Number of possible haplotypes
  int np;         // Number of phases, diploid
  int haploid_np; // As above, haploid
  int test_hap;   // To be imputed haplotype
  string hname;   // Name of haplotype locus
  int cnt_f;      // Number of founders to be phased
  int current;    // Number of current haplotype being tested

  vector<vector<bool> > hap; // Coding for each haplotype
  vector<int> S;             // List of predictor SNP numbers

  // Haploid chromosome codes
  bool X;
  bool haploid;
  
  // Estimated haplotype frequencies
  vector<double> f;
  
  // Individual posterior probabilities
  vector<vector<double> > pp;
  
  // Individual x phase 
  vector<vector<int> > hap1;
  vector<vector<int> > hap2;
  
  // Phase markers and frequencies
  vector<int> ph_hap1;
  vector<int> ph_hap2;
  vector<double> ph_freq;

  vector<int> haploid_ph_hap1;
  vector<double> haploid_ph_freq;

  // Lookup table for haplotype number given SNPs
  map<vector<bool>,int> hapmap;
  
  // Whether individual has ambiguous phase
  // i.e. hap1[i].size() == 1 
  vector<bool> ambig;
  
  // Transmission/untransmission counts
  vector<double> trans;
  vector<double> untrans;

  // Should we skip this person?
  vector<bool> include;

  // Tag SNP imputation
  vector<vector<int> >  new_pred_locus;
  vector<string>        new_pred_allele;
  vector<map<string,double> > new_pred_weighted_allele;

  vector<Locus*>        new_map;
  vector<vector<bool> > new_one;
  vector<vector<bool> > new_two;

  // Output files: haplotype frequencies
  ofstream HFRQ;     
  ofstream HTEST;     
  ofstream HPHASE;     
  ofstream WGT;

  // Temporary storage for chi-sqs from haplotype tests
  // and odds ratio (haplotype specific tests)
  double result;
  double pvalue;
  double odds;

  HaploPhase(Plink & P_) : P(P_)  
    {
      ambig.resize(P.n,false);
      include.resize(P.n,true);
      pp.resize(P.n);
      hap1.resize(P.n);
      hap2.resize(P.n);
      X=haploid = false;
    }  

  // Read list of tests/tags
  void readTagFile();

  // Make sliding window list of tests
  void makeSlidingWindow(int);

  // Make set of local haplotype proxies
  void makeProxySet(int);

  // Display haplotype frequencies
  void calculateHaplotypeFrequencies();

  // Make test set for haplotype tests
  map<int,int> makeTestSet(boolvec_t &, boolvec_t &);

  // Return subhaplotype name, formatted
  string getSubHaplotypeName(boolvec_t &, boolvec_t &, int);

  // Perform haplotype tests
  vector_t performHaplotypeTests();

  // Impute all haplotypes
  void imputeAllHaplotypes();
  
  // Display haplotype phases
  void calculatetHaplotypePhases();


  void reset()
    {
      ns = nh = np = 0;
      test_hap = -1;
      for (int i=0; i<hap.size(); i++) hap[i].clear();
      hap.clear();
      S.clear();
      f.clear();
      validN = 0;
      for (int i=0; i<pp.size(); i++)
	pp[i].clear();
      for (int i=0; i<hap1.size(); i++)
	hap1[i].clear();
      for (int i=0; i<hap2.size(); i++)
	hap2[i].clear();
			
      pp.resize(P.n);
      hap1.resize(P.n);
      hap2.resize(P.n);
      ambig.resize(P.n,false);
      include.resize(P.n,true);
    } 
  
  void name(string n) { hname = n; } 
  
  string haplotypeName(int i)
    {
      string str;
      for (int s=0; s<ns; s++)
	{
	  if ( i == -1 ) str += "-"; // haploid gap
	  else if ( hap[i][s] ) str += P.locus[S[s]]->allele1;
	  else str += P.locus[S[s]]->allele2;	      
	}
      return str;
    }


  int nHap() { return nh; } 

  double testHaplotypeFreq() { if (test_hap>=0) return f[test_hap]; else return -1; }  
   
  // Main routine to driver phasing
  void phaseAllHaplotypes();

  // Given list of SNP numbers, set up all possible haplotypes (hap)
  void  enumerateHaplotypes(vector<int>&);
  
  // Give list of haplotypes in each phase
  void enumerateAllPhases();

  // Set test haplotype (query hap with allele string)
  void  setTestHaplotype(string);

  // Return possible haplotype list
  vector<string> returnHaplotypes(vector<int>&);
  
  // Determine possible haplotype phases for an individual
  void enumeratePhase(int);

  // Get rid of unlikely phases
  void prunePhase(int);

  // Phase non-founders, given we have haplotype frequencies
  // and score/rescore for TDT
  void phaseAndScoreNonfounder(int);

  // Use offspring information to help resolve parental phase
  void resolveWithKids(int);

  // Report haplotype phase for an individual
  void reportPhase();

  // Report haplotype phase for an individual, alternate format
  void reportPhaseWideFormat();
  

  // E-M algorithm
  void performEM();
  void performEM_original();

  // Convenience function to return r^2 given two SNPs
  double rsq(int,int);

  // Return pairwise r^2
  double rsq_internal(int,int);
  double rsq_internal(boolvec_t &,
		      boolvec_t &,
		      boolvec_t &,
		      boolvec_t &);

  double freq(boolvec_t &, boolvec_t &);

  
  void imputeThisHaplotype(int);
  double imputeHaplotypes(int,bool&,bool&);
  vector_t imputeGenotype(int,int);

  void reportHaplotypeFrequencies();

  void haplotypicCC(map<int,int> &, int,bool);
  void haplotypicWeightedCC();

  void haplotypicTDT(map<int,int> &, int,bool);
  void haplotypicWeightedTDT();

  void haplotypicQTL(map<int,int> &, int,bool);

  int validN;     // Number of non-missing founders
  
  // Perform non-founder fill-in phasing?
  bool nonfounders;

};

#endif
