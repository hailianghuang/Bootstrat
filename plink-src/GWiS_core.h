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


#ifndef __GWiS_CORE_H__
#define __GWiS_CORE_H__

//define a number that is essentially '0'
#define EPS 1E-7
//define the R2 corresponding to the variance inflation factor of 5
#define VIF_R2 0.2
//the maximum model size
#define MAX_K 20

using namespace std;

class SNP_SUMMARY{

 public:
  SNP_SUMMARY();
  
  double var;
  double mean;
  double pheno_cov; 
  
  int gene_id;
};

class GWiS_math{
  
  //get the lower element of the LD matrix
  static double & getLDlower(vector<vector<double> > &, int, int );

 public:
  static bool UpdateCovP(int, vector <SNP_SUMMARY>&, vector<vector<double> > &);
  static bool UpdateCovMat(int,  vector<vector<double> > &);
  static double getSSM(SNP_SUMMARY &, double);
  static int getBestSNP(vector <SNP_SUMMARY> &, vector<vector<double> > &, vector <int> &);
  static double SSM2BIC(double, double, int, int, double);
  static double findBestModel(double, int, double, int, vector<int> &, vector <SNP_SUMMARY> &, vector<vector<double> > &);
  static double calT(int , vector<vector<double> > & );


};

#endif
