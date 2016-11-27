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


#ifndef __GWiS_H__
#define __GWiS_H__

using namespace std;

class GWiS{
  
  double getGenoVar(CSNP *, double, vector<Individual*>, int*);
  double getGenoMean(CSNP *, vector<Individual*>,int* );

  double getPhenoVar(vector<Individual*>, double, int*);
  double getPhenoMean(vector<Individual*>, int*);

  double getSNPCov(CSNP *, CSNP *, double, double, vector<Individual*>, int*);
  double getSNPPhenoCov(CSNP *, vector<Individual*>, double, double, int*);
  
 public:
  GWiS(Plink &);


};

#endif
