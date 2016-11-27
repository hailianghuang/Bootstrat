

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
#include <cmath>
#include "options.h"
#include "plink.h"

#define EPS  0.00001
#define EPS2 0.0000001

using namespace std;

void Plink::preCalcSinglePoint()
{
  m1.resize(0);
  m2.resize(0);
  pos.resize(0);
    
  // All single markers
  for (int i=par::run_start;i<=par::run_end;i++)
    {
      m1.push_back(i);
      m2.push_back(i);
      pos.push_back(0);
    }

  // Final analysis is genome-wide phenotype on IBD 
  m1.push_back(-1);
  m2.push_back(-1);
  pos.push_back(0);


}

void Plink::preCalcMultiPoint()
{
  
  ////////////////////////////////////////////////////////////////////
  // Mltipoint marker map to span range pos1-fringe to pos(nl)+fringe

  // Uniform grid in cM space; Uniform grid in marker space (i.e. X
  // positions between each marker pair).
  
  // Reset map

  m1.resize(0); m2.resize(0); pos.resize(0);
  
  // For each position on the cM map, determine the two flanking
  // markers and the proportional distance between the two.

  // Uniform map in cM space
  if (par::inter_grid==0)
    {
      for (double p=locus[par::run_start]->pos - par::fringe;
	   p <= locus[par::run_end]->pos + par::fringe;
	   p += par::grid)
	{
	  
	  // find marker that comes before this position...

	  if (p < locus[par::run_start]->pos)
	    { 
	      m1.push_back(-1);
	      m2.push_back(par::run_start);
	      double p2 = (p- (locus[par::run_start]->pos - par::fringe))
		/ par::fringe;
	      pos.push_back(p2); 
	    }
	  
	  // ... in the normal range...
	  for (int i=par::run_start; i<par::run_end; i++)
	    {	  
	      if (p >= locus[i]->pos && 
		  p < locus[i+1]->pos )
		{
		  double p2 = (p-locus[i]->pos)
		    / (locus[i+1]->pos - 
		       locus[i]->pos);
		  pos.push_back(p2);
		  m1.push_back(i);
		  m2.push_back(i+1);		  
		}
	    }
	  
	  // ... or after the last marker
	  if (p >= locus[par::run_end]->pos)
	    {
	      m1.push_back(par::run_end);
	      m2.push_back(-1);
	      double p2 = (p-locus[par::run_end]->pos)
		/ par::fringe;
	      pos.push_back(p2); 
	    }    
	}
    }
  else
    {
      // Uniform map in marker-space
      
      // ...before the first
      for (int j=0;j<par::inter_grid;j++)
	{
	  m1.push_back(-1);
	  m2.push_back(par::run_start);
	  pos.push_back((double)j/(double)par::inter_grid); 
	}
      
      // normal range
      for (int i=par::run_start; i<par::run_end; i++)
	{	  
	  for (int j=0;j<par::inter_grid;j++)
	    {
	      m1.push_back(i);
	      m2.push_back(i+1);
	      pos.push_back((double)j/(double)par::inter_grid); 
	    }
	}
      
      // after last marker
      for (int j=0;j<=par::inter_grid;j++)
	{
	  m1.push_back(par::run_end);
	  m2.push_back(-1);
	  pos.push_back((double)j/(double)par::inter_grid); 
	}
    }
  
  printLOG("Multipoint map has "+int2str(m1.size())+" positions\n");

  // Final analysis is genome-wide phenotype on IBD 
  m1.push_back(-1);
  m2.push_back(-1);
  pos.push_back(0);

  if (par::verbose)
    {
      cout << "Multipoint marker map information\n";
      for (int i=0; i<pos.size();i++)
	{
	  cout << i << "\t" 
	       << m1[i] << "\t" 
	       << m2[i] << "\t" 
	       << pos[i] << "\n";
	}
    }
 
}

vector<double> Plink::calcSinglePoint(vector<Z> & IBD, Z IBDg)
{
    vector<double> pihat;

    // Singlepoint individual loci
    for (int l=0; l<nl; l++)
    {	  
      
      //pihat.push_back( (IBD[l].z1*0.5+IBD[l].z2) - (IBDg.z1*0.5 + IBDg.z2)  );      
      pihat.push_back( (IBD[l].z1*0.5+IBD[l].z2) );      
      
      if (par::multi_output) 
	{
 	  cout << "S " 
	       << pairid << " " 
 	       << locus[m1[l]]->name << "\t"
	       << IBD[l].z0 << " "
	       << IBD[l].z1 << " "
	       << IBD[l].z2 << " => "
	       << IBD[l].z1*0.5+IBD[l].z2 << " "
	       << (IBDg.z1*0.5 + IBDg.z2) << "\n";	  
	}
    }

    // Final analysis is genome-wide IBD
    if (!par::done_global_pihat) pihat_G.push_back( (IBDg.z1*0.5 + IBDg.z2) );
    
    return pihat;
}


vector<double> Plink::calcMultiPoint(vector<Z> & IBD, Z IBDg, ofstream & MP)
{



  // Hidden Markov Model to estimate IBD states given IBS states,
  // alleles frequencies and a genetic map

  // Q = { IBD0, IBD1, IBD2 }           hidden states
  // S = { IBS0, IBS1, IBS2 }           observed alphabet

  // A = transition probabilities matrix  ( genetic map, genome-wide IBD)
  // B = P(IBS|IBD)                     emission/output probabilities 

  // Question: find the generating path, i.e. sequence of hidden IBD states

  // Use forward/backward algorithm for posterior decoding problem
  
  bool is0zero = IBDg.z0 < EPS ? true : false;
  bool is1zero = IBDg.z1 < EPS ? true : false;
  bool is2zero = IBDg.z2 < EPS ? true : false;
  
  ///////////////////////////
  // Z0, Z1 and Z2
  
  double z0 = IBDg.z0;
  double z1 = IBDg.z1;
  double z2 = IBDg.z2;


  /////////////////////////////////////////
  // Calculate two IBD probabilities

  // For haploid genome pair
  
  // 00 ( 1 - (m-1)*G ) + (m-1)*G*(1-(1/(2^(m-1)-1)) ) 
  // 01 1 - 00
  // 10 (m-1)*G
  // 11 1 - (m-1)*G

  // P(haploid genome is IBD) = x1 and x2 

  // x1 * x2 = z2
  // (1-x1) * (1-x2) = z0
  // x^2 - x + 0.25 = 0 

  // Solution to this quadratic equation gives:
  // ax^2 + bx + c = 0
  // a = -1; b = 1-z0+z2; c = -z2

  double mA, mB;

  if (is2zero)
    {
      mA = 1 - log(z1) / log(2.0);
      mB=0;
    }
  else
    {
      double b = 1 - z0 + z2;

      // After 'nudging' IBD probabilies, this should always be
      // positive -- but allow for rounding errors with fabs
      double b2 = fabs( b*b - 4*z2 );
      double x = sqrt( b2 );
            
      double t1 = (-b + x ) / -2;
      double t2 = (-b - x ) / -2;

      mA = 1 - log(t1) / log(2.0);
      mB = 1 - log(t2) / log(2.0);

    }

  
  // L = 1R Q1 T1 Q2 T2 .. T(M-1) Q(M) 1C 
  
  // 1R = 1x3 vector
  // 1C = 3x1 vector
  
  // Return value: vector of pi-hats
  vector<double> pihat;
  
  // Working matrices
  vector<Z> left(nl);
  vector<Z> right(nl);
 
  // Q matrix diagonal elements
  double Q00, Q11, Q22;
  
  // Left conditional begins with first locus on diagonal

  left[0] = IBD[0];
  
  left[0].z0 /= IBDg.z0;
  left[0].z1 /= IBDg.z1;
  left[0].z2 /= IBDg.z2;
  
  if (is0zero) left[0].z0 = 0;
  if (is1zero) left[0].z1 = 0;
  if (is2zero) left[0].z2 = 0;
  
  
  // Scaling factor
  double S = 1.0/(left[0].z0 + left[0].z1 + left[0].z2);
  left[0].z0 *= S;
  left[0].z1 *= S;
  left[0].z2 *= S;
  
  // Right conditional initial point
  right[nl-1] = IBD[nl-1];
  
  right[nl-1].z0 /= IBDg.z0;
  right[nl-1].z1 /= IBDg.z1;
  right[nl-1].z2 /= IBDg.z2;
  
  if (is0zero) right[nl-1].z0 = 0;
  if (is1zero) right[nl-1].z1 = 0;
  if (is2zero) right[nl-1].z2 = 0;
  
  // Scaling factor
  S = 1.0 / ( right[nl-1].z0 + right[nl-1].z1 + right[nl-1].z2 );
  right[nl-1].z0 *= S;
  right[nl-1].z1 *= S;
  right[nl-1].z2 *= S;
  


  ///////////////////////
  // Left conditional 
  
  int l=1;

  for (int l2=par::run_start+1; l2<=par::run_end; l2++)
    {	  
      
      // Build transition matrix
      
      // 1 x 3 . 3x3 . 3x3 . ... 

      Z prev;

      prev = left[l-1];

      buildT(locus[l2]->pos - locus[l2-1]->pos,
	     is2zero,
	     mA,mB);
	     
             
      // Q diagonal elemnts in IBD[l]
      // Left conditional [ R(n+1) * T * R(n) ] 
      left[l].z0 = ( ( prev.z0 * T00
		       + prev.z1 * T10
		       + prev.z2 * T20 ) * IBD[l].z0 );
      
      left[l].z1 = ( ( prev.z0 * T01
		       + prev.z1 * T11
		       + prev.z2 * T21 ) * IBD[l].z1 );
      
      left[l].z2 = ( ( prev.z0 * T02
		       + prev.z1 * T12
		       + prev.z2 * T22 ) * IBD[l].z2 );
      
      left[l].z0 /= IBDg.z0;
      left[l].z1 /= IBDg.z1;
      left[l].z2 /= IBDg.z2;
      
      if (is0zero) left[l].z0 = 0;
      if (is1zero) left[l].z1 = 0;
      if (is2zero) left[l].z2 = 0;
      
      // Scaling factor (sum to 1)
      double S = left[l].z0 + left[l].z1 + left[l].z2;
      
      // TODO: remove this check
      if (S < EPS) left[l] = IBDg;
      else
	{
	  left[l].z0 /= S;
	  left[l].z1 /= S;
	  left[l].z2 /= S;
	}
      
      // Shift left
      l++;
      
    } // Next marker interval
  
  if (par::verbose)
    {
      cout << "LEFT CONDITIONAL\n";
      for (int i=0;i<left.size();i++)
	{
	  cout << left[i].z0 << "  " 
	       << left[i].z1 << "  " 
	       << left[i].z2 << "\n";
	}
    }
  
  
  
   ///////////////////////
   // Right conditional 
  
   l = nl-1-1;
 
   for (int l2=par::run_end-1; l2>=par::run_start; l2--)
     {	  
      
       buildT(locus[l2+1]->pos - locus[l2]->pos,
	      is2zero,
	      mA,mB);

       // Right conditional [ R(n+1) * T * R(n) ] 
       right[l].z0 = ( ( right[l+1].z0 * T00
 			+ right[l+1].z1 * T10
 			+ right[l+1].z2 * T20 ) * IBD[l].z0 );
      
       right[l].z1 = ( ( right[l+1].z0 * T01
 			+ right[l+1].z1 * T11
 			+ right[l+1].z2 * T21 ) * IBD[l].z1 );
      
       right[l].z2 = ( ( right[l+1].z0 * T02
 			+ right[l+1].z1 * T12
 			+ right[l+1].z2 * T22 ) * IBD[l].z2 );
       
       // COndition on genome-wide IBD
       right[l].z0 /= IBDg.z0;
       right[l].z1 /= IBDg.z1;
       right[l].z2 /= IBDg.z2;
       
       if (is0zero) right[l].z0 = 0;
       if (is1zero) right[l].z1 = 0;
       if (is2zero) right[l].z2 = 0;
       
       // Scaling factor (sum to 1)
       double S = right[l].z0 + right[l].z1 + right[l].z2;
       if (S < EPS) right[l] = IBDg;
       else
	 {
	   right[l].z0 /= S;
	   right[l].z1 /= S;
	   right[l].z2 /= S;
	 }

       // Shift right
       l--;
      
     } // Next marker interval
  
  
   if (par::verbose)
     {
       cout << "RIGHT CONDITIONAL\n";
       for (int i=0;i<right.size();i++)
 	cout << right[i].z0 << "  " 
 	     << right[i].z1 << "  " 
 	     << right[i].z2 << "\n";
       cout <<"\n";
     }
  
    
   ///////////////////////////
   // Multipoint calculation
   
   if (par::verbose) cout << "FULL CONDITIONAL\n";
	 
   // skip last position (gIBD)
   for (int i=0; i<pos.size()-1; i++) 
     {
      
       double p1, p2;
     
       if (m1[i]==-1) 
 	p1 = locus[par::run_start]->pos - par::fringe;
       else p1 = locus[m1[i]]->pos;
      
       if (m2[i]==-1) 
 	p2 = locus[par::run_end]->pos + par::fringe;
       else p2 = locus[m2[i]]->pos;
      
       double d1 = pos[i] * (p2-p1);
       double d2 = (1 - pos[i]) * (p2-p1);
      
       // P0 = L * TL * Q0 * TR * R;
       // P1 = L * TL * Q1 * TR * R;
       // P2 = L * TL * Q2 * TR * R;

       //  L * TL

       // Left T matrix
       buildT(d1,is2zero,mA,mB);
       
       // Left & right conditional
       Z L;
       
       if (m1[i]==-1) { 
 	 L.z0 = L.z1 = L.z2 = 1;
       }
       else L = left[m1[i] - par::run_start];
       
       // * Q{0/1/2} [ Q 3x3 matrix -- just extracts elements ]
       
       Z M0;
       Z M1;
       Z M2;
       
       M0.z0 = ( L.z0 * T00
		 + L.z1 * T10
		 + L.z2 * T20 ) ;
       
       M1.z1 = ( L.z0 * T01
		 + L.z1 * T11
		 + L.z2 * T21 ) ;
       
       M2.z2 = ( L.z0 * T02
		 + L.z1 * T12
		 + L.z2 * T22 ) ;
       
       
       // * TR (transpose)
       buildT(d2,is2zero,mA,mB);
       
       // Finally, 1x3 . 3x1 = 1x1
       
       Z R;      
       
       if (m2[i]==-1) { 
	 R.z0 = R.z1 = R.z2 = 1;
       }
       else R = right[m2[i] - par::run_start];
       
       double P0 = M0.z0 * (( T00*R.z0 + T10*R.z1 + T20*R.z2)); 
       double P1 = M1.z1 * (( T01*R.z0 + T11*R.z1 + T21*R.z2)); 
       double P2 = M2.z2 * (( T02*R.z0 + T12*R.z1 + T22*R.z2)); 

       // Standardized P
       double S = 1.0/(P0+P1+P2);
       P0 *= S;
       P1 *= S;
       P2 *= S;
       
       // Apply Bayes Theorem 
       S = 1.0 / (P0*IBDg.z0 + P1*IBDg.z1 + P2*IBDg.z2 );
       P0 = (P0*IBDg.z0) * S;
       P1 = (P1*IBDg.z1) * S;
       P2 = (P2*IBDg.z2) * S;
       
       if (par::verbose)
	 cout << P0 << "\t" << P1 << "\t" << P2 << "\n";

       if (par::multi_output)
	 {
	   double p1, p2;
	   string n1, n2;
	  
 	  if (m1[i]==-1) 
 	    {
 	      p1 = locus[par::run_start]->pos - par::fringe;
 	      n1 = "fringe";
 	    }
 	  else 
 	    {
 	      p1 = locus[m1[i]]->pos;
 	      n1 = locus[m1[i]]->name;
 	    }
	  
 	  if (m2[i]==-1) 
 	    {
 	      p2 = locus[par::run_end]->pos + par::fringe;
 	      n2 = "fringe";
 	    }
 	  else 
 	    {
 	      p2 = locus[m2[i]]->pos;
 	      n2 = locus[m2[i]]->name;
 	    }
	  
 	  MP << par::run_chr << " "
	     << pairid << " " 
	     << p1 + pos[i] * (p2-p1) << " " 
	     << P0 << "  " 
	     << P1 << "  " 
	     << P2 << "  "
	     << (P1*0.5+P2) << "  " 
	     << (IBDg.z1*0.5 + IBDg.z2) << "\n";
	 }
       
       
       // Record pi-hat estimate       
       pihat.push_back( (P1*0.5+P2) );	  
       
     }
  
   
   // Final analysis is genome-wide IBD
   if (!par::done_global_pihat) 
     pihat_G.push_back( (IBDg.z1*0.5 + IBDg.z2) );
   
   return pihat;

}



void Plink::buildT(double G, bool z2zero, double mA, double mB)
{  
  

  /////////////////////////////////////////
  // Build 2x2 haploid transition matrices

  //                  m
  //                 ------------
  //                 0        1
  //                 ------------
  //
  //  m-1 |  0  |    X        t
  //      |  1  |    1-X      1-t

  // 
  // where t = mA * G = recombination fraction theta
  //      mA = number of meioses separating the haploid genomes
  //       G = genetic distance in Morgans
  // 
  
  double A01 = (1- pow(1-G,mA-2) * (G*G+(1-G)*(1-G)))/(pow(2,mA-1)-1);
  double A00 = 1 - A01;
  double A11 = pow((1-G),(mA-2)) * (G*G+(1-G)*(1-G) );
  double A10 = 1 - A11;
  
  double B00;
  double B01;
  double B10;
  double B11;

  
  if (z2zero) 
    {
      B00 = B11 = 1;
      B01 = B10 = 0;
    }
  else
    {
      B01 = (1- pow(1-G,mB-2) * (G*G+(1-G)*(1-G)))/(pow(2,mB-1)-1);
      B00 = 1 - B01;
      B11 = pow((1-G),(mB-2)) * (G*G+(1-G)*(1-G) );
      B10 = 1 - B11;
    }
 

//   if (par::verbose)
//     {
//       cout << "A[i,j] = \n";
//       cout << "\t" << A00 << "\t" << A01 << "\n"
// 	   << "\t" << A10 << "\t" << A11 << "\n\n";
      
//       cout << "B[i,j] = \n";
//       cout << "\t" << B00 << "\t" << B01 << "\n"
// 	   << "\t" << B10 << "\t" << B11 << "\n\n";
//     }
  
  ///////////////////////////
  // Build transition matrix
  
  T00 = A00*B00;              
  T01 = A00*B01 + A01*B00;       
  T02 = A01*B01;
  
  T10 = A00*B10 + A10*B00; 
  T11 = A00*B11 + A01*B10 + A10*B01 + A11*B00;
  T12 = A01*B11 + A11*B01;
  
  T10 /= 2;
  T11 /= 2;
  T12 /= 2;

  T20 = A10*B10; 
  T21 = A10*B11 + A11*B10;
  T22 = A11*B11;
  

//   if (par::verbose)
//     {
//       cout << "cM = " << G <<  "\n";
//       cout << "transition matrix\n"
// 	   << T00 << "\t" << T01 << "\t" << T02 << "\n"
// 	   << T10 << "\t" << T11 << "\t" << T12 << "\n"
// 	   << T20 << "\t" << T21 << "\t" << T22 << "\n";
//     }
  
}

