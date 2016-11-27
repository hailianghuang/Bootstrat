#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;

int main()
{

  const int MaxFam = 1000000;
  const int MaxSibSize = 10; // 11 sibs

// PARAMETERS...  

  // each family 50:50 [A]/[B]
  // each parent 50:50 [1]:[2], other parent is opposite

  // allele freqs
  double p_A;
  double p_B;

  // strata effects 
  double strat_A; 
  
  // basic kid VCs
  double Vpoly, Vc, Vn;
  double a, d;
  
  // parental values
  double VpolyP, VcP, VnP; 
  double aP, dP;
  
  int N, SibSize;
  int pat1, pat2, mat1, mat2;
  double g22, g12, g11;
  double mg22, mg12, mg11;
  double pat_pheno, mat_pheno;

  int g[MaxSibSize][3];
  int g2[MaxSibSize][3];
  double gv[MaxSibSize];

  double resC, resPoly[MaxSibSize], resN[MaxSibSize];
  double pheno[MaxSibSize];

  int i,j;
  double u1, u2, z1, z2;

  int mode_within, mode_kid;

  // RANDOMIZE
  srand(time(0));
  
  // allele freq set, strata
  cin >> p_A >> p_B >> strat_A;
  
  // more strata
  cin >> mode_within >> mode_kid;

  // kid parameters
  cin >> a >> d;
  cin >> Vpoly >> Vc >>  Vn;

  // parent parameters
  cin >> aP >> dP;
  cin >> VpolyP >> VcP >> VnP ; 

  // sample sizes
  cin >> N >> SibSize;

  // genotypic scores
  g22 =  a;
  g12 =  d;
  g11 = -a;


// loop for each sibshp

  for (j=1; j <= N; j++)
  {
      
      // Step 1. Generate family strata type
      
      int patAB = 1;
      if (rand()/double(RAND_MAX) < 0.5 ) patAB = 2;

      int matAB = 1;
      if (mode_within==0) matAB = patAB;
      else if (rand()/double(RAND_MAX) < 0.5 ) matAB = 2;

                        
      // Step 1. Generate Parental Mating Types
      
      // q > 1
      // p > 2 : increaser allele

      pat1 = 1;
      pat2 = 1;
      mat1 = 1;
      mat2 = 1;

      double ppat, pmat;

      if (patAB==1) ppat = p_A;
      else ppat = p_B;
      
      if (matAB==1) pmat = p_A;
      else pmat = p_B;

      if (rand()/double(RAND_MAX) < ppat) pat1 = 2;
      if (rand()/double(RAND_MAX) < ppat) pat2 = 2;

      if (rand()/double(RAND_MAX) < pmat) mat1 = 2;
      if (rand()/double(RAND_MAX) < pmat) mat2 = 2;


      // Step 2. Generate Offspring Genotypes

      for (i=1; i<=SibSize; i++)
	{
	  
	  if (rand()/double(RAND_MAX) < 0.5) 
	    {
	      g[i][1] = pat1;
	      g2[i][1] = 1;
	    }
	  else
	    {
	      g[i][1] = pat2;
	      g2[i][1] = 2;
	    }


	  if (rand()/double(RAND_MAX) < 0.5) 
	    {
	      g[i][2] = mat1;
	      g2[i][2] = 3;
	    }
	  else
	    {
	      g[i][2] = mat2;
	      g2[i][2] = 4;
	    }
    }


      // paternal genotypic score
      double gv_pat, gv_mat;
      if ( (pat1 == 2) && (pat2 == 2) ) gv_pat = g22;
      if ( (pat1 == 1) && (pat2 == 2) ) gv_pat = g12;
      if ( (pat1 == 2) && (pat2 == 1) ) gv_pat = g12;
      if ( (pat1 == 1) && (pat2 == 1) ) gv_pat = g11;
      
      // paternal genotypic score
      if ( (mat1 == 2) && (mat2 == 2) ) gv_mat = g22;
      if ( (mat1 == 1) && (mat2 == 2) ) gv_mat = g12;
      if ( (mat1 == 2) && (mat2 == 1) ) gv_mat = g12;
      if ( (mat1 == 1) && (mat2 == 1) ) gv_mat = g11;



      // Step 3. Generate Offspring genetic values conditional on genotypes
      
      for (i=1; i<=SibSize; i++)
      {
	  if ( (g[i][1] == 2) && (g[i][2] == 2) ) gv[i] = g22;
	  if ( (g[i][1] == 1) && (g[i][2] == 2) ) gv[i] = g12;
	  if ( (g[i][1] == 2) && (g[i][2] == 1) ) gv[i] = g12;
	  if ( (g[i][1] == 1) && (g[i][2] == 1) ) gv[i] = g11;
      }

      
      // Generate PATERNAL & MATERNAL POLYGENES
      u1 = rand() / double(RAND_MAX);
      u2 = rand() / double(RAND_MAX);
      z1 = sqrt(-2*log(u1)) * cos(2*M_PI*u2);
      z2 = sqrt(-2*log(u1)) * sin(2*M_PI*u2);
      

      double patA = z1;
      double matA = z2;

      double patPoly = patA * sqrt(VpolyP);
      double matPoly = matA * sqrt(VpolyP);

      // Generate offspring Polygenic influence
      for (i=1; i<=SibSize; i++)
      {
	  u1 = rand() / double(RAND_MAX);
	  u2 = rand() / double(RAND_MAX);
	  z1 = sqrt(-2*log(u1)) * sin(2*M_PI*u2);
	  
	  resPoly[i]  = (z1 * sqrt(Vpoly))/sqrt((double)2);
	  resPoly[i] += (((patA+matA)/2) * sqrt(Vpoly)) / sqrt((double)2);
	}

      
      // Generate FAMILY WIDE SHARED ENV..

      u1 = rand() / double(RAND_MAX);
      u2 = rand() / double(RAND_MAX);
      z1 = sqrt(-2*log(u1)) * cos(2*M_PI*u2);
      
      resC = z1 * sqrt(Vc);
 


// Step 5. Generate nonshared variable Vn
//         mean 0 variance Vn

      u1 = rand() / double(RAND_MAX);
      u2 = rand() / double(RAND_MAX);
      z1 = sqrt(-2*log(u1)) * cos(2*M_PI*u2);
      z2 = sqrt(-2*log(u1)) * sin(2*M_PI*u2);

      double patN = z1 * sqrt(Vn);
      double matN = z2 * sqrt(Vn);
      
      for (i=1; i<=SibSize; i++)
      {
	  u1 = rand() / double(RAND_MAX);
	  u2 = rand() / double(RAND_MAX);
	  z1 = sqrt(-2*log(u1)) * cos(2*M_PI*u2);
	  
	  resN[i] = z1 * sqrt(Vn);
      }



  
// Step 6. Calculate trait :  Pheno = Step 3 + Step 4 + Step 5
      
      pat_pheno = gv_pat + patPoly + resC + patN;
      mat_pheno = gv_mat + matPoly + resC + matN;

      // between-family stratification
      if (patAB==1) pat_pheno += strat_A;
      if (matAB==1) mat_pheno += strat_A;
      
      double rtmp = rand() / double(RAND_MAX);

      // what to do with the kids??? vis a vis stratificaiton effects??
      for (i=1; i<=SibSize; i++)
	 {
	     pheno[i] = gv[i] + resPoly[i] + resC + resN[i];
	     
	     // kid-stratification mode
	     // 1 = average of pat/mat
	     // 2 = random of pat/mat, same for each kid
	     // 3 = random of pat/mat, diff. for each kid
	     // 4 = largest of pat/mat
	     double kid_strat=0;
	     
	     // average
	     if (mode_kid==1) 
	     {
		 if (patAB==1)
		     kid_strat = 0.5;
		 if (matAB==1)
		     kid_strat += 0.5;
		 
		 pheno[i] += kid_strat * strat_A;
	     }
		     

	     // random parent, within family
	     if (mode_kid==2)
	     {
		 
		 if (rtmp > 0.5) kid_strat = patAB;
		 else kid_strat = matAB;
		 if (kid_strat==1) pheno[i] += strat_A;
	     }
	     

	     // random parent, random kid
	     if (mode_kid==3) 
	     {
		 if (rand()/double(RAND_MAX) > 0.5) kid_strat = patAB;
		 else kid_strat = matAB;
		 if (kid_strat==1) pheno[i] += strat_A;
	     }
	     
	     // largest
	     if (mode_kid==4) 
	     {
		 if ( patAB==1 || matAB==1 ) kid_strat = 1;
		 else kid_strat = 2;
		 if (kid_strat==1) pheno[i] += strat_A;
	     }
	     	     
	 }
      
      // Output
      
  // parents

  cout << j 
       << " 1 " // individual ID
       << " 0 0 " // founder
       << " 1 " // male
       << " " << pat1 << " " << pat2 // genotypes
       << " " << pat_pheno 
       << "\n";

  cout << j 
       << " 2 " // individual ID
       << " 0 0 " // founder
       << " 2 " // male
       << " " << mat1 << " " << mat2 // genotypes
       << " " << mat_pheno 
       << "\n";
  
  for (i=1; i<=SibSize; i++)
  {
      cout << j << " "
	   << i+2 // individual ID
	   << "  1 2 " // parents
	   << " 1 " // male
	   << " " << g[i][1] << " " << g[i][2] // genotypes
	  //	  << "  " << g2[i][1] << " " << g2[i][2] // genotypes
	   << "   " << pheno[i] // no quantitative trait ;
	   << "\n";
  }
  
    } // end J loop
  
  

  
}
