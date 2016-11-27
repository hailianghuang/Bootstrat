
#include <iostream>
#include <cmath>
#include "crandom.h"

using namespace std;


int main()
{
  
  // number of individuals
  int n;

  // number of alleles
  int na;

  // residual variance
  double res;

  int seed;

  cin >> seed;

  cin >> n // sample size
      >> na; // number alleles;
  
  vector<string> name;
  vector<int> dose;
  vector<double> f;//freqs
  vector<double> a;//effect size;

  for (int i=0; i<na; i++)
    {
      string s;
      cin >> s;
      name.push_back(s);
    }

  for (int i=0; i<na; i++)
    {
      int d;
      cin >> d;
      dose.push_back(d);
    }

  for (int i=0; i<na; i++)
    {
      double fq;
      cin >> fq;
      f.push_back(fq);
    }

  for (int i=0; i<na; i++)
    {
      double t;
      cin >> t;
      a.push_back(t);
    }
  
  
  // residual variance

  cin >> res;
  
  CRandom::srand(seed);

  //# genotypes 
  
  for (int i=0; i<n; i++)
    {
      double r1 = CRandom::rand();
      double r2 = CRandom::rand();      
      int g1 = 0;
      int g2 = 0;
      
      double c = f[0];
      for (int j=1; j<na; j++)
	{
	  if ( r1 > c )
	    g1 = j;
	  if ( r2 > c )
	    g2 = j;
	  c += f[j];
	}
      
      // Residual variance
      double u1 = CRandom::rand();
      double u2 = CRandom::rand();
      double z = sqrt(-2*log(u1)) * cos(2*M_PI*u2);

      // Gene effect
      z += a[g1] + a[g2];
      
      cout << i+1 << " 1 var1 ";
      
      cout << name[g1] << " " 
	   << dose[g1] << " ";

      cout << name[g2] << " " 
	   << dose[g2] << " ";
	
      cout << "\n";
      
      cout << "PHENO ";
      cout << i+1 << " 1 0 0 1 " << z << "\n";
      
    }
  
  
  exit(0);
}
