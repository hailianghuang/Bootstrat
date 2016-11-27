#ifndef __RTEST_H__
#define __RTEST_H__

#include <iostream>

enum rTestType { RTEST_VANILLA = 0 , 
		 RTEST_FREQWEIGHT = 1 ,
		 RTEST_VARIABLE_THRESHOLD = 2 };
		 
enum rTestHyp { RTEST_1SIDED_HIGH = 1 , 
		RTEST_1SIDED_LOW  = 2 , 
		RTEST_2SIDED      = 3 };

class RVData {
 public: 

  int n; // variant #
  int l; // ind #
  int c; // allele count (A1)
  int tc; // total allele count 

  // Diploid, individual data
  RVData(int n, int l, int c ) : n(n), l(l), c(c) { tc = 2; }

  // Pooled level data
  RVData(int n, int l, int c, int tc ) : n(n), l(l), c(c) , tc(tc) { }
  
  friend std::ostream & operator<<(std::ostream & out, RVData & g)
  { 
    out << g.n << " " << g.l << " " << g.c;
    return out;
  }
  
  bool operator< (const RVData & b) const
    {
      if ( n < b.n ) return true;
      if ( n > b.n ) return false;
      return l < b.l;
    }


};

#endif
