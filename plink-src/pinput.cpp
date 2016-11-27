

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
#include <map>
#include <algorithm>
#include <bitset>
#include <limits>
#include <errno.h>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "sets.h"

extern ofstream LOG;

void Plink::readDataPooledFormat()
{

  
  //////////////////////
  // Check files exist
  
  checkFileExists(par::poolfile);
  checkFileExists(par::mapfile);
  checkFileExists(par::famfile);


  ///////////////////////////////////////////////
  // .map file
  
  vector<bool> include;
  vector<int> include_pos(0);
  int nl_actual=0;
  
  // Read in MAP file: this function also allocates 
  // slots for SNP-major mode

  readMapFile(par::mapfile,
	      include,
	      include_pos,
	      nl_actual);
  
  

  //////////////////////////////////////////////
  // First read a reference file? 

  map<string,string> refallele;
  map<string,string> refallele2;
  
  if ( ! par::ref_file )
    error("You need to supply a reference allele file, --reference");

  if ( par::ref_file )
    {
      set<string> mset;
      for (int l=0; l< locus.size(); l++)
	mset.insert( locus[l]->name );

      int notfound = 0;
      checkFileExists( par::ref_file_name );
      ifstream REF( par::ref_file_name.c_str(), ios::in );
      while ( ! REF.eof() )
	{
	  vector<string> tok = tokenizeLine( REF );
	  if ( tok.size() == 0 ) 
	    continue;
	  
	  if ( tok.size() != 2 && tok.size() != 3 )
	    error("Problem with line in [ " 
		  + par::ref_file_name 
		  + " ] :\n " + displayLine( tok ) );


	  if ( mset.find( tok[0] ) == mset.end() )
	    {
	      ++notfound;
	      continue;
	    }
	  
	  refallele.insert( make_pair( tok[0], tok[1] ));

	  // A second allele also specified?
	  if ( tok.size() == 3 ) 
	    refallele2.insert( make_pair( tok[0], tok[2] ));
	  
	}

      printLOG("Read reference alleles for " 
	       + int2str( refallele.size() ) 
	       + " sites\n");

      if ( notfound>0 )
	printLOG(int2str( notfound ) 
		 + " SNPs in reference but not in map file\n");
      if ( refallele.size() < locus.size() ) 
	printLOG(int2str( locus.size() - refallele.size() ) 
		 + " SNPs in map file but not in reference\n");
      REF.close();
    }


  ///////////////////////////////////////////////
  // .fam
  
  readFamFile(par::famfile);
  

  // Allocate space for individual-major mode, set to 
  // missing by default...

  // Either missing (TF) by default; or reference allele (FF)
  bool code = ! par::ref_file;

  if ( par::ref_file )
    {
      for (int l=0; l< locus.size(); l++)
	{	  
	  map<string,string>::iterator i = refallele.find( locus[l]->name );
	  map<string,string>::iterator i2 = refallele2.find( locus[l]->name );
	  
	  // If we cannot find, we need to set genotypes to missing instead
	  
	  if ( i != refallele.end() )
	    {
	      locus[l]->allele1 = i->second;
	      
	      if ( i2 != refallele2.end() )
		locus[l]->allele2 = i2->second;
	      else if ( par::lfile_allele_count )
		locus[l]->allele2 = i->second + "v";
	    }
	  else
	    {
	      for (int i=0; i<sample.size(); i++)
		{
		  if ( par::SNP_major )
		    {
		      SNP[l]->one[i] = true;
		      SNP[l]->two[i] = false;
		    }
		  else
		    {
		      sample[i]->one[l] = true;
		      sample[i]->two[l] = false;
		    }
		}
	    }
	}
    }


  ///////////////////////////////////////////////
  // Any sets to read?
  
  if ( par::read_set ) 
    {
      cout << "abother111..\n";
      //readSet();
      if ( pS == NULL ) cout << "IS NULL!\n";
      cout << "hmm..\n";
      pS->initialiseSetMapping();
      cout << "abother..\n";
      rdata.resize(setname.size());
    }
  else
    rdata.resize(1);

  cout << "step1\n";

  ///////////////////////////////////////////////
  // .lgen
  
  FILE * PED;
  
  PED = fopen64(par::poolfile.c_str(),"r");
  if ( PED == NULL )
    error("Problem opening PGEN file, errno = "+int2str(errno));
    
  // We can now read any number of individual/genotype lines, in any
  // order; we also do not assume that all genotypes are given --
  // these will be missing by default


  map<string,int> imap;
  map<string,int> iperson;

  for (int i=0; i<include.size(); i++)
    {
      if ( include[i] ) 
	{
	  int k = include_pos[i];
	  imap.insert( make_pair ( locus[k]->name , k ) );
	}
    }
  
  for (int i=0; i<sample.size(); i++)
    {
      iperson.insert( make_pair 
		      ( sample[i]->fid + "_" + sample[i]->iid , 
			i ) );
    }
  

  // Format is 
  // PID (SET) VAR A1CNT TOT
  
  bool fatal = false;
  string fmsg = "";
  
  while( ! feof(PED) )
    {
      
      string pid = "";
      string snp = "";
      string one = "";
      string two = "";

      int f = 0;
      
      if ( readString( PED , pid ) ) f++;
      
      if ( pid == "" ) 
	continue;
      
      if ( readString( PED , snp ) ) f++;
      if ( readString( PED , one ) ) f++;
      if ( readString( PED , two ) ) f++;
      
      map<string,int>::iterator im = imap.find(snp);
	
      int k = 
	im != imap.end() ?
	im->second : -1;
      
      // Need to read second allele?
      
      int a1,a2;

      bool prob = false;
      if ( ! from_string<int>( a1, one , std::dec ) )
	prob = true;
      
      if ( ! from_string<int>( a2, two , std::dec ) )
	prob = true;
      
      if ( a1 < 0 || a2 < 0 )
	prob = true;
     
		
      map<string,int>::iterator peri 
	= iperson.find( pid+"_"+pid );      
      
      Individual * person = 
	peri != iperson.end() ?
	sample[peri->second] : NULL ; 
      
      
      // Ignore this genotype?
      
      if ( ( ! person ) || k < 0 ) 
	continue;
      
      vector<int> is;
      
      if ( par::read_set )
	{
	  for (int i=0;i<setname.size();i++)
	    {
	      map<int, set<int> >::iterator si = pS->setMapping.find(k);
	      set<int>::iterator si2 = si->second.find(i);
	      if ( si2 != si->second.end() )
		is.push_back(i);
	    }
	}
      else
	is.push_back(0);


      cout << "step2\n";
      int ip = peri->second;	
      Locus * loc = locus[k];
      cout << "ip , k " << ip << " " << k << "\n";

      // Add genotype count to pool store
	for (int j=0; j<is.size(); j++)
	  {
	    cout << " is[j] == " << j << " " << is[j] << "\n";
	    // person, locus, allele count, total # of alleles
	    rdata[is[j]].push_back( RVData(ip,k,a1,a1+a2) );
	  }

	// Give an error message
	if ( prob )
	  error("Problem with PGEN file format");

      }
    

    
    fclose(PED);
  
}


