

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

extern ofstream LOG;

void Plink::readHaplotypeFormat()
{

  
  //////////////////////
  // Check files exist
  
  checkFileExists(par::hapfilename);
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
  
  
  ///////////////////////////////////////////////
  // HAP file map mode
  
  // 0 = default = 0 1  no missing data
  // 1 =           A/C/G/T  and 0 = missing

  bool mode01 = par::hapfile_coding == 0;

  //////////////////////////////////////////////
  // First read a reference file? 

  map<string,string> refallele;
  map<string,string> refallele2;

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

  if ( ! par::SNP_major)
    {
      for (int i=0; i<sample.size(); i++)
	{
	  sample[i]->one.resize(nl_actual,code);
	  sample[i]->two.resize(nl_actual,false);
	}
    }
  else
    {
      for (int l=0; l<SNP.size(); l++)
	{
	  SNP[l]->one.resize(sample.size(),code);
	  SNP[l]->two.resize(sample.size(),false);
	}
    }

  
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
  // .haps file
  
  FILE * PED;
  
  PED = fopen64(par::hapfilename.c_str(),"r");

  if ( PED == NULL )
    error("Problem opening HAPS file, errno = "+int2str(errno));
  
  printLOG("Reading haplotypes from [ " + par::hapfilename + " ]\n");

  // Assume format is ID 0101010010101010
  //                  ID 1010110101010100
  
  // if (mode01 == F) then we have alleles, plus possible missing data
  // otherwise, assume 01 are two alleles, no missing data

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


  nl_all = locus.size();
    
  
  // Do we need to supply allele names?
  if ( mode01 && ! par::ref_file ) 
    {
      for (int l=0; l<nl_all; l++)
	{
	  locus[l]->allele1 = "A";
	  locus[l]->allele2 = "B";
	}
    }


  // Whether or not we want to look at a locus is in the include[] vector
  // The genomic position of locus i is k=include_pos[i] -> locus[k]
  
  bool fatal = false;
  string fmsg = "";

  while( ! feof(PED) )
    {
      
      string fid1 = "";
      string iid1 = "";
      string hap1 = "";
      
      string fid2 = "";
      string iid2 = "";
      string hap2 = "";

      int f = 0;
      
      if ( readString( PED , fid1 ) ) f++;
      if ( fid1 == "" ) continue;
      if ( readString( PED , iid1 ) ) f++;
      if ( readString( PED , hap1 ) ) f++;
 
      if ( readString( PED , fid2 ) ) f++;
      if ( fid1 == "" ) continue;
      if ( readString( PED , iid2 ) ) f++;
      if ( readString( PED , hap2 ) ) f++;
      
      // Which person, if any?

      if ( fid1 != fid2 ) continue;
      if ( iid1 != iid2 ) continue;
      
      // Do we have the correct amount of data?
      
      if ( hap1.size() != nl_all || 
	   hap2.size() != nl_all )
	continue;

      // Otherwise, find the person, and map
      // the genotypes

      map<string,int>::iterator peri 
	= iperson.find( fid1+"_"+iid1 );      
      
      Individual * person = 
	peri != iperson.end() ?
	sample[ peri->second ] : NULL ; 
      
      int ip = peri->second;
      
      if ( person == NULL ) 
	continue;      

      for (int l = 0; l < nl_all; l++)
	{
	  
	  char one = hap1[l];
	  char two = hap2[l];
	  
	  // Assume mode01 for now, and SNP-major reading mode
	  
	  if ( mode01 )
	    {
	      Locus * loc = locus[l];
	      
	      // Assume if (par::SNP_major)
	      
	      // 00 hom
	      
	      if ( one=='0' ) 
		{
		  if ( two=='0' )
		    {
		      SNP[l]->one[ip] = false;
		      SNP[l]->two[ip] = false;
		    }	      
		  else // 01 het
		    {
		      SNP[l]->one[ip] = false;
		      SNP[l]->two[ip] = true;
		    }
		}
	      else
		{
		  if ( two == '1' )
		    {
		      SNP[l]->one[ip] = true;
		      SNP[l]->two[ip] = true;
		    }
		  else  // 10 het als, under mode01
		    {
		      SNP[l]->one[ip] = false;
		      SNP[l]->two[ip] = true;
		    }		  
		}
	    }
	  
	  // Next genotype in this pair of haplotypes
	}
      
      // Next pair of haplotypes
    }
  
    
    fclose(PED);
  
}


