

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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include "plink.h"
#include "options.h"
#include "helper.h"

//////////////////////////////////////
// Get allele frequencies for sample
// 

void Plink::filterSNPs_individual()
{

  printLOG("Applying filters (individual-major mode)\n");

  // vector to record which SNPs to be deleted
  vector<bool> del(locus.size(),false);

  // vector to record which individuals to be deleted
  vector<bool> indel(sample.size(),false);



  /////////////////////////////////////////////////
  // First display number of founders/nonfounders 

  cnt_f=0;
  for (int i=0; i<n; i++)
     if ( sample[i]->founder ) cnt_f++;	    
  printLOG(int2str(cnt_f)+" founders and "+int2str(n-cnt_f)+
	   " non-founders found\n");

  if (cnt_f<n) par::has_nonfounders = true;

  /////////////////////////////////////////////////
  // Remove individuals with too many missing calls

  double total_genotyping = 0;

  if (par::MAX_IND_MISSING < 1)
    {

      int n_removed = 0;
      int n_orig = n;
     
      for (int i=0;i<sample.size();i++)
	{

	  // Sum missingness over all SNPs
	  int m=0;

	  for (int l=0; l<locus.size();l++)
	    if ( sample[i]->one[l] && (!sample[i]->two[l]) ) m++;

	  // Too much missingness?
	  if ( (double)m/(double)locus.size() > par::MAX_IND_MISSING )
	    {
	      indel[i] = true;
	      n_removed++;
	    }
	}
    
      ////////////////////////////////////////
      // Save list of any removed individuals
      
      if (n_removed>0)
	{
	  string f = par::output_file_name + ".irem";
	  printLOG("Writing list of removed individuals to [ " + f + " ]\n");
	  ofstream REM;
	  REM.open(f.c_str(), ifstream::out);
	  for (int i=0;i<sample.size();i++)
	    if (indel[i])
	      REM << sample[i]->fid << "\t" << sample[i]->iid << "\n";
	  REM.close();      
	  
         
	  // And now remove these individuals, so that 
	  // SNP-based statistics are calculated with 
	  // these samples already excluded
	  
	  n_removed = deleteIndividuals(indel);
	  
	}

  
      stringstream s2;
      s2 << n_removed << " of " << n_orig 
	 << " individuals removed for low genotyping ( MIND > " 
	 << par::MAX_IND_MISSING << " )\n";
      printLOG(s2.str());
         

    } // end remove individuals

 
  /////////////////////////////////
  // Calculate or read from file? 

  if (par::af_read)
    {

      checkFileExists(par::af_file);
      printLOG( "Reading allele frequencies from [ " + par::af_file + " ] \n");
      
      // Make hash of original SNP names
      map<string,int> mlocus;
      for (int l=0;l<nl_all;l++)
	mlocus.insert(make_pair(locus[l]->name,l));
      map<string,int>::iterator ilocus;

      ifstream FRQ;
      FRQ.open(par::af_file.c_str());
      FRQ.clear();
      
      string dum1, dum2, dum3, dum4, dum5, dum6;
      string snpname;
      double freq;
      int nm;

      for (int l=0; l<nl_all; l++)
	{
	  locus[l]->freq = -1;
	  locus[l]->nm = 0;
	}

      // Skip header line
      FRQ >> dum1 >> dum2 >> dum3 >> dum4 >> dum5 >> dum6;
      
      while(!FRQ.eof())
	{
	  char cline[256];
	  FRQ.getline(cline,256,'\n');
	  string sline = cline;
	  if (sline=="") continue;
	  stringstream ss(sline); 
	  string buf; 
	  vector<string> tokens; 
	  while (ss >> buf)
	    tokens.push_back(buf);
	  if (tokens.size() == 0)
	    continue;
	  else if (tokens.size() != 6)
	    error("Problem with allele frequency line: 6 fields required:\n"+sline+"\n");

	  ilocus = mlocus.find(tokens[1]);
	  if (ilocus != mlocus.end())
	    {
	      locus[ilocus->second]->freq = atof(tokens[4].c_str());
	      locus[ilocus->second]->nm = atoi(tokens[5].c_str());
	    }
	}
      FRQ.clear();
      FRQ.close();

    }
  


  /////////////////////////////////
  // Consider each locus at a time
  
  vector<string> hetlist(0);

  int exc_miss=0;
  int exc_maf=0;

  vector<Locus*> no_founders_found_list; 

  for (int l=0; l<nl_all; l++)
    {

      if (!par::af_read)
	{
	  locus[l]->freq = 0;
	  // count 1 per allele, for frequency
	  locus[l]->nm = 0; 
	}

      // count 1 per genotype, for missingness
      int geno_nm = 0; 

      bool X = false;
      bool haploid = false;
      
      if (par::chr_sex[locus[l]->chr]) X=true;
      else if (par::chr_haploid[locus[l]->chr]) haploid=true;
      
      ///////////////////////////////
      // Iterate over each individual

      for (int i=0; i<n; i++)
	{
	  
	  Individual * person = sample[i];	    
	  
	  // Check female Y genotypes
	  if ( par::chr_Y[locus[l]->chr] && ! person->sex )
	    {
	      // Set to missing, unless in a RECODE mode
	      if ( ! par::preserve_all_genotypes )
		{
		  person->one[l] = true;
		  person->two[l] = false;
		}		  
	    }
	  
	  // For haploid heterozygosity check, also consider all individuals
 	  if ( haploid || ( X && person->sex==1 ) )
	    {
	      if ( (!person->one[l]) && person->two[l] )
	 	{
		  hetlist.push_back( person->fid + "\t" 
				     + person->iid + "\t" 
				     + locus[l]->name );
		  
		  // Set to missing, unless in RECODE mode
		  if ( ! par::preserve_all_genotypes )
		    {
		      person->one[l] = true;
		      person->two[l] = false;
		    }
		}
	    } 


	  // For missing genotypes, consider all individuals
	  if ( ! ( person->one[l] && (!person->two[l]) ) ) geno_nm++;

	  if (!par::af_read)	  
	    {
	      // For allele frequencies
	      // only consider founders	
	      if ( par::summ_nonfounders || person->founder ) 
		{
		  
		  if ( haploid || ( X && person->sex==1 ) )
		    {
		      
		      //////////////////
		      // Haploid counts
		      
		      // "1" allele count
		      
		      if ( (!person->one[l]) && (!person->two[l]) )   //  FF = hom(11)
			{
			  locus[l]->freq++;
			  locus[l]->nm++;
			}	
		      else if ( person->one[l] && person->two[l] )   //  TT = hom(22)
			{
			  locus[l]->nm++;
			}
		      
		    }
		  else
		    {
		      
		      //////////////////
		      // Autosomal count
		      
		      // "1" allele count
		      
		      if (!person->one[l])
			{ 
			  if (!person->two[l])  //   00 = hom(11)
			    {
			      locus[l]->freq+=2;
			      locus[l]->nm+=2;
			    }	
			  else                  //   01 = het(12)
			    {
			      locus[l]->freq+=1;
			      locus[l]->nm+=2;
			    }
			}
		      else if ( person->two[l] ) // 11 = hom(22)
			{
			  locus[l]->nm+=2;
			}
		    }
		  
		}
	    }
	}


      ///////////////////////////////////////	      
      // Calculate frequencies, if required
      
      if (!par::af_read)
	{
	  
	  if (par::af_count)
	    {
	      // Use freq to store count (keep as is)	
	      
	      // Use "bp" to store number of allele 2
	      locus[l]->bp = (long int)(locus[l]->nm - locus[l]->freq);
 
	      // Use "pos" to store number missing
	      locus[l]->pos = n - geno_nm;
	    }
	  else
	    {
	      if (locus[l]->nm>0)
		locus[l]->freq /= (double)locus[l]->nm;
	      else
		{
		  locus[l]->freq = 1;
		  // If we aren't getting rid of it anyway
		  if ( (double)geno_nm/(double)n >= (1-par::MAX_GENO_MISSING))
		    no_founders_found_list.push_back(locus[l]);
		}	  
	    }
	}
      
      /////////////////////////////////////////////////
      // Record total proportion of missingness

      double snp_genotyping = n>0 ? (double)geno_nm/(double)n : 0;
      total_genotyping += snp_genotyping;


      /////////////////////////////////////////////////
      // Exclude if SNP has too many missing genotypes

      if ( (double)geno_nm/(double)n < (1-par::MAX_GENO_MISSING))
	{
	  del[l] = true;
	  exc_miss++;
	}
      
      ////////////////////////////////////////////////
      // Make allele1 always the least common allele

      if ( (!par::af_count) && locus[l]->freq > 0.5 ) 
	{
	  // then we need to swap alleles

	  locus[l]->freq = 1 - locus[l]->freq;
	  
	  string tmp = locus[l]->allele2;
	  locus[l]->allele2 = locus[l]->allele1;
	  locus[l]->allele1 = tmp;
	  
	  for (int i=0; i<n; i++)
	    if ( sample[i]->one[l] == sample[i]->two[l] ) 
	      {
		sample[i]->one[l] = !sample[i]->one[l];
		sample[i]->two[l] = !sample[i]->two[l];
	      }
	  
	}

    }
  



  /////////////////////////////////////////////////
  // Save list of any heterozygous haploid alleles
  
  if (hetlist.size()>0)
    {
      printLOG(int2str( hetlist.size()) + " heterozygous haploid genotypes called\n");
      string f = par::output_file_name + ".hh";
      printLOG( "Writing list of heterozygous haploid genotypes to [ " + f + " ]\n");
      ofstream REM;
      REM.open(f.c_str(), ifstream::out);
      for (int i=0; i<hetlist.size(); i++)
	REM << hetlist[i] << "\n";
      REM.close();      
    }
  hetlist.clear();


  /////////////////////////////////////////////////
  // Save list of SNPs with no founders observed

  if (no_founders_found_list.size()>0)
    {
      printLOG(int2str( no_founders_found_list.size()) + " SNPs with no founder genotypes observed\n");
      printLOG("Warning, MAF will be set to 0 for these SNPs (see --nonfounders)\n");
      string f = par::output_file_name + ".nof";
      printLOG( "Writing list of these SNPs to [ " + f + " ]\n");
      ofstream NOF;
      NOF.open(f.c_str(), ifstream::out);
      for (int i=0; i<no_founders_found_list.size(); i++)
        NOF << no_founders_found_list[i]->name << "\n";
      NOF.close();
    }
  no_founders_found_list.clear();


  
  //////////////////////////
  // Write allele freq file

  if (par::af_write)
    {
      if (par::include_cluster_from_file)
	calcStratifiedAlleleFreqs();
      else
	{

	  ofstream FRQ;
	  string f = par::output_file_name + ".frq";
	  if (par::af_count) f += ".count";
	  
	  if (par::summ_nonfounders)
	    printLOG("Writing allele frequencies (all individuals) to [ " + f + " ] \n");
	  else
	    printLOG("Writing allele frequencies (founders-only) to [ " + f + " ]\n");
	  
	  if (par::af_count)
	    printLOG("Display counts rather than frequencies\n");

	  FRQ.open(f.c_str(), ifstream::out);
	  FRQ << setw(4) << "CHR" << " "
	      << setw(12) << "SNP" << " "
	      << setw(4) << "A1" << " "
	      << setw(4) << "A2" << " ";
	  if (par::af_count)
	    FRQ << setw(6) << "C1" << " "
		<< setw(6) << "C2" << " "
		<< setw(6) << "G0" << "\n";
	  else
	    FRQ << setw(12) << "MAF" << " " 
		<< setw(8) << "NCHROBS" 
		<< "\n";	  
	  
	  for (int l=0; l<nl_all; l++)
	    {
	      string a1 = locus[l]->allele1;
	      if (a1=="") a1="0";
	      FRQ << setw(4)  << locus[l]->chr  << " "
		  << setw(12) << locus[l]->name  << " "
		  << setw(4)  << a1  << " "
		  << setw(4)  << locus[l]->allele2  << " ";
	      
	      if (par::af_count)
		{
		  FRQ << setw(6) << int( locus[l]->freq ) << " "
		      << setw(6) << int( locus[l]->bp  ) << " "
		      << setw(6) << int( locus[l]->pos ) << "\n";
		}
	      else
		{
		  if ( locus[l]->freq >= 0 )
		    FRQ << setw(12) << locus[l]->freq << " ";
		  else
		    FRQ << setw(12) << "NA" << " ";
		  FRQ << setw(8) << locus[l]->nm << "\n";
		}
	    }
	  FRQ.close();
	  
	}

      shutdown();
    }



  /////////////////////////
  // Write HWE statistics
  
  if (par::HWD_test || par::HWD_report)
    {
      
      ofstream HWD;
      if (par::HWD_report)
	{

	  if (par::summ_nonfounders)
	    printLOG("Writing Hardy-Weinberg tests (all individuals) to [ " +
		     par::output_file_name + ".hwe ] \n");
	  else
	    printLOG("Writing Hardy-Weinberg tests (founders-only) to [ " +
		     par::output_file_name + ".hwe ] \n");
	 
	  string f = par::output_file_name + ".hwe";
	  HWD.open(f.c_str(), ifstream::out);
	  
	  HWD.precision(4);
      
	    HWD << setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(8) << "TEST" << " "
		<< setw(20) << "GENO" << " "
		<< setw(8) << "O(HET)" << " "
		<< setw(8) << "E(HET)" << " "
		<< setw(12) << "P_HWD" << " "
		<< "\n";
	}
	  
      int cnt=0, cnt_a=0, cnt_u=0;
      
      // Consider each locus
      
      for (int l=0; l<nl_all; l++)      
	{
	  
	  // Compute p-values for HWE test in cases, controls & all
	  
	  // Only consider founders
	  
	  int a11, a12, a22;
	  int u11, u12, u22;
	  int b11, b12, b22;
	  
	  a11=a12=a22=0;
	  u11=u12=u22=0;
	  b11=b12=b22=0;
	  
	  bool X = false;
	  bool haploid = false;
	  
	  if (par::chr_sex[locus[l]->chr]) X=true;
	  else if (par::chr_haploid[locus[l]->chr]) haploid=true;
	  
	  // Iterate over each individual, founders only
	  
	  for (int i=0; i<n; i++)
	    {
	      Individual * person = sample[i];	    
	      
	      ///////////////////////////////////////////////
	      // Only consider founders, & diploid genotypes
	      
	      if ( par::summ_nonfounders || person->founder )
		if ( ! ( haploid || ( X && person->sex==1 ) ) )
		  {	
		    
		    // Consider everybody, irrespective of phenotype
		    // (QT, C/C or missing)
		    
		    if (!person->one[l])
		      { 
			if (!person->two[l]) b11++;   //   00 = hom(11)
			else b12++;                   //   01 = het(12)
		      }
		    else if ( person->two[l] ) b22++; // 11 = hom(22)
		    
		    
		    if (par::bt)  // for binary trait, separately for cases/controls
		      {
			if (person->phenotype == 1)
			  {
			    
			    if (!person->one[l])
			      { 
				if (!person->two[l]) u11++;   //   00 = hom(11)
				else u12++;                   //   01 = het(12)
			      }
			    else if ( person->two[l] ) u22++; //   11 = hom(22)
			    
			  }
			else if (person->phenotype == 2)
			  {
			    if (!person->one[l])
			      { 
				if (!person->two[l]) a11++;   //   00 = hom(11)
			      else a12++;                     //   01 = het(12)
			      }
			    else if ( person->two[l] ) a22++; //   11 = hom(22)
			  }
			
		      }
		}
	     
	    } // next individual
	  

	  // Allele frequencies
	  double afreq = 0, ufreq = 0, freq = 0;
	  
	  bool include_cases = true;
	  bool include_controls = true;
	  if (par::qt)
	    freq = ( a11 + a12/2 ) / (double)( a11+a12+a22 );
	  else
	    {
	      afreq = ( a11 + (double)a12/2.0 ) / (double)( a11+a12+a22 );
	      ufreq = ( u11 + (double)u12/2.0 ) / (double)( u11+u12+u22 );
	      freq =  ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );

	      if ( a11+a12+a22 == 0 ) include_cases = false;
	      if ( u11+u12+u22 == 0 ) include_controls = false;		
	    }
	  
	  double exp_a11, exp_a12, exp_a22;
	  double exp_u11, exp_u12, exp_u22;

	  
	  
	  if (par::qt)
	    {
	      
	      double p;

	      if (par::HWD_standard)
		{
	      
		  double exp_11 = freq * freq * (a11+a12+a22);
		  double exp_12 = 2 * freq * (1-freq) * (a11+a12+a22);
		  double exp_22 = (1-freq) * (1-freq) * (a11+a12+a22);
		  
		  double chisq = ( (a11-exp_11)*(a11-exp_11) ) / exp_11 
		    + ( (a12-exp_12)*(a12-exp_12) ) / exp_12 
		    + ( (a22-exp_22)*(a22-exp_22) ) / exp_22 ;
		  
		  p = chiprobP(chisq,1);
		}
	      else
		p = SNPHWE( b12, b11, b22 );


	      if (par::HWD_report)
		{
		  HWD << setw(par::pp_maxsnp) << locus[l]->name << " "
		      << setw(8) << "ALL(QT)" << " "
		      << setw(20) 
		      << int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22) << " "
		      << setw(8) << (double)b12/(double)(b11+b12+b22) << " "
		      << setw(8) << 2 * freq * (1-freq)  << " "
		      << setw(12) << p << " "
		      << "\n";
		}
	      
	      if ( p <= par::HWD_limit ) 
		{
		  cnt++;
		  del[l] = true;
		}
	    }
	  else
	    {
	      // For case/control data
	      
	      double p, p_a, p_u;

	      if (par::HWD_standard)
		{
	      
		  double exp_a11 = afreq * afreq * (a11+a12+a22);
		  double exp_a12 = 2 * afreq * (1-afreq) * (a11+a12+a22);
		  double exp_a22 = (1-afreq) * (1-afreq) * (a11+a12+a22);
		  
		  double exp_u11 = ufreq * ufreq * (u11+u12+u22);
		  double exp_u12 = 2 * ufreq * (1-ufreq) * (u11+u12+u22);
		  double exp_u22 = (1-ufreq) * (1-ufreq) * (u11+u12+u22);
		  
		  double exp_11 = freq * freq * (b11+b12+b22);
		  double exp_12 = 2 * freq * (1-freq) * (b11+b12+b22);
		  double exp_22 = (1-freq) * (1-freq) * (b11+b12+b22);
		  
		  double chisq_a = ( (a11-exp_a11)*(a11-exp_a11) ) / exp_a11 
		    + ( (a12-exp_a12)*(a12-exp_a12) ) / exp_a12 
		    + ( (a22-exp_a22)*(a22-exp_a22) ) / exp_a22 ;
		  
		  double chisq_u = ( (u11-exp_u11)*(u11-exp_u11) ) / exp_u11 
		    + ( (u12-exp_u12)*(u12-exp_u12) ) / exp_u12 
		    + ( (u22-exp_u22)*(u22-exp_u22) ) / exp_u22 ;
		  
		  double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11 
		    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12 
		    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
		  
		  p = chiprobP(chisq,1);
		  p_a = chiprobP(chisq_a,1);
		  p_u = chiprobP(chisq_u,1);
		}
	      else
		{
		  p = SNPHWE( b12, b11, b22 );
		  p_a = SNPHWE( a12, a11, a22 );
		  p_u = SNPHWE( u12, u11, u22 );
		}

	      if (par::HWD_report)
		{
		  HWD << setw(par::pp_maxsnp) << locus[l]->name  << " "
		      << setw(8) << "ALL" << " "
		      << setw(20) 
		      << int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22) << " "
		      << setw(8) << (double)b12/(double)(b11+b12+b22) << " "
		      << setw(8) << 2 * freq * (1-freq)  << " "
		      << setw(12) << p  << "\n";

		  HWD << setw(par::pp_maxsnp) << locus[l]->name  << " "
		      << setw(8) << "AFF" << " "
		      << setw(20) 
		      << int2str(a11)+"/"+int2str(a12)+"/"+int2str(a22) << " "
		      << setw(8) << (double)a12/(double)(a11+a12+a22) << " "
		      << setw(8) << 2 * afreq * (1-afreq)  << " ";
	  
		  if (include_cases)
		    HWD << setw(12) << p_a  << "\n";
		  else
		    HWD << setw(12) << "NA" << "\n";

		  HWD << setw(par::pp_maxsnp) << locus[l]->name  << " "
		      << setw(8) << "AFF" << " "
		      << setw(20) 
		      << int2str(u11)+"/"+int2str(u12)+"/"+int2str(u22) << " "
		      << setw(8) << (double)u12/(double)(u11+u12+u22) << " "
		      << setw(8) << 2 * ufreq * (1-ufreq)  << " ";
		  
		  if (include_cases)
		    HWD << setw(12) << p_u  << "\n";
		  else
		    HWD << setw(12) << "NA" << "\n";		  

		}


	      // Increase counts: in cases
	      if ( include_cases && 
		   p_a < par::HWD_limit && p>-1 ) cnt_a++;

	      // Controls (and, if possible, exclude on this value)
	      if ( include_controls && 
		   p_u < par::HWD_limit && p>-1 ) 
		{
		  cnt_u++;	      

		  if ( ! par::HWD_filter_on_all ) 
		    {
		      del[l] = true;
		      cnt++;
		    }
		}

	      // In total sample, and if needed, exclude here
	      if ( p < par::HWD_limit && p>-1 ) 
		{ 
		  if ( par::HWD_filter_on_all || ! include_controls )
		    {
		      del[l] = true;		  
		      cnt++;
		    }

		}
	      
	    }

	} // next locus
      

      // Finish the report...
      if (par::HWD_report)
	{
	  HWD.close();
	}

      // ...or finish pruning
      stringstream s2;
      s2 <<  cnt 
	 << " markers to be excluded based on HWE test ( p <= " 
	 << par::HWD_limit 
	 << " )\n";
      if (par::bt)
	{
	  s2 << "\t" << cnt_a << " markers failed HWE test in cases\n"
	     << "\t" << cnt_u << " markers failed HWE test in controls\n";
	}
      printLOG(s2.str());
      
    }

  


  ///////////////////////////////////////////////////
  // Summary statistics for genotyping/missing rates

  if (par::report_missing)
    {
      
      printLOG( "Writing individual missingness information to [ " +
		par::output_file_name + ".imiss ] \n");
      
      ofstream MIS;
      string f = par::output_file_name + ".imiss";
      MIS.open(f.c_str(), ifstream::out);
      MIS << setw(par::pp_maxfid) << "FID" << " "
	  << setw(par::pp_maxiid) << "IID" << " " 
	  << setw(10) << "MISS_PHENO" << " "
	  << setw(10) << "N_MISS" << " "
	  << setw(10) << "F_MISS" << "\n";
      
      // By individual
      
      for (int i=0;i<sample.size();i++)
	{
 	  MIS << setw(par::pp_maxfid) << sample[i]->fid << " " 
	      << setw(par::pp_maxiid) << sample[i]->iid << " ";
	  if (sample[i]->missing) MIS << setw(10) << "Y" << " ";
	  else MIS << setw(10) << "N" << " " ;
	  int m=0;
	  for (int l=0; l<locus.size();l++)
	    if ( sample[i]->one[l] && (!sample[i]->two[l]) ) m++;
	  MIS << setw(10) << m << " " 
	      << setw(10) << (double)m/(double)locus.size() << "\n";
	}
      MIS.close();
      
      // By locus
      
      printLOG("Writing locus missingness information to [ " +
	       par::output_file_name +".lmiss ] \n");
      f = par::output_file_name + ".lmiss";
      MIS.open(f.c_str(), ifstream::out);
      MIS.clear();
      
      MIS << setw(4) << "CHR" << " "
	  << setw(par::pp_maxsnp) << "SNP" << " ";
       
      if (par::include_cluster_from_file)
	MIS << setw(10) << "CLST" << " ";

      MIS << setw(8) << "N_MISS" << " ";

      if (par::include_cluster_from_file)
	MIS << setw(8) << "N_CLST" << " ";

      MIS << setw(10) << "F_MISS" << "\n";
      
      vector<Locus*>::iterator loc = locus.begin();          
      int l=0;
      while ( loc != locus.end() )
	{
	  
	  for (int k=0; k<nk; k++)
	    {
	      
	      MIS << setw(4) << (*loc)->chr << " "
		  << setw(par::pp_maxsnp) << (*loc)->name << " ";
	      
	      if (par::include_cluster_from_file)
		MIS << setw(10) << kname[k] << " ";
	      
	      vector<Individual*>::iterator person = sample.begin();
	      int m=0, c=0;
	      while ( person != sample.end() )
		{
		  if (par::include_cluster_from_file)
		    {
		      if ( (*person)->sol == k ) 
			{
			  if ( (*person)->one[l] && ! (*person)->two[l] ) m++;
			  c++;
			}
		    }
		  else
		    {
		      if ( (*person)->one[l] && ! (*person)->two[l] ) m++;
		      c++;
		    }

		  // Next individual
		  (*person)->i1++;
		  (*person)->i2++;
		  person++;
		}

	      MIS << setw(8) << m << " ";
	      if (par::include_cluster_from_file)
		MIS << setw(8) << c << " ";
	      MIS << setw(10) << (double)m / (double)c << "\n";
	    
	    }
	  
	  // Next SNP
	  loc++;
	  l++;
	}
      
      MIS.close();

    }




  /////////////////////////////////
  // Check for being monomorhpic

  // Because this is a vector, it is quicker
  // to copy to required subset into a new vector
  // rather than perform many element erases on 
  // a large vector
  

  // First, remove each individual's genotypes 
  // for to-be-removed markers


  /////////////////////////////////
  // Remove rare SNPs 
  
  vector<Locus*>::iterator loc = locus.begin();
  vector<bool>::iterator d = del.begin();

  while ( loc != locus.end() )
    {
     
      if ( (*loc)->freq<0 || 
	   (*loc)->freq < par::min_af ||
	   (*loc)->freq > par::max_af )        	 
	{
	  *d = true;
	  exc_maf++;
	}
      d++;
      loc++;
    }
  
  
  /////////////////////////////////////////
  // Remove SNPs based on thresholds
    
  if ( locus.size() > 0 ) 
    printLOG("Total genotyping rate in remaining individuals is " + dbl2str(total_genotyping/(double)locus.size())+"\n");
    
  printLOG(int2str(exc_miss)+" SNPs failed missingness test ( GENO > "+dbl2str(par::MAX_GENO_MISSING)+" )\n");
  printLOG(int2str(exc_maf)+" SNPs failed frequency test ( MAF < "+dbl2str(par::min_af));
  if (par::max_af < 0.5 ) printLOG(" or MAF > " + dbl2str(par::max_af));
  printLOG(" )\n");
  
  int tmp = deleteSNPs(del);

}

