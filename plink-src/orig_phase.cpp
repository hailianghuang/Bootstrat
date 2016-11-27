

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
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
#include <cmath>
#include <vector>
#include <map>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"

extern ofstream LOG;

using namespace std;


void HaploPhase::imputeAllHaplotypes()
{

  ///////////////////////////////////////////////
  // Impute missing SNPs -- create a new datafile
  
  // Imputation rules:
  
  // Missing predictor allele -> missing haplotype
  // P(H|G) < 0.8 (default) -> missing haplotype
  

  /////////////////////////////////
  // Phase all specified haplotypes
  
  phaseAllHaplotypes();


  ///////////////////////////
  // Write new PED file

  string filename = par::output_file_name + ".impute.ped";
  
  P.printLOG("Writing imputed ped file to [ " + filename + " ] \n");
  ofstream PED(filename.c_str(), ios::out);
  PED.clear();

  for (int i=0;i<P.n;i++)
    {
      
      Individual * person = P.sample[i];
      PED << person->fid << " "
	  << person->iid << " "
	  << person->pat << " "
	  << person->mat << " "
	  << person->sexcode << " ";
      
      if (par::bt)
	PED << (int)person->phenotype;
      else
	PED << person->phenotype;
      
      for (int l=0;l<new_one[i].size();l++)
	{
	  if ( (!new_one[i][l]) && (!new_two[i][l]) )  
	    PED << par::recode_delimit 
		<< new_map[l]->allele1 << " " 
		<< new_map[l]->allele1;
	  else if ( (!new_one[i][l]) && new_two[i][l]) 
	    PED << par::recode_delimit 
		<< new_map[l]->allele1 << " " 
		<< new_map[l]->allele2;
	  else if (  new_one[i][l]   && new_two[i][l]) 
	    PED << par::recode_delimit 
		<< new_map[l]->allele2 << " " 
		<< new_map[l]->allele2;
	  else PED << par::recode_delimit 
		   << par::missing_genotype << " "
		   << par::missing_genotype;
	}
      PED << "\n";
    }
  
  PED.close();
  
  
  //////////////////
  // Write new map
  
  filename = par::output_file_name + ".impute.map";
  P.printLOG("Writing imputed map file to [ " + filename + " ] \n");
  ofstream MAP(filename.c_str(),ios::out);
  MAP.clear();
  for (int l=0; l<new_map.size(); l++)
    {
      if (new_map[l]->bp >= 0)  // rare SNPs have been removed
	MAP << new_map[l]->chr << "\t"
	    << new_map[l]->name << "\t"
	    << new_map[l]->pos << "\t"
	    << new_map[l]->bp  << "\n";
    }
  MAP.close();
  
}






void HaploPhase::calculateHaplotypeFrequencies()
{

  string f = par::output_file_name + ".frq.hap";
  
  if (par::display_hap_freqs)
    {
      P.printLOG("Writing haplotype frequencies to [ " + f + " ]\n");
      HFRQ.open(f.c_str(), ios::out);
      HFRQ.precision(4);
      
      HFRQ << setw(10) << "LOCUS" << " " 
	   << setw(12) << "HAPLOTYPE" << " "
	   << setw(10) << "F" << "\n";
    }

  // Phase all SNPs (with frequency flag set, this routine
  // will write haplotype frequencies to HFRQ

  P.haplo->phaseAllHaplotypes();
  
  if (par::display_hap_freqs)
    HFRQ.close();
  
  // So we do not re-write them
  par::display_hap_freqs = false;
  
}





void HaploPhase::imputeThisHaplotype(int l)
{
  

  ////////////////////////////////////////////////////
  // Impute for all individiuals, if common haplotype
      
  double w = 0;
  double c = 0;
      
  if (testHaplotypeFreq() >= par::min_hf && 
      testHaplotypeFreq() <= par::max_hf ) 
    {
      
      for (int i=0; i<P.n; i++) 
	{
	  bool b1, b2;	
	  
	  // Haplotype-weighting
	  double t = imputeHaplotypes(i, b1, b2);
	  if (t>=0)
	    {
	      w += t;
	      c++;	  
	    };
	  
	  
	  // And set new imputed genotypes
	  new_one[i].push_back(b1);
	  new_two[i].push_back(b2);
	  
	}
      
      
      // TEMPORAILY REMOVE WEIGHTING FUNCTION
      if (false)
	WGT << new_map[l]->name << "\t"
	    << w/c << "\n";
      
    }
  else
    {    
      // Remove rare haplotypes from map
      new_map[l]->bp = -1;
    }
  
}

void HaploPhase::phaseAllHaplotypes()
{

  //////////////////////////////
  // Begin phasing (and testing)

  int nms = new_map.size();
  
  for (int l=0; l<nms; l++)
    {
      
      if (!par::silent) 
	cout << l+1 << " out of " 
	     << nms << " haplotypes phased       \r";
      
      
      // Set current haplotype
      current = l;
      
      // Get number of predictors
      int s = new_pred_locus[l].size();
      
      // X or haploid markers?
      X = par::chr_sex[new_map[l]->chr];
      haploid = par::chr_haploid[new_map[l]->chr]; 

      // Same window as previous, unless in weighted 
      // multimarker test mode, in which case do all
      bool same = true;
      if (l==0 || par::weighted_mm) same = false;
      else
	{
	  if (new_pred_locus[l].size() != new_pred_locus[l-1].size() ) 
	    { same = false; }
	  else for (int j=0;j<s;j++)
	    {
	      if (new_pred_locus[l][j] != new_pred_locus[l-1][j]) 
		{ same = false; break; }
	    }
	}
            
            
      ///////////////////////////////////////////////////////////////
      // Calculate haplotype frequencies and posterior probabilities,
      // if different from previous round
      
      if (!same)  
	{
	  
	  ///////////////////////////////////////////////////
	  // Initialise
	  
	  reset();

	  for (int i=0; i<P.n; i++)
	    {
	      Individual * person = P.sample[i];
	      if ( person->founder ) 
		{
		  if ( haploid || ( X && person->sex ) ) 
		    validN++;
		  else
		    validN+=2;
		}
	    }

	  // Is this a wildcard haplotype?
	  string nm = new_map[l]->name;
	  if (nm.substr(nm.size()-1) == "_") 
	    name(new_map[l]->name.substr(0,new_map[l]->name.find("_")));
	  else
	    name(new_map[l]->name);

	  	  
	  ///////////////////////////////////////////////////
	  // Generate binary list of all possible haplotypes
	  
	  enumerateHaplotypes(new_pred_locus[l]);

	  /////////////////////////////////////////////////
	  // Enumerate all possible phases for each founder
	  
	  for (int i=0; i<P.n; i++) 
	    if (P.sample[i]->founder)
	      enumeratePhase(i);
	  
	  
	  /////////////////////////////////////
	  // Use offspring to reduce ambiguity
	  
  // To be added in; not yet implemented //
// 	  if (nonfounders)
//  	    for (int i=0; i<P.n; i++) 
// 	      if (P.sample[i]->founder && ambig[i])
// 		resolveWithChildren(i)
		  
	  

	  ////////////////////////////////
	  // E-M phasing based on founders
	  
	  if ( par::phase_old_EM)
	    performEM_original();
	  else
	    performEM();

	  ////////////////////////
	  // Prune unlikely phases 

 	  for (int i=0; i<P.n; i++) 
 	    if (P.sample[i]->founder)
 	      prunePhase(i);


	  ///////////////////////////////////
	  // Fill-in phasing any non-founders
	  
	  if (nonfounders)
	  {

	    if ( haploid ) 
	      error("Family-based haplotyping only for autosomes and X chromosome");

	    // List all possible non-rare phases
//  	    enumerateAllPhases();
	    
	    for (int i=0; i<P.n; i++) 
	      if (!P.sample[i]->founder)
		{
		  // Only phase affecteds, if performing TDT
		  if ( par::test_hap_TDT || par::proxy_TDT )
		    if ( !P.sample[i]->aff )
		      continue;
		  
		  phaseAndScoreNonfounder(i);
		  prunePhase(i);
		}
	  } 	      

	  ////////////////////////
	  // Post-phasing actions?
	  
	  if (par::display_hap_freqs)
	    reportHaplotypeFrequencies();
	  
	  if (par::test_hap_CC || par::test_hap_TDT || par::test_hap_QTL ) 
	    {
	      
	      if ( ! par::phase_hap_all )
		setTestHaplotype(new_pred_allele[l]);
	      
	      vector_t tmp = performHaplotypeTests();
	    }
	  
	  if (par::display_phase_probs) 
	    {
	      if ( par::display_phase_probs_wide )
		reportPhaseWideFormat();
	      else
		reportPhase();
	    }
	  
	} // next unique window



      if (par::impute_tags) 
	{
	  setTestHaplotype(new_pred_allele[l]);
	  imputeThisHaplotype(l);
	}

    }
  
  if (!par::silent) 
    cout << "\n";


}



void HaploPhase::enumerateHaplotypes(vector<int> & s)
{  
  
  // Make list of haplotypes, and code as +/-
  
  S = s;
  ns = s.size();
  nh = (int)pow((double)2,ns);
  f.resize(nh);

  ph_hap1.clear();
  ph_hap2.clear();
  ph_freq.clear();
  
  // Optionally, for haplotype-based TDT
  if (par::test_hap_TDT || par::proxy_TDT )
    {
      trans.clear();
      untrans.clear();

      trans.resize(nh,0);
      untrans.resize(nh,0);
    }
  
  unsigned int h=0;
  
  while(h<nh)
    {
      vector<bool> m1;
      
      unsigned int p=1;
      for (int s=0;s<ns;s++)
	{
	  if ( h & p ) m1.push_back(false);
	  else m1.push_back(true);
	  p <<= 1;
	}
      
      // Add to list of haplotypes
      hap.push_back(m1);
      
      // Add to HapMap
      hapmap.insert(make_pair(m1,h));
      
      // Consider next haplotype
      h++;     

   }

  // Uniform starting values
  // for (int i=0;i<nh;i++) 
  //   f[i] = 1/(double)nh;

  // Product of allele frequencies
  double psum = 0;
  for (int h=0;h<nh;h++) 
    {
      f[h] = 1;
      for (int s=0; s<ns; s++)
	{
	  if ( !hap[h][s] ) 
	    f[h] *= P.locus[S[s]]->freq;
	  else
	    f[h] *= ( 1 - P.locus[S[s]]->freq );
	}
      psum += f[h];
  }

}


void HaploPhase::enumerateAllPhases()
{  

  // Note: this function is no longer called
  
  // For individuals w/out parents: make a list of all possible
  // phases. Note: currently, we do not use this (i.e. we always
  // require father and mother to be 'observed' (i.e. genotyped and
  // adequately phased).

  // Also note: issue with representing heterozygote haplotypes twice
  // in list: previously we did not, but we now change this (seeing as
  // it is never used in any case...)

  // Also: now we build separate lists for diploid and haploid
  // chromosomes, not that we use either.
  
  // Diploid possible phases

  if ( ! haploid  ) 
    {
      for (int h1=0;h1<nh;h1++)
	for (int h2=0;h2<nh;h2++)
	  {
	    double freq = f[h1] * f[h2];
	    if (h1!=h2) freq *= 2;
	    if ( freq >= par::hap_min_phase_prob )
	      {
		ph_freq.push_back( freq );
		ph_hap1.push_back( h1 );
		ph_hap2.push_back( h2 );		
	      }
	  }        
    }

  // Haploid possible phases

  if ( haploid || X ) 
    {
      for (int h1=0;h1<nh;h1++)
	{
	  double freq = f[h1];
	  if ( freq >= par::hap_min_phase_prob )
	    {
	      haploid_ph_freq.push_back( freq );
	      haploid_ph_hap1.push_back( h1 );
	    }
	} 
    }


  // Original total number of phases 

  np = ph_hap1.size();
  haploid_np = haploid_ph_hap1.size();

}



vector<string> HaploPhase::returnHaplotypes(vector<int> & slist)
{
  vector<string> str;
  
  enumerateHaplotypes(slist);
  
  for (int h=0; h<hap.size(); h++)     
    {
      string hstr;
      for (int s=0;s<ns;s++)
        if (!hap[h][s]) hstr += P.locus[S[s]]->allele1;
        else if (P.locus[S[s]]->allele2=="") hstr += "0";
	else hstr += P.locus[S[s]]->allele2;
      str.push_back(hstr);
    }
  return str;
}



void HaploPhase::setTestHaplotype(string t)
{

  // Create match template
  vector<bool> tmp(ns,false); 
  for (int s=0;s<ns;s++)
    if ( P.locus[S[s]]->allele1 == t.substr(s,1) )
      tmp[s] = true;

  // Consider each haplotype 
  
  test_hap = -1; 
  
  for (int h=0; h<hap.size(); h++)     
    {
      bool match = true;
      
      for (int s=0;s<ns;s++)
	{
	  if ( hap[h][s] != tmp[s] )
	    {
	      match = false;	  	  
	      break;
	    }
	}
      
      if (match) 
	{
	  test_hap = h;
	  break;
	}
    }

}


void HaploPhase::prunePhase(int i)
{
  
  // pp[i][z]
  // hap1[i][z]
  // hap2[i][z]
  
  if ( (!include[i]) || 
       (!ambig[i])   ) return;

  double psum = 0;
  
  vector<double> new_pp(0);
  vector<int> new_h1(0);
  vector<int> new_h2(0);
  
  for (int z=0; z < hap1[i].size(); z++)
    {
      if (pp[i][z] >= par::hap_min_phase_prob)
	{
	  new_pp.push_back( pp[i][z] );
	  psum += pp[i][z];
	  new_h1.push_back( hap1[i][z] );
	  new_h2.push_back( hap2[i][z] );
	}
    }
  
  // Normalise?
  if ( pp[i].size() > new_pp.size() )
    {
      for (int z=0; z < new_pp.size(); z++)
	new_pp[z] /= psum;
    }
  
  // Update
  pp[i] = new_pp;
  hap1[i] = new_h1;
  hap2[i] = new_h2;

}



void HaploPhase::enumeratePhase(int i)
{  

  vector<bool> s1(ns);
  vector<bool> s2(ns);

  hap1[i].clear();
  hap2[i].clear();
 
  // Flipping allele-coding for homozygotes

  for (int s=0; s<ns; s++)
    {
      if (par::SNP_major)
	{	
	  s1[s] = P.SNP[S[s]]->one[i];
	  s2[s] = P.SNP[S[s]]->two[i];
	}
      else
	{
	  s1[s] = P.sample[i]->one[S[s]];
	  s2[s] = P.sample[i]->two[S[s]];
	}
      
      if ( s1[s] == s2[s] )
	{
	  s1[s] = !s1[s];
	  s2[s] = !s2[s];
	}
      
    }
  
  
  //////////////////////////////////////////////////////////
  // Count amount of missing genotype data at this position

  int nm = 0;
  for (int s=0; s<ns; s++)
    if ( s1[s] && !s2[s] )
      nm++;

	
  // If any missing genotypes, this person counts 
  // as ambiguous

  if (nm>0) ambig[i] = true;  


  // But if too much missing genotype data, then 
  // we should not even try to phase this individual
  // for this region; note -- females should always be
  // missing all genotypes for Y, so we don't need to 
  // worry about allowing for a special case here.

  if ( (double)nm/(double)ns >= par::hap_missing_geno ) 
    {
      include[i] = false;
      
      if ( P.sample[i]->founder) 
	{
	  if ( haploid || ( X && P.sample[i]->sex ) ) 
	    validN--;
	  else
	    validN-=2;
	}

      return;
    }


  ///////////////////////////////////////////////
  // 2 or more hets at any loci -> ambiguous
  // Haploid genotypes should never be heterozygous, 
  // so we are okay here w.r.t X chromosome

  int het=0;
  for (int s=0; s<ns; s++)
    if ( (!s1[s]) && s2[s] ) het++;
  if (het>1) ambig[i] = true;
  

  //////////////////////////////////////
  // Construct list of consistent phases
  
  if (!ambig[i])
    {		
      // Unambiguous means all no missing genotypes
      // and less than 2 hets
      
      // Match haplotype alleles: haploid individuals
      // will just be coded as homozygous here (but 
      // when considering phases, frequencies, etc
      // we will take care of this downstream)

      hap1[i].push_back( hapmap.find(s1)->second );
      hap2[i].push_back( hapmap.find(s2)->second );
            
    }
  else
    {
      
      // For individuals with ambiguity 
      
      // Which are the ambiguous sites
      // (missing or heterozygous)

      // We will not observe any hets for haploid 
      // individuals: but we do need to make sure 
      // that missing haploid genotypes are not 
      // allowed to be heterozygous

      vector<bool> het_site(ns,false);
      vector<bool> mis_site(ns,false);

      int firstHeterozygote = ns;
            
      int ambig_cnt=0;
      for (int s=0; s<ns; s++)
	{
	  // het
	  if ( (!s1[s]) && s2[s] )
	    {
	      het_site[s] = true;
	      
	      if ( firstHeterozygote == ns )
		firstHeterozygote = s;
	      
	      ambig_cnt++;
	      
	    }
	  
	  // missing
	  if ( s1[s] && (!s2[s]) )
	    {

	      mis_site[s] = true;
	      		      
	      // haploid:  0- or 1-
	      // diplod :  00 or 01 or 11
	      
	      if ( haploid || ( X && P.sample[i]->sex ) )
		ambig_cnt++;
	      else
		ambig_cnt+=2;
	      
	    }
	}


      int ambig_nh = (int)pow((double)2,ambig_cnt);

      vector<bool> h1(ns);
      vector<bool> h2(ns);
 
      int original_firstHeterozygote = firstHeterozygote;
	  
      int h=0;
      while(h<ambig_nh)
	{

	  vector<bool> m1;
	  
	  unsigned int p=1;
	  for (int s=0;s<ambig_cnt;s++)
	    {
	      if ( h & p ) m1.push_back(false);
	      else m1.push_back(true);
	      p <<= 1;
	    }

	  // Splice m1-variant into h1, h2 to 
	  // reconstruct ambiguous sites

	  int ac=0;
	  bool skip = false;

	  // If missing homozygote imputed, the next 
	  // het can be fixed possibly; also, if not 
	  // already done, the first missing genotype
	  // can be fixed
	  
	  int firstHeterozygote = original_firstHeterozygote ;

	  for (int s=0; s<ns; s++)
	    {
	      
	      // Reconstruct heterozygous sites
	      // (except the first one)
	      
	      if (het_site[s])
		{
		  
		  if (m1[ac]) 
		    {
		      h1[s] = true;
		      h2[s] = false;
		    }
		  else
		    {
		      if ( s != firstHeterozygote ) 
 			{
			  h1[s] = false;
			  h2[s] = true;
			}
		      else
			{
			  h1[s] = false;
			  h2[s] = true;
			  
			  skip = true;			  
			}
		    }
		  
		  ac++;

		}
	      else if (mis_site[s])
		{
		  // Treat haploid and diploid differently for missing
		  // genotype data
		  
		  if ( haploid || ( X && P.sample[i]->sex ) )
		    {
		      // select a hemizygote/homozygote
		      if (m1[ac])
			h1[s] = h2[s] = false;
		      else 
			h1[s] = h2[s] = true;
		      
		      ac++;
		    }
		  else
		    {
		      // Make het
		      if (m1[ac]) 
			{			  			  
			  
			  if (m1[ac+1])
			    {
			      h1[s] = false;
			      h2[s] = true;
			    }
			  else
			    {
			      if ( s < firstHeterozygote ) 
				{				  
				  skip = true;
				}				  
			      else
				{
				  h1[s] = true;
				  h2[s] = false;
				}
			    }

			  firstHeterozygote = s;
			  
			}
		      else 
			{
			  
			  // otherwise, select a homozygote
			  if (m1[ac+1])
			    h1[s] = h2[s] = false;
			  else 
			    h1[s] = h2[s] = true;
			}
		      
		      ac+=2;

		    }
		}
	      else
		{
		  // Maintain unambigous site
		  // (which might be 1st het)
		  h1[s] = s1[s];
		  h2[s] = s2[s];
		}
	    }

	  if ( ! skip ) 
	    {

	      // Add to (non-redundant) list?
	      
	      int n1 = hapmap.find( h1 )->second;
	      int n2 = hapmap.find( h2 )->second;
	  
	      hap1[i].push_back( n1 );
	      hap2[i].push_back( n2 );	      
	      
	    }

	  // Consider next haplotype pair
	  h++;     
	  
	}
    }
  
 
  // Make space for posterior probabilities
  
  if (ambig[i]) pp[i].resize(hap1[i].size());

}

void HaploPhase::reportPhase()
{
  string fn = par::output_file_name+".phase-"+hname;
  ofstream PHASE(fn.c_str(), ios::out);
  
  P.printLOG("Writing phased haplotypes for " 
	     + hname + " to [ " + fn + " ]\n");
  
  PHASE << setw(12) << "FID" << " "
	<< setw(12) << "IID" << " "
	<< setw(4)  << "PH" << " "
	<< setw(10) << "HAP1" << " "
	<< setw(10) << "HAP2" << " "
	<< setw(12) << "POSTPROB" << " "
    //  << setw(12) << "WEIGHT" << " "
	<< setw(6) << "BEST" << " "
	<< "\n";
  
  PHASE.precision(4);
  
  for (int i = 0 ; i < P.n ; i++ ) 
    {
      
      if (include[i])
	{
	  
	   for (int z = 0 ; z < hap1[i].size(); z++) 
	     {
	       
	       PHASE << setw(12) << P.sample[i]->fid << " " 
		     << setw(12) << P.sample[i]->iid << " "
		     << setw(4) << z << " "
		     << setw(10) << haplotypeName(hap1[i][z])  << " ";
	       if ( haploid || ( X && P.sample[i]->sex ) )
		 PHASE << setw(10) << haplotypeName( -1 ) << " ";
	       else
		 PHASE << setw(10) << haplotypeName(hap2[i][z]) << " ";

	       if (ambig[i])	
		 {	
		   PHASE << setw(12) << pp[i][z] << " ";
		   int max_z = 0;
		   for (int z2=0; z2<hap1[i].size(); z2++)
		     max_z = pp[i][z2] > pp[i][max_z] ? z2 : max_z ;

 // 		  int fac=1;
 // 		  if ( hap1[i][max_z] != hap2[i][max_z] ) fac=2;
 // 		  double w = ( pp[i][max_z] - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] ) 
 // 		    / ( 1 - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] );

		   // Do not output weight for now
		   // if (max_z==z) PHASE << setw(12) << w << " " << setw(6) << 1 << " " << "\n";
		   // else PHASE << setw(12) << "."  << " " << setw(6) << 0  << " " << "\n";

		   if (max_z == z) PHASE << setw(6) << 1  << " " << "  ";
		   else PHASE << setw(6) << 0  << " " << "  ";
		 }
	       else
		 PHASE << setw(12) << 1  << " "
		       << setw(6) << 1  << " "
		       << "  "; 

	       // Genotypes
	       for (int s=0; s<ns; s++)
		 PHASE << genotype(P,i,S[s]) << " ";
	       PHASE << "\n";


	     }
	 }

       // Report also on excluded individuals
       // (Should be 0-size phase-set)
      else 
	{
	  
	  PHASE << setw(12) << P.sample[i]->fid << " " 
		<< setw(12) << P.sample[i]->iid << " "
		<< setw(4) << "NA" << " "
		<< setw(10) << "NA"  << " "
		<< setw(10) << "NA" << " "
		<< setw(12) << "NA"  << " "
		<< setw(6) << "NA"  << "   ";
	  
	  // genotypes
	  for (int s=0; s<ns; s++)
	    PHASE << genotype(P,i,S[s]) << " ";
	  PHASE << "\n";
	  
	  
	}
    }
  
  PHASE.close();
  
}  


void HaploPhase::reportPhaseWideFormat()
{
  string fn = par::output_file_name+".wphase-"+hname;
  ofstream PHASE(fn.c_str(), ios::out);
  
  P.printLOG("Writing wide-format phased haplotypes for " 
	     + hname + " to [ " + fn + " ]\n");
  
  PHASE << setw(par::pp_maxfid) << "FID" << " "
	<< setw(par::pp_maxiid) << "IID" << " ";
  
  for (int h=0; h<nh; h++)
    if (f[h] >= par::min_af)
      PHASE << setw(8) << "H_"+haplotypeName(h) << " ";
  
  PHASE << "\n";
  
  PHASE.precision(4);
  
  
  for (int i = 0 ; i < P.n ; i++ ) 
    {
      if (include[i])
	{
	  
	  PHASE << setw(par::pp_maxfid) << P.sample[i]->fid << " " 
		<< setw(par::pp_maxiid) << P.sample[i]->iid << " ";
	  
	  vector_t hcnt(nh,0);
	  
	  for (int z = 0 ; z < hap1[i].size(); z++) 
	    {
	      
	      if (ambig[i])
		{
		  hcnt[hap1[i][z]] += pp[i][z];
		  if ( ! ( haploid || ( X && P.sample[i]->sex ) ) )
		    hcnt[hap2[i][z]] += pp[i][z];
		}
	      else
		{
		  hcnt[hap1[i][z]] ++;
		  if ( ! ( haploid || ( X && P.sample[i]->sex ) ) )
		    hcnt[hap2[i][z]] ++;
		}
	      
	    }
	  
	  for (int h=0; h<nh; h++)      
	    if (f[h] >= par::min_af)
	      PHASE << setw(8) << hcnt[h] << " ";
	  PHASE << "\n";
	  
	}
      
      // Report also on excluded individuals
      // (Should be 0-size phase-set)
      
      else 
	{
	  PHASE << setw(par::pp_maxfid) << P.sample[i]->fid << " " 
		<< setw(par::pp_maxiid) << P.sample[i]->iid << " ";
	  
	  for (int h=0; h<nh; h++)      
	    if (f[h] >= par::min_af)
	      PHASE << setw(8) << "NA" << " ";
	  PHASE << "\n";
	}
    }
  
   PHASE.close();

 }  


class MultiLocusGenotype
{
public:

  vector<bool> g;
  int count;
  int reference;
  vector<bool> skip;

  bool operator< (const MultiLocusGenotype & b) const
  {
    for (int i=0; i<g.size(); i++)
      if ( g[i] != b.g[i] ) 
	return g[i];
    return false;
  }

  bool operator== (const MultiLocusGenotype & b) const
  {
    for (int i=0; i<g.size(); i++)
      if ( g[i] != b.g[i] ) 
	return false;
    return true;
  }

};

namespace std
{
  template<>
  class less<MultiLocusGenotype*> {
  public:
    bool operator()(MultiLocusGenotype const* p1, MultiLocusGenotype const* p2)
    {     
      if(!p1)
	return true;
      if(!p2)
	return false;
      
      if (p1->g < p2->g)
	return true;
      
      return false;
    }
  };  
};


void HaploPhase::performEM()
{
  

  // Store count, reference individual
  set<MultiLocusGenotype*> genotypes;
  
  // Store which genoGroup a person belongs to
  vector<MultiLocusGenotype*> genoGroup(P.n, (MultiLocusGenotype*)0);
  
  // Make groups
  for (int i=0; i<P.n; i++)
    {

      // Do we want to consider this individual?
      
      if ( ! ( P.sample[i]->founder && include[i]) )
	continue;
	
      // Build a new multilocus genotype set

      MultiLocusGenotype * m = new MultiLocusGenotype;

      // Include sex here for X chr SNPs

      if ( X ) 
	m->g.push_back( P.sample[i]->sex );

      // Genotypes

      for (int s=0; s<ns; s++)
	{
	  m->g.push_back( P.SNP[ S[s] ]->one[i]); 
	  m->g.push_back( P.SNP[ S[s] ]->two[i]); 
	}
      
      // One individual, this individual

      m->count = 1;
      m->reference = i;
      m->skip.resize( hap1[i].size() , false );

      // But have we already seen a similar genoGroup?

      set<MultiLocusGenotype*>::iterator im = genotypes.find(m);
      
      if ( im  == genotypes.end() )
	{	  
	  genoGroup[i] = m;
	  genotypes.insert( m );
	}
      else
	{
	  delete m;
	  (*im)->count++;	  
	  genoGroup[i] = *im;
	}
    }
  
  
  // Use the Individual 'flag' to indicate whether or not this individual
  // is assigned to a genoGroup

  for (int i=0; i<P.n; i++)
    {

      if ( ! ( P.sample[i]->founder && include[i]) )
	continue;
      
      if ( genoGroup[i]->count >= 1 ) 
 	P.sample[i]->flag = false; // no individual assignment
      else
 	P.sample[i]->flag = true;
    }
  
//   // Remove small groups

//   set<MultiLocusGenotype*>::iterator im = genotypes.begin();
//   while ( im != genotypes.end() )
//     {      
//       if( (*im)->count < 10 )
// 	{
// 	  delete *im;	  
// 	  ++im;
// 	  genotypes.erase(im-1);
// 	}
//       else
// 	++im;
//     }


   vector<double> uc(nh,0);  // unambigous counts
   vector<double> ac(nh,0);  // ambigous counts
   

   ///////////////////////////////////////////
   // Count numbers of unambigous haplotypes
   // as these stay constant throughout EM
   
   //////////////
   // genoGroups 

   set<MultiLocusGenotype*>::iterator im = genotypes.begin();
   while ( im != genotypes.end() )
     {
       int i = (*im)->reference; 
       if (!ambig[i])
	 {
	   uc[hap1[i][0]] += (*im)->count;
	   if ( ! ( haploid || ( X && P.sample[i]->sex ) ) )
	     uc[hap2[i][0]] += (*im)->count;
	 }
       ++im;
     }

   ///////////////
   // Individuals 
   
//    for (int i=0; i<P.n; i++)
//      {
//        if ( P.sample[i]->flag ) 
// 	 {
// 	   if (P.sample[i]->founder && include[i])
// 	     {
// 	       if (!ambig[i])
// 		 {
// 		   uc[hap1[i][0]]++;
// 		   if ( ! ( haploid || ( X && P.sample[i]->sex ) ) )
// 		     uc[hap2[i][0]]++;
// 		 }
// 	     } 
// 	 }
//      }


   //////////////////  
   // Begin E-M
   
   vector<bool> zero( nh, false );
   int zeroed = 0;

   double sampleLogLikelihood = 0;
   int it = 0;

   for (int j=0; j<=2000; j++)
     {

       //////////////////////////
       // E-step for genoGroups

       im = genotypes.begin();
       while ( im != genotypes.end() )
	 {
	   int i = (*im)->reference; 
	   if (ambig[i])
	     {
	       double s=0;
	       // Haploid phases...
	       if ( haploid || ( X && P.sample[i]->sex ) )
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     {
		       pp[i][z] = f[hap1[i][z]];
		       s += pp[i][z];
		     }
		 }
	       else // ... or diploid
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     {
		       
// 		       if ( (*im)->skip[z] )
// 			 continue;
		       
		       int h1 = hap1[i][z];
		       int h2 = hap2[i][z];
		       
//  		       if ( zero[h1] || zero[h2] )
//  			 {
//  			   (*im)->skip[z] = true;
//  			   continue;
// 			 }
		       
		       pp[i][z] = f[h1] * f[h2];
		       if (h1 != h2) pp[i][z] *= 2;
		       
		       s += pp[i][z];
		     }
		 }
	       
	       for (int z=0; z<hap1[i].size(); z++) pp[i][z] /= s;	       

	     }
	   ++im;
	 }
       
       //////////////////////////
       // E-step for individuals
       
//        for (int i=0; i<P.n; i++)
// 	 {
// 	   if (P.sample[i]->flag && P.sample[i]->founder && include[i])
// 	     {
// 	       if (ambig[i])
// 		 {

// 		   double s=0;

// 		   // Haploid phases...
// 		   if ( haploid || ( X && P.sample[i]->sex ) )
// 		     {
// 		       for (int z=0; z<hap1[i].size(); z++)
// 			 {
// 			   pp[i][z] = f[hap1[i][z]];
// 			   s += pp[i][z];
// 			 }
// 		     }
// 		   else // ... or diploid
// 		     {
// 		       for (int z=0; z<hap1[i].size(); z++)
// 			 {
// 			   pp[i][z] = f[hap1[i][z]] * f[hap2[i][z]];
// 			   if (hap1[i][z] != hap2[i][z]) pp[i][z] *= 2;
// 			   s += pp[i][z];
// 			 }
// 		     }
		   
// 		   for (int z=0; z<hap1[i].size(); z++) pp[i][z] /= s;
// 		 }
// 	     }
//	 }


       /////////////////////////////////////
       // M-step for pre-counted haplotypes

       // unambiguous counts
       for (int h=0; h<nh; h++)
	 f[h] = uc[h];

       ////////////////////////////////////
       // M step for ambiguous genoGroups

       im = genotypes.begin();
       while ( im != genotypes.end() )
	 {
	   int i = (*im)->reference; 
	   if (ambig[i])
	     {
	       if ( haploid || ( X && P.sample[i]->sex ) )
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     f[hap1[i][z]] += pp[i][z] * (*im)->count;
		 }
	       else
		 for (int z=0; z<hap1[i].size(); z++)
		   {

 		     if ( (*im)->skip[z] )
 		       continue;

		     f[hap1[i][z]] += pp[i][z] * (*im)->count;
		     f[hap2[i][z]] += pp[i][z] * (*im)->count;		  
		   }
	     }
	   ++im;
	 }
       

       ////////////////////////////////////
       // M step for ambiguous individuals

//        for (int i=0; i<P.n; i++)
// 	 if (P.sample[i]->flag && P.sample[i]->founder && include[i])
// 	   if (ambig[i])
// 	     {

// 	       if ( haploid || ( X && P.sample[i]->sex ) )
// 		 {
// 		   for (int z=0; z<hap1[i].size(); z++)
// 		     f[hap1[i][z]] += pp[i][z];
// 		 }
// 	       else
// 		 for (int z=0; z<hap1[i].size(); z++)
// 		   {
// 		     f[hap1[i][z]] += pp[i][z];
// 		     f[hap2[i][z]] += pp[i][z];		  
// 		   }
//	     }

       // validN is the total number of *chromosomes*
       for (int h=0; h<nh; h++)      
	 f[h] /= (double)validN;


       /////////////////////
       // Update likelihood (every 5th iteration)

       if ( j % 5 == 0 )
	 {
	   
//  	   for (int h=0; h<nh; h++) 
//  	     if ( ! zero[h] )
//  	       {
//  		 if ( f[h] == 0 ) {
//  		   zeroed++;
//  		   zero[h] = true;
//  		 }
// 	       }

	   double lnl = 0;
	   
	   // genoGroups

	   im = genotypes.begin();
	   while ( im != genotypes.end() )
	     {
	       int i = (*im)->reference; 
	       
	       double lk = 0;
	       
	       if ( haploid || ( X && P.sample[i]->sex ) )
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     lk += f[hap1[i][z]];
		 }
	       else
		 for (int z=0; z<hap1[i].size(); z++)
		   {
 		     if ( (*im)->skip[z] )
 		       continue;

		     lk += f[hap1[i][z]] * f[hap2[i][z]] ;
		     if ( hap1[i][z] != hap2[i][z] )
		       lk += f[hap1[i][z]] * f[hap2[i][z]] ;
		   }

	       lnl -= log(lk) * (*im)->count ;
	       
	       ++im;
	     }
	   
	   
	   // And individuals
	   
// 	   for (int i=0; i<P.n; i++)
// 	     {
	       
// 	       if (P.sample[i]->flag && P.sample[i]->founder && include[i])
// 		 {
// 		   double lk = 0;
		   
// 		   if ( haploid || ( X && P.sample[i]->sex ) )
// 		     {
// 		       for (int z=0; z<hap1[i].size(); z++)
// 			 lk += f[hap1[i][z]];
// 		     }
// 		   else
// 		     for (int z=0; z<hap1[i].size(); z++)
// 		       {
// 			 lk += f[hap1[i][z]] * f[hap2[i][z]] ;
// 			 if ( hap1[i][z] != hap2[i][z] )
// 			   lk += f[hap1[i][z]] * f[hap2[i][z]] ;
// 		       }
// 		   lnl -= log(lk);
// 		 }
//	     }
	   
	   if ( j > 0 && sampleLogLikelihood - lnl < 0.0001 )
	     {	   
	       it = j;
	       break;
	     }
	   
	   sampleLogLikelihood = lnl;
	 }


     }


   for (int j=0; j<P.n; j++)
     {

       if ( ! ( P.sample[j]->founder && include[j]) )
	 continue;

       if ( !P.sample[j]->flag ) 
	 {

	   if (ambig[j])
	     {
	       
	       MultiLocusGenotype * m = genoGroup[j];
	       int i = m->reference;
	       if ( i != j )
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     {
		       pp[j][z] = pp[i][z];	       
		     }
		 }
	     }
	 }
     }


   // Free genoGroups
   
   im = genotypes.begin();
   while ( im != genotypes.end() )
     {
       delete *im;
       ++im;
     }

   
}


///////////////////////////////////////////////////
// Legacy version of algorithm without genoGrouping

void HaploPhase::performEM_original()
 {

   vector<double> uc(nh,0);  // unambigous counts
   vector<double> ac(nh,0);  // ambigous counts

   // Count numbers of unambigous haplotypes
   // as these stay constant throughout EM

   for (int i=0; i<P.n; i++)
     {
       if (P.sample[i]->founder && include[i])
	 {
	   if (!ambig[i])
	     {
	       uc[hap1[i][0]]++;
	       if ( ! ( haploid || ( X && P.sample[i]->sex ) ) )
		 uc[hap2[i][0]]++;
	     }
	 } 
     }



   //////////////////  
   // Begin E-M


   double sampleLogLikelihood = 0;

   for (int j=0; j<=1000; j++)
     {

	     
       ///////////
       // E-step

       for (int i=0; i<P.n; i++)
	 {
	   if (P.sample[i]->founder && include[i])
	     {
	       if (ambig[i])
		 {

		   double s=0;

		   // Haploid phases...
		   if ( haploid || ( X && P.sample[i]->sex ) )
		     {
		       for (int z=0; z<hap1[i].size(); z++)
			 {
			   pp[i][z] = f[hap1[i][z]];
			   s += pp[i][z];
			 }
		     }
		   else // ... or diploid
		     {
		       for (int z=0; z<hap1[i].size(); z++)
			 {
			   pp[i][z] = f[hap1[i][z]] * f[hap2[i][z]];
			   if (hap1[i][z] != hap2[i][z]) pp[i][z] *= 2;
			   s += pp[i][z];
			 }
		     }
		   
		   for (int z=0; z<hap1[i].size(); z++) pp[i][z] /= s;
		 }
	     }
	 }



       ///////////
       // M-step

       // unambiguous counts
       for (int h=0; h<nh; h++)
	 f[h] = uc[h];

       // then add the fractional ones
       for (int i=0; i<P.n; i++)
	 if (P.sample[i]->founder && include[i])
	   if (ambig[i])
	     {

	       if ( haploid || ( X && P.sample[i]->sex ) )
		 {
		   for (int z=0; z<hap1[i].size(); z++)
		     f[hap1[i][z]] += pp[i][z];
		 }
	       else
		 for (int z=0; z<hap1[i].size(); z++)
		   {
		     f[hap1[i][z]] += pp[i][z];
		     f[hap2[i][z]] += pp[i][z];		  
		   }
	     }

       // validN is the total number of *chromosomes*
       for (int h=0; h<nh; h++)      
	 f[h] /= (double)validN;


       /////////////////////
       // Update likelihood (every 5th iteration)

       if ( j % 5 == 0 )
       if ( true )
	 {

	   double lnl = 0;
	   
	   for (int i=0; i<P.n; i++)
	     {
	       
	       if (P.sample[i]->founder && include[i])
		 {
		   double lk = 0;
		   
		   if ( haploid || ( X && P.sample[i]->sex ) )
		     {
		       for (int z=0; z<hap1[i].size(); z++)
			 lk += f[hap1[i][z]];
		     }
		   else
		     for (int z=0; z<hap1[i].size(); z++)
		       {
			 lk += f[hap1[i][z]] * f[hap2[i][z]] ;
			 if ( hap1[i][z] != hap2[i][z] )
			   lk += f[hap1[i][z]] * f[hap2[i][z]] ;
		       }
		   lnl -= log(lk);
		 }
	     }
	   
	   if ( j > 0 && sampleLogLikelihood - lnl < 1e-4 )
	     {	   
	       break;
	     }
	   
	   sampleLogLikelihood = lnl;
	 }
       
     }
}


void HaploPhase::reportHaplotypeFrequencies()
{ 
  for (int h=0; h<nh; h++)      
    {
      if (f[h] >= par::min_af)
	{
	  HFRQ << setw(10) << hname << " "
	       << setw(12) << haplotypeName(h) << " " 
	       << setw(10) << f[h] << "\n";
	}
    }
}


map<int,int> HaploPhase::makeTestSet(boolvec_t & mask, boolvec_t & allele)
{
  map<int,int> tests;
	      
  for (int h2=0; h2 < nh; h2++)
    {
      bool is_A = true;
      for (int s = 0; s < ns ; s++)
	{
	  if ( mask[s] && hap[h2][s] != allele[s] )
	    is_A = false;
	}
      
      if ( is_A ) 
	tests.insert(make_pair(h2,0));
      else 
	tests.insert(make_pair(h2,1));
    }

  return tests;
}

string HaploPhase::getSubHaplotypeName(boolvec_t & mask, boolvec_t & allele, int blank)
{
  string str = "";

  for (int s=0; s < ns; s++)
    {
      if ( s == blank )
	str += " ";
      else if ( mask[s] )
	{
	  if ( allele[s] ) str += P.locus[ S[s] ]->allele1;
	  else str += P.locus[ S[s] ]->allele2;	      
	}
      else
	str += ".";
    }
  return str;
}

vector_t HaploPhase::performHaplotypeTests()
{

  vector_t v(0);
  
  // Weighted multimarker test (test a single Haplotype)
  
  if ( par::weighted_mm )
    {
      
      if (par::test_hap_CC)
	haplotypicWeightedCC();
      else if (par::test_hap_TDT)
	haplotypicWeightedTDT();
      
      return v;
    }
  
  
  // Standard haplotype tests
  
  if ( !par::phase_hap_all ) 
    {
      
      map<int,int> tests;

      // Of the specific, prespecified haplotype? 
      // Test only the specific haplotype
      
      if (f[test_hap] >= par::min_af)
	{
	  
	  for (int h2=0; h2 < nh; h2++)
	    {
	      if (f[h2] >= par::min_af)
		{
		  if (test_hap==h2) 
		    tests.insert(make_pair(h2,0));
		  else 
		    tests.insert(make_pair(h2,1));
		}
	    }
	  
	  if (par::test_hap_CC)
	    haplotypicCC(tests,2,true);
	  else if (par::test_hap_QTL)
	    haplotypicQTL(tests,2,true);
	  else if (par::test_hap_TDT)
	    haplotypicTDT(tests,2,true);
	  
	}
      
    }
  else
    {
      // Perform an omnibus test and all haplotype-specific
      // tests, if the --hap-all flag has been set
  
      // Omnibus test
      
      // tests[0] = 0
      // tests[1] = 1
      // tests[2] = 2
      // ...
      // tests[h] = h
      
      map<int,int> tests;
      int nch=0;
      for (int h=0; h < nh; h++)
	if (f[h] >= par::min_af)
	  {
	    tests.insert(make_pair(h,nch++));
	  }
      
      if (nch>2)
	{
	  if (par::test_hap_CC)
	    haplotypicCC(tests,nch,true);
	  else if (par::test_hap_QTL)
	    haplotypicQTL(tests,nch,true);
	  else if (par::test_hap_TDT)
	    haplotypicTDT(tests,nch,true);
	}
      
      // Haplotype-specific test
      
      // tests[0] = 0
      // tests[1] = 1,2,...,h
      
      // tests[0] = 1
      // tests[1] = 0,2,3,...,h
      
      // etc
      
      for (int h=0; h < nh; h++)
	{
	  if (f[h] >= par::min_af)
	    {
	      tests.clear();
	      
	      for (int h2=0; h2 < nh; h2++)
		{
		  if (f[h2] >= par::min_af)
		    {
		      if (h==h2) 
			{
			  tests.insert(make_pair(h2,0));
			}
		      else tests.insert(make_pair(h2,1));
		    }
		}
	      
	      if (par::test_hap_CC)
		haplotypicCC(tests,2,true);
	      else if (par::test_hap_QTL)
		haplotypicQTL(tests,2,true);
	      else if (par::test_hap_TDT)
		haplotypicTDT(tests,2,true);
	    }
	}
      
    } // end of --hap-all routine
 
 

  // Dummy return vector (for now)
  return v;
  
  
}
  

vector_t HaploPhase::imputeGenotype(int i, int l)
{
  // Probability of AA, AB and BB for position 'l' 
  // (of ns SNPs) for individual 'i'

  vector_t g(3);
  
  if ( X || haploid ) 
    {
      g[0] = g[1] = g[2] = 0;
      return g;

      error("HaploPhase::imputeGenotypess() not yet set up for X \n");  
    }

  // Not able to be imputed?

  if (! include[i] )
    {
      g[0] = g[1] = g[2] = 0;
      return g;
    }
  
  
  // Unambiguous imputation?
  
  if (!ambig[i])
    {

      int h1 = hap1[i][0];
      int h2 = hap2[i][0];
      
      bool s1 = hap[h1][l];
      bool s2 = hap[h2][l];
      
      if ( s1 != s2 ) 
	g[1] = 1;
      else if ( s1 ) 
	g[0] = 1;
      else
	g[2] = 1;
      
      return g;
    }


  // Weighted, ambiguous imputation?

  for (int z=0; z<hap1[i].size(); z++)
    {
      
      // ?? include?? if (pp[i][max_z] >= par::hap_post_prob)
      
      int h1 = hap1[i][z];
      int h2 = hap2[i][z];
      
      bool s1 = hap[h1][l];
      bool s2 = hap[h2][l];
      
      if ( s1 != s2 ) 
	g[1] += pp[i][z];
      else if ( s1 ) 
	g[0] += pp[i][z];
      else
	g[2] += pp[i][z];
            
    } // next possible phase
  
  return g;
    
}  
  

double HaploPhase::imputeHaplotypes(int i, bool & n1, bool & n2)
{


  if ( X || haploid ) 
    error("HaploPhase::imputeHaplotypes() not yet set up for X \n");

  
  //////////////////////////////////////////////
  // Based on P(H|G) impute inferred haplotypes
  // for above-threshold individuals  
  
  // Not imputed
  double w = -1;
  
  // First for individuals of unambiguous phase

  if (include[i])
    {    
      if (!ambig[i])
	{
	  if (hap1[i][0] == test_hap) n1 = false;
	  else n1 = true;
	  
	  if (hap2[i][0] == test_hap) n2 = false;
	  else n2 = true;
	  
	  // Resolve potential het/missing coding confusion
	  
	  if ( n1 && (!n2) )
	    {
	      n1 = false;
	      n2 = true;	      
	    }
	  
	  // Unambiguous weighting (0, 1 or 2 copie of test_hap)
	  if (!n1)
	    {
	      if (!n2) w = 2;
	      else w = 1;
	    }
	  else w = 0;
	}
      else
	{
	  
	  // Second, for ambiguous individuals impute and assign weight
	  
	  int max_z = 0;
	  for (int z=0; z<hap1[i].size(); z++)
	    max_z = pp[i][z] > pp[i][max_z] ? z : max_z ;
	  
	  // Set missing by default
	  
	  n1 = true;
	  n2 = false;	      
	  
	  // Consider each phase z
	  
	  // Above threshold? 
	  if (pp[i][max_z] >= par::hap_post_prob)
	    {
	      
	      // Do we match 'test_hap' ( '1' allele ) 
	      // or not? ( '2' allele )
	      
	      if (hap1[i][max_z] == test_hap) n1 = false;
	      else n1 = true;
	  
	      if (hap2[i][max_z] == test_hap) n2 = false;
	      else n2 = true;
	  
	      // Resolve potential het/missing coding confusion
	  
	      if ( n1 && (!n2) )
		{
		  n1 = false;
		  n2 = true;	      
		}		  
	      
	  
	      // Unambiguous weighting (1 or 2 copies of test_hap)

	      // We are saying either 0, 1 or 2 copies
	      // Consider each haplotype

	      // Number imputed / Actual number
	      for (int z=0; z<hap1[i].size(); z++)
		{
		  if (hap1[i][z] == test_hap) w += pp[i][z];
		  if (hap2[i][z] == test_hap) w += pp[i][z];
		}
	  
	      w = pp[i][max_z] / w;

	      // 	  if (!n1)
	      //  	    {
	      //  	      if (!n2) w = pp[i][test_hap] * 1;
	      //  	      else w = pp[i][test_hap] * 0.5;
	      //  	    }
	      
	      // 	  int fac=1;
	      // 	  if ( hap1[i][max_z] != hap2[i][max_z] ) fac=2;
	      // 	  w = ( pp[i][max_z] - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] ) 
              //              / ( 1 - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] );
	      
	    }
	}     
    }
  else
    {
      // Excluded inviduals
      n1 = true;
      n2 = false;
    }
  
  
  return w;
}


void HaploPhase::resolveWithKids(int i)
{
  
  // Consider the founders in each family, who 
  // have at least 1 child, and a genotyped spouse

  // We require a full family, with two parents
  //  if ( ! f->parents ) return;

//   int pati = pat->

//   Individual * pat =  f->pat;
//   Individual * mat =  f->mat;

//   A/a B/b A/a B/b -> A/A B/B
//                      AB / AB 

//   for (int i=0; i< P.family[f].size(); i++)
//     cout << P.family[f]->fid << "\t"
// 	 << P.family[f]->iid << "\t"
// 	 << P.family[f]->pat->iid << "\t"
// 	 << P.family[f]->mat->iid << "\n";
         
}


void HaploPhase::phaseAndScoreNonfounder(int i)
{


  //////////////////////////////////////////////
  // Link this individual up with their parents

  int father = P.sample[i]->ip;
  int mother = P.sample[i]->im;
  
  bool nofather = false;
  bool nomother = false;

  if (father==-1) 
    nofather = true;
  else if (!include[father])
    nofather = true;

  if (mother==-1) 
    nomother = true;
  else if (!include[mother])
    nomother = true;
  
  // For TDT purposes, we require both parents to be 'observed'
  // i.e. so we never we to consider the "AllPhases" list (so
  //      we now do not bother generating it, i.e. enumerateAllPhases()
  //      function call is commented out in the main loop above

  if ( nofather || nomother ) 
  {
     include[i] = false;
     return;
  }


  int pat_phases = hap1[father].size();
  int mat_phases = hap1[mother].size();


  // Too much ambiguity?

  if ( pat_phases * mat_phases >= par::hap_max_nf_phases )
    {
      include[i] = false;
      return;
    }


  // Keep track of transmitted and non-transmitted 
  // haplotypes if performing a TDT-type analysis

  vector<vector<int> > trans1(0);
  vector<vector<int> > untrans1(0);
  

  /////////////////////////////////////////
  // Perform fill-in phasing for offspring
  
  // Step 1. Enumerate possible offspring phases
  //         or set to not include if too much missing

  enumeratePhase(i);
  

  //////////////////////////////////////////////
  // Do we want to attempt to reconstruct phase?
  
  if ( !include[i] ) 
    {     
      return;
    }



  ////////////////////////////////////////////////
  // Step 2. Joint distribution of parental phases
  
  vector<double> prob(0);
  double psum = 0;
  
  int pcnt=1;

  for (int z=0; z < hap1[i].size() ; z++)
    {
      
      // Is this phase possible given parental haplotypes?

      int h1 = hap1[i][z];
      int h2 = hap2[i][z];

      // Reset posterior probability for this phase
      pp[i][z] = 0;
      
      for (int z1=0; z1 < pat_phases ; z1++)
	for (int z2=0; z2 < mat_phases ; z2++)
	  {	
	    
	    // If no mother or father exists, we are 
	    // using the standard ph[] enumeration of
	    // all possible haplotypes
	    
	    int p1, p2, m1, m2;
	    
	    if (nofather)
	      {
		p1 = ph_hap1[z1];
		p2 = ph_hap2[z1];
	      }
	    else
	      {
		p1 = hap1[father][z1];
		p2 = hap2[father][z1];
	      }
	    
	    if (nomother)
	      {
		m1 = ph_hap1[z2];
		m2 = ph_hap2[z2];
	      }
	    else
	      {
		m1 = hap1[mother][z2];
		m2 = hap2[mother][z2];
	      }
	    

	    // Is this particular parental phase genotypically
	    // consistent with the offspring?  As haploid genotypes
	    // are coded as homozygotes, we do not need to worry about 
	    // this here

	    // We can add in MT / Y awareness here (i.e. such that 
	    // we can fill-in offspring haplotype with the appropriate
	    // parental haplotype, if the offspring has missing genotype
	    // data -- but, given the context, this is not a priority 
	    // and we will not worry about this for now -- i.e. only 
	    // consider autosomal case here -- the haploid cases will
	    // just get ignored.  But we should consider the X here

	    bool consistent = false;
	    
	    if ( X && P.sample[i]->sex ) 
	      { 
		// enumeratePhase() should only have specified
		// possible homozygous phases -- if male X chr, 
		// we only need to check that it is consistent 
		// with at least one maternal X
		
		if ( h1 == m1 || 
		     h1 == m2  ) consistent = true;
		
	      }
	    else if ( ( h1 == p1 && h2 == m1 ) || 
		      ( h1 == p1 && h2 == m2 ) || 
		      ( h1 == p2 && h2 == m1 ) || 
		      ( h1 == p2 && h2 == m2 ) || 
		      ( h1 == m1 && h2 == p1 ) || 
		      ( h1 == m1 && h2 == p2 ) || 
		      ( h1 == m2 && h2 == p1 ) || 
		      ( h1 == m2 && h2 == p2 ) ) consistent = true;
	    
	    if ( consistent ) 
	      {		  

		double p = 1;
		
		if (nofather)
		  p *= ph_freq[z1];
		else
		  if (ambig[father]) p *= pp[father][z1];
		
		if (nomother)
		  p *= ph_freq[z2];
		else
		  if (ambig[mother]) p *= pp[mother][z2];
		
// We explicitly consider both phases, we remove this line
//		if (h1!=h2) p *= 2;
		
		psum += p;
		

		// Track all consistent phases
		
		prob.push_back(p);
		
		if (ambig[i])
		  pp[i][z] += p;		
		

		///////////////////////////////
		// Now consider transmissions
		// (if we are performing a TDT)

		if (par::test_hap_TDT || par::proxy_TDT)
		  {
		    		    
		    vector<int> t1(nh);
		    vector<int> u1(nh);

		    // Father heterozygous? 
		    if ( p1 != p2 )
		      {
			// Mother homozygous?
			if ( m1 == m2 ) 
			  {
			    // then select a kid allele that matches, 
			    // and score the other one for transmission
			    if ( h1 == m1 ) 
			      {
				t1[h2]++;
				if (p1==h2)
				  u1[p2]++;
				else
				  u1[p1]++;
			      }
			    else
			      {
				t1[h1]++;
				if (p1==h1)
				  u1[p2]++;
				else
				  u1[p1]++;
				
			      }
			  }
			else
			  {
			    // Both parents are heterozygous, 
			    
			    // Transmitted alleles are unambiguous
			    t1[h1]++;
			    t1[h2]++;
			    
			    // Untransmitted alleles
			    // i.e. which two are left over 
			    // after accounting for the two
			    // transmitted alleles
			    
			    bool pat_accounted = false;
			    bool mat_accounted = false;
			    if ( p1 != h1 && p1 != h2 )
			      {
				u1[p1]++;
				pat_accounted = true;
			      }
			    else if ( p2 != h1 && p2 != h2 )
			      {
				u1[p2]++;
				pat_accounted = true;
			      }
			    
			    if ( m1 != h1 && m1 != h2 )
			      {
				u1[m1]++;
				mat_accounted = true;
			      }
			    else if  ( m2 != h1 && m2 != h2 )
			      {
				u1[m2]++;
				mat_accounted = true;
			      }
			    
			    // This only happens with AB x AB -> AB
			    if ( ! ( pat_accounted || mat_accounted ) ) 
			      {
				u1[h1]++;
				u1[h2]++;				    
			      }
			    else if ( ( ( !pat_accounted ) &&   mat_accounted   ) || 
				      (    pat_accounted   && (!mat_accounted)  ) )
			      {
				// If only 1 untransmitted allele accounted, for, it must
				// be the doubled allele that is untransmitted)
				
				if ( p1 == m1 || p1 == m2 ) 
				  u1[p1]++;
				else
				  u1[p2]++;
			      }
			    
			    
			    // AB AB    BB -- 2 acconuted for
			    //  AA
			    
			    // AB AB    AB -- 0 accounted for
			    //  AB 
			    
			    // AB AC    BA -- 1 accounted for
			    //  AC
			    
			    // AC AB    BA -- 1 accounted for
			    //  AC
			    
			    // AB CD    DB  -- 2 accounted for 
			    //  AC
			    
			  }
		      }
		    else if ( m1 != m2 ) 
		      {
			// Mother heterozygous, father homozygous
			if ( h1 == p1 )
			  {
			    t1[h2]++;
			    if (m1==h2)
			      u1[m2]++;
			    else
			      u1[m1]++;
			  }
			else
			  {
			    t1[h1]++;
			    if (m1==h1)
			      u1[m2]++;
			    else
			      u1[m1]++;
			  }
		      }
		    
		    
		    trans1.push_back(t1);
		    untrans1.push_back(u1);			
		    

		    if (par::verbose)
		      {
			Individual * pat = P.sample[father];
			Individual * mat = P.sample[mother];
			Individual * kid = P.sample[i];
			cout << "PHASE " << pcnt++ << "\t";
			cout << pat->fid << " : " 
			     << pat->iid << " x " 
			     << mat->iid << " = " 
			     << kid->iid << " : ";
			
			for (int i=0; i< t1.size(); i++)
			  cout << t1[i] << " ";
			cout << "\t";
			for (int i=0; i< u1.size(); i++)
			  cout << u1[i] << " ";
			cout << "\n";
		      }
		    
		    
		  } // End: score transmissions
	      } // End: if genotypically consistent
	  } // Next parental phase 
    } // Next offspring possible phase
  

  // Normalise probabilities

  for (int z=0; z < prob.size(); z++)
    prob[z] /= psum;
  
  if (ambig[i])
    for (int z=0; z < pp[i].size(); z++)
      pp[i][z] /= psum;
  
  

  // Finished counting: add to sample transmission counts
  
  if (par::test_hap_TDT || par::proxy_TDT)
    {      
      // For each haplotype
      for (int h=0; h<nh; h++)
	{
	  // For each possible phase for this non-founder
	  for (int z=0; z<prob.size(); z++)
	    {
	      trans[h] += trans1[z][h] * prob[z];
	      untrans[h] += untrans1[z][h] * prob[z];
	    }
	}
    }
  
  return; 
}


double HaploPhase::rsq_internal(int s1, int s2)
{
  // A convenience function for SNP x SNP r^2
  // i.e. here it does not matter which allele 
  // we consider, so just re-use mask

  if ( s1 > ns || s2 > ns ) 
    error("Problem in rsq_internal(int,int)");

  boolvec_t m1(ns,false);
  boolvec_t m2(ns,false);

  m1[s1] = true;
  m2[s2] = true;

  return rsq_internal(m1,m1,m2,m2);
}


double HaploPhase::rsq_internal(boolvec_t & mask1,
				boolvec_t & alleles1,
				boolvec_t & mask2,
				boolvec_t & alleles2)
{
  
  // Assume f[] has been populated with sensible values
  // and hap[][] contains alleles
  
  if ( mask1.size() != ns ||
       mask2.size() != ns ||
       alleles1.size() != ns ||
       alleles2.size() != ns )
    error("Internal error in Phase::rsq");
  
  //  ---X-X-  mask1
  //     0-0   alleles1
  
  //  ----X--  mask2
  //      1    alleles2
  
  // i.e. find r^2 between 00 haplotype made of SNPs 4 & 6 
  // from 7 SNP haplotype with allele 1 of SNP 5
  
  // Calculate frequency of first haplotype (fA)

  double fA = 0;
  double fB = 0;
  double fAB = 0, fAb = 0, faB = 0, fab = 0;

  for (int h = 0; h < nh; h++)
    {

      bool is_A = true;
      bool is_B = true;
      
      bool is_AB = true;
      bool is_Ab = true;
      bool is_aB = true;

      for (int s = 0; s < ns ; s++)
	{

	  if ( mask1[s] && hap[h][s] != alleles1[s] )
	    is_A = false;
	
	  if ( mask2[s] && hap[h][s] != alleles2[s] )
	    is_B = false;
	  
	  if ( ( mask1[s] && hap[h][s] != alleles1[s] ) ||
	       ( mask2[s] && hap[h][s] != alleles2[s] ) )
	    is_AB = false;
	
	  if ( ( mask1[s] && hap[h][s] != alleles1[s] ) ||
	       ( mask2[s] && hap[h][s] == alleles2[s] ) )
	    is_Ab = false;

	  if ( ( mask1[s] && hap[h][s] == alleles1[s] ) ||
	       ( mask2[s] && hap[h][s] != alleles2[s] ) )
	    is_aB = false;
	  
      }
      
      if ( is_A ) 
	fA += f[h];
      
      if ( is_B ) 
	fB += f[h];
      
      if ( is_AB ) 
	fAB += f[h];
      else if ( is_aB ) 
	faB += f[h];
      else if ( is_Ab ) 
	fAb += f[h];
      else
	fab += f[h];
      
      // Next haplotype
    } 

  double fa = 1 - fA;
  double fb = 1 - fB;

  double D = fAB - fA * fB;
  double denom = fA * fa * fB * fb;
  if ( denom == 0 )
    return -1;
  return (D*D) / denom;
  
}


double HaploPhase::freq(boolvec_t & mask1,
			boolvec_t & alleles1)
  
{
  
  // Assume f[] has been populated with sensible values
  // and hap[][] contains alleles
  
  if ( mask1.size() != ns ||
       alleles1.size() != ns )
    error("Internal error in Phase::rsq");
  
  //  ---X-X-  mask1
  //     0-0   alleles1
  
  double fA = 0;

  for (int h = 0; h < nh; h++)
    {

      bool is_A = true;

      for (int s = 0; s < ns ; s++)
	{
	  if ( mask1[s] && hap[h][s] != alleles1[s] )
	    is_A = false;
	}
      
      if ( is_A ) 
	fA += f[h];
      
      // Next haplotype
    } 

  return fA;
}


double HaploPhase::rsq(int l1, int l2)
{
  cout << " is HERE -2\n";
  reset();
  cout << " is HERE -2b\n";
  new_pred_locus.resize(1);
  new_map.resize(1);
  cout << " is HERE -1\n";

  vector<int> twoSNPs(2);
  twoSNPs[0] = l1;
  twoSNPs[1] = l2;

  new_pred_locus[0] = twoSNPs;
  new_map[0] = P.locus[l1];

  bool old_silent = par::silent;
  par::silent = true;
  cout << " is HERE 0\n";
  new_pred_allele = listPossibleHaplotypes(P, new_pred_locus[0]);
  cout << "about to phase\n";
  phaseAllHaplotypes();
  cout << "about to phase DONE\n";

  // hname = locus[l]->name;
  
  par::silent = old_silent;
  
  return rsq_internal(0,1);

}



void Plink::calcPairwiseLD()
{

  int l1 = getMarkerNumber(*this,par::ld_SNP1);
  int l2 = getMarkerNumber(*this,par::ld_SNP2);
  
  if ( l1 == l2 ) 
    error("Cannot compute LD with self");
  
  if ( l1 == -1 ) 
    error("--ld {marker} {marker}: first marker not found");

  if ( l2 == -1 ) 
    error("--ld {marker} {marker}: second marker not found");

  cout << "l1,l2 = " << l1 << " " << l2 << "\n";

printLOG("\nLD information for SNP pair [ " 
	   + par::ld_SNP1 + " " + par::ld_SNP2 + " ]\n");
  
  printLOG("\nr-sq = " + dbl2str( haplo->rsq(l1,l2) ) + "\n\n");

  return;
}
