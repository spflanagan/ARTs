#pragma once

//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//this contains the class population and the operations contained within that class. 
//most of the work for the ARTs program happens here

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include "random_numbers.h"
#include "classes.h"


using namespace std;

class population
{
public:
	int population_size, num_mal, num_fem, num_progeny, sex_trait;
	double sex_theta, sex_omega, courter_thresh, parent_thresh,pref_thresh;
	bool extinct;
	vector<double> theta, mean_fem_traits, mean_mal_traits, d_parentfreq, d_courterfreq;
	vector<individual> adults;
	vector<individual> progeny;
	vector<tracker> courter_qtls, parent_qtls, pref_qtls, courter_env_qtls, parent_env_qtls, pref_env_qtls, courter_thresh_qtls, parent_thresh_qtls, cthresh_env_qtls, pthresh_env_qtls;
	vector<tracker> maf, hs, maj_allele, ref_allele;


	population()
	{
		population_size = num_mal = num_fem = num_progeny = sex_trait = int();
		theta = mean_fem_traits = mean_mal_traits = vector<double>();
		sex_theta = sex_omega = courter_thresh = parent_thresh = pref_thresh = double();
		adults = progeny = vector<individual>();
		courter_env_qtls = courter_qtls = parent_env_qtls = parent_qtls = pref_env_qtls = pref_qtls = parent_thresh_qtls = courter_thresh_qtls =cthresh_env_qtls = pthresh_env_qtls = maf = hs = vector<tracker>();
		extinct = false;
		sgenrand(time(0));
	}
	
	//initialize
	void determine_pref_thresh(parameters gp)
	{
		if (gp.court_trait)
			pref_thresh = courter_thresh;
		else
			pref_thresh = parent_thresh;
	}
	void set_thresholds(parameters gp)
	{//this function is run after traits have been set
		int j, adult_count;
		vector<double> tempallele1, tempallele2;
		
		courter_thresh = parent_thresh = adult_count = 0;
		//if traits are conditional
		if (gp.courter_conditional || gp.parent_conditional)
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive)
				{
					adults[j].assign_conditional_traits(gp);
					if (gp.courter_conditional)
						courter_thresh = courter_thresh + adults[j].courter_trait;
					if (gp.parent_conditional)
						parent_thresh = parent_thresh + adults[j].parent_trait;
					adult_count++;
				}
			}
			courter_thresh = courter_thresh / adult_count;
			parent_thresh = parent_thresh / adult_count;
		}
		//if traits have a genetic basis
		double courter_sd, parent_sd;
		if (gp.court_trait)
		{
			courter_thresh = calc_mean_courter_ae(gp);
			courter_sd = calc_courter_sd(gp, courter_thresh);
		}
		if (gp.parent_trait)
		{
			parent_thresh = calc_mean_parent_ae(gp);
			parent_sd = calc_parent_sd(gp, parent_thresh);
		}
		//if thresholds evolve then we need alleles
		if (gp.thresholds_evolve)
		{
			for (j = 0; j < gp.num_alleles; j++) //set up specific alleles.
			{
				if(gp.courter_conditional || gp.court_trait)
					tempallele1.push_back(randnorm((courter_thresh/gp.num_qtl), courter_sd/gp.num_qtl));
				if(gp.parent_conditional || gp.parent_trait)
					tempallele2.push_back(randnorm((parent_thresh/gp.num_qtl), parent_sd/gp.num_qtl));
			}
		}
		//assign morph to each ind
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (gp.thresholds_evolve)
			{
				if (gp.courter_conditional || gp.court_trait)
					adults[j].assign_threshold_gt(gp, j, tempallele1, true);
				if (gp.parent_conditional || gp.parent_trait)
					adults[j].assign_threshold_gt(gp, j, tempallele2, false);
			}
			adults[j].determine_threshold(gp, courter_thresh, parent_thresh);
			if (!adults[j].female)
			{
				if (gp.courter_conditional || gp.court_trait)
				{
					if (adults[j].courter_trait < adults[j].ind_cthresh)
						adults[j].courter = false;
					else
						adults[j].courter = true;
				}
				if (gp.parent_conditional || gp.parent_trait)
				{
					if (adults[j].parent_trait < adults[j].ind_pthresh)
						adults[j].parent = false;
					else
						adults[j].parent = true;
				}
			}
		}
		determine_pref_thresh(gp);//initialize the preference threshold
	}
	void initialize_supergene(parameters & gp)
	{
		int m, mm, sg_start;
		double supergene_nmarkers = gp.supergene_prop*gp.num_markers;
		if (gp.qtl_per_chrom[0] == (gp.num_qtl / gp.num_chrom))//otherwise the supergene info was specified by user
		{
			//chose a chromosome for the supergene
			int sg_chrom = randnum(gp.num_chrom);
			for (m = 0; m < gp.num_chrom; m++)
			{
				if (sg_chrom == m)
					gp.qtl_per_chrom[m] = gp.num_qtl;
				else
					gp.qtl_per_chrom[m] = 0;
			}
		}
		
		if (gp.court_trait && gp.parent_trait)
		{
			if (gp.thresholds_in_supergene)
			{
				supergene_nmarkers = 4 * supergene_nmarkers;
				if (supergene_nmarkers <= 4 * gp.num_qtl)//make sure there's enough space in the supergene region.
					supergene_nmarkers = 4 * gp.num_qtl;
				sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
				for (m = 0; m < gp.num_chrom; m++)
				{
					courter_qtls.push_back(tracker());
					parent_qtls.push_back(tracker());
					courter_thresh_qtls.push_back(tracker());
					parent_thresh_qtls.push_back(tracker());
					for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
					{
						courter_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						parent_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						courter_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						parent_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					}
				}
			}
			else
			{
				supergene_nmarkers = 2 * supergene_nmarkers;
				if (supergene_nmarkers <= 2 * gp.num_qtl)//make sure there's enough space in the supergene region.
					supergene_nmarkers = 2 * gp.num_qtl;
				sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
				for (m = 0; m < gp.num_chrom; m++)
				{
					courter_qtls.push_back(tracker());
					parent_qtls.push_back(tracker());
					for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
					{
						courter_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						parent_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					}
				}
			}
		}
		else
		{
			if ((gp.court_trait || gp.parent_trait) && gp.thresholds_in_supergene)
			{
				if (supergene_nmarkers <= 2*gp.num_qtl)//make sure there's enough space in the supergene region.
					supergene_nmarkers = 2*gp.num_qtl;
				sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
			}
			else
			{
				if (supergene_nmarkers <= gp.num_qtl)//make sure there's enough space in the supergene region.
					supergene_nmarkers = gp.num_qtl;
				sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
			}
			for (m = 0; m < gp.num_chrom; m++)
			{
				if (gp.court_trait)
				{
					if (gp.thresholds_in_supergene)
						courter_thresh_qtls.push_back(tracker());
					courter_qtls.push_back(tracker());
				}
				if (gp.parent_trait)
				{
					parent_qtls.push_back(tracker());
					if (gp.thresholds_in_supergene)
						parent_thresh_qtls.push_back(tracker());
				}
				for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
				{
					if (gp.court_trait)
					{
						courter_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						if(gp.thresholds_in_supergene)
							courter_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					}
					if (gp.parent_trait)
					{
						parent_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
						if(gp.thresholds_in_supergene)
							parent_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					}
				}
			}
		}
		if (gp.courter_conditional && gp.parent_conditional && gp.thresholds_in_supergene)
		{
			supergene_nmarkers = 2 * supergene_nmarkers;
			if (supergene_nmarkers <= 2 * gp.num_qtl)//make sure there's enough space in the supergene region.
				supergene_nmarkers = 2 * gp.num_qtl;
			sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
			for (m = 0; m < gp.num_chrom; m++)
			{
				courter_thresh_qtls.push_back(tracker());
				parent_thresh_qtls.push_back(tracker());
				for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
				{
					courter_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					parent_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
				}
			}
		}
		if (gp.courter_conditional && gp.thresholds_in_supergene && !gp.parent_conditional)
		{
			if (supergene_nmarkers <= gp.num_qtl)//make sure there's enough space in the supergene region.
				supergene_nmarkers = gp.num_qtl;
			sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
			for (m = 0; m < gp.num_chrom; m++)
			{
				courter_thresh_qtls.push_back(tracker());
				for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
				{
					courter_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
				}
			}
		}
		if (!gp.courter_conditional && gp.parent_conditional && gp.thresholds_in_supergene)
		{
			if (supergene_nmarkers <= gp.num_qtl)//make sure there's enough space in the supergene region.
				supergene_nmarkers = gp.num_qtl;
			sg_start = randnum(gp.num_markers - supergene_nmarkers); //choose marker start
			for (m = 0; m < gp.num_chrom; m++)
			{
				parent_thresh_qtls.push_back(tracker());
				for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
				{
					parent_thresh_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
				}
			}
		}
}
	void initialize_adults(parameters gp)
	{
		int j, jj, jjj;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			adults.push_back(individual());
			if (gp.court_trait || gp.ind_pref || gp.cor_prefs || !gp.random_mating)
				adults[j].courter_trait = 0;
			if (gp.parent_trait)
				adults[j].parent_trait = 0;
			adults[j].pot_rs = gp.max_fecund;
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				adults[j].maternal.push_back(chromosome());
				adults[j].paternal.push_back(chromosome());
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci.push_back(int());
					adults[j].paternal[jj].loci.push_back(int());
				}
				if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].courter_ae.push_back(double());
						adults[j].paternal[jj].courter_ae.push_back(double());
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].parent_ae.push_back(double());
						adults[j].paternal[jj].parent_ae.push_back(double());
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].pref_ae.push_back(double());
						adults[j].paternal[jj].pref_ae.push_back(double());
					}
				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
					{
						for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
						{
							adults[j].maternal[jj].courter_thresh.push_back(double());
							adults[j].paternal[jj].courter_thresh.push_back(double());
						}
					}
					if (gp.parent_conditional || gp.parent_trait)
					{
						for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
						{
							adults[j].maternal[jj].parent_thresh.push_back(double());
							adults[j].paternal[jj].parent_thresh.push_back(double());
						}
					}
				}
			}
			if (gp.gene_network)//set up interactions
			{
				if (gp.court_trait == true)
				{
					adults[j].initialize_network(gp, adults[j].courter_int, adults[j].courter_Y, adults[j].courter_x, adults[j].courter_Z);
				}
				if (gp.parent_trait == true)
				{
					adults[j].initialize_network(gp, adults[j].parent_int, adults[j].parent_Y, adults[j].parent_x, adults[j].parent_Z);
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					adults[j].initialize_network(gp, adults[j].pref_int, adults[j].pref_Y, adults[j].pref_x, adults[j].pref_Z);
				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
						adults[j].initialize_network(gp, adults[j].cthresh_int, adults[j].cthresh_Y, adults[j].cthresh_x, adults[j].cthresh_Z);
					if (gp.parent_conditional || gp.parent_trait)
						adults[j].initialize_network(gp, adults[j].pthresh_int, adults[j].pthresh_Y, adults[j].pthresh_x, adults[j].pthresh_Z);
				}
			}
		}
	}
	void initialize_progeny(parameters gp)
	{
		int j, jj, jjj;
		//determine the maximum possible fecundity
		int max_potential_fecund = max(gp.max_fecund, gp.rs_c);
		max_potential_fecund = max(max_potential_fecund, gp.rs_nc);
		max_potential_fecund = max(max_potential_fecund, gp.rs_np);
		max_potential_fecund = max(max_potential_fecund, gp.rs_p);
		for (j = 0; j < (max_potential_fecund*gp.carrying_capacity); j++)
		{
			progeny.push_back(individual());
			if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
				progeny[j].courter_trait = 0;
			if (gp.parent_trait)
				progeny[j].parent_trait = 0;
			progeny[j].pot_rs = gp.max_fecund;
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				progeny[j].maternal.push_back(chromosome());
				progeny[j].paternal.push_back(chromosome());
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					progeny[j].maternal[jj].loci.push_back(int());
					progeny[j].paternal[jj].loci.push_back(int());
				}
				if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						progeny[j].maternal[jj].courter_ae.push_back(double());
						progeny[j].paternal[jj].courter_ae.push_back(double());
						if (gp.thresholds_evolve)
						{
							progeny[j].maternal[jj].courter_thresh.push_back(double());
							progeny[j].paternal[jj].courter_thresh.push_back(double());
						}
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						progeny[j].maternal[jj].parent_ae.push_back(double());
						progeny[j].paternal[jj].parent_ae.push_back(double());
						if (gp.thresholds_evolve)
						{
							progeny[j].maternal[jj].parent_thresh.push_back(double());
							progeny[j].paternal[jj].parent_thresh.push_back(double());
						}
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						progeny[j].maternal[jj].pref_ae.push_back(double());
						progeny[j].paternal[jj].pref_ae.push_back(double());
					}
				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
					{
						for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
						{
							progeny[j].maternal[jj].courter_thresh.push_back(double());
							progeny[j].paternal[jj].courter_thresh.push_back(double());
						}
					}
					if (gp.parent_conditional || gp.parent_trait)
					{
						for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
						{
							progeny[j].maternal[jj].parent_thresh.push_back(double());
							progeny[j].paternal[jj].parent_thresh.push_back(double());
						}
					}
				}
			}
			if (gp.gene_network)//set up interactions
			{
				if (gp.court_trait == true)
				{
					progeny[j].initialize_network(gp, progeny[j].courter_int, progeny[j].courter_Y, progeny[j].courter_x, progeny[j].courter_Z);
				}
				if (gp.parent_trait == true)
				{
					progeny[j].initialize_network(gp, progeny[j].parent_int, progeny[j].parent_Y, progeny[j].parent_x, progeny[j].parent_Z);
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					progeny[j].initialize_network(gp, progeny[j].pref_int, progeny[j].pref_Y, progeny[j].pref_x, progeny[j].pref_Z);
				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
						progeny[j].initialize_network(gp, progeny[j].cthresh_int, progeny[j].cthresh_Y, progeny[j].cthresh_x, progeny[j].cthresh_Z);
					if (gp.parent_conditional || gp.parent_trait)
						progeny[j].initialize_network(gp, progeny[j].pthresh_int, progeny[j].pthresh_Y, progeny[j].pthresh_x, progeny[j].pthresh_Z);
				}
			}
		}
	}
	void assign_env_qtls(vector<tracker>&env_qtls, parameters gp)
	{
		int j, jj, jjj, count;
		//assign qtls to be environmental.
		for (j = 0; j < gp.num_chrom; j++)
		{
			env_qtls.push_back(tracker());
			for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
			{
				env_qtls[j].per_locus.push_back(-1);
			}
		}
		count = 0;
		while(count < gp.num_env_qtl)
		{
			int index = 0;
			int this_qtl = randnum(gp.num_qtl);
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
				{
					index++;
					if (index == this_qtl && env_qtls[jj].per_locus[jjj] < 0) 
					{
						env_qtls[jj].per_locus[jjj] = this_qtl;
						count++;
					}
				}
			}
		}
	}
	void assign_env_qtls(vector<tracker>&env_qtls, parameters gp, int global_qtl_per_chrom)//for female prefs with supergene
	{
		int j, jj, jjj;
		//assign qtls to be environmental.
		for (j = 0; j < gp.num_chrom; j++)
		{
			env_qtls.push_back(tracker());
			for (jj = 0; jj < global_qtl_per_chrom; jj++)
			{
				env_qtls[j].per_locus.push_back(-1);
			}
		}
		for (j = 0; j < gp.num_env_qtl; j++)
		{
			int index = 0;
			int this_qtl = randnum(gp.num_qtl);
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				for (jjj = 0; jjj < global_qtl_per_chrom; jjj++)
				{
					index++;
					if (index == this_qtl)
						env_qtls[jj].per_locus[jjj] = this_qtl;
				}
			}
		}
	}
	void initialize(parameters & gp)
	{
		std::cout << "Initializing a population.\n";
		population_size = gp.carrying_capacity;
		int j, jj, jjj;
		num_fem = num_mal = 0;
		sex_theta = -100;
		//set up the qtls
		if (gp.supergene)
		{//if it's a supergene, male traits QTLs are together
			initialize_supergene(gp);
			//set up females
			//female preference is not a part of a supergene if it's independent
			if (gp.cor_prefs == true)
			{//then it's the same as the courter trait
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						pref_qtls[j].per_locus.push_back(courter_qtls[j].per_locus[jj]);
				}
			}
			if (gp.ind_pref == true)
			{//otherwise it's its own thing.
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < (gp.num_qtl / gp.num_chrom); jj++)
						pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
			}
			if ((gp.cor_prefs || gp.ind_pref) && pref_qtls.size() == 0)//sanity check!
			{
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < (gp.num_qtl / gp.num_chrom); jj++)
						pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
				gp.ind_pref = true;
			}
		}
		else
		{//otherwise they're distributed along the chromosomes.
			if (gp.court_trait == true)
			{
				for (j = 0; j < gp.num_chrom; j++)
				{
					courter_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						courter_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
			}
			if (gp.parent_trait == true)
			{
				for (j = 0; j < gp.num_chrom; j++)
				{
					parent_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						parent_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
			}
			
			//set up females
			if (gp.cor_prefs == true)
			{//then it's the same as the courter trait
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						pref_qtls[j].per_locus.push_back(courter_qtls[j].per_locus[jj]);
				}
			}
			if (gp.ind_pref == true)
			{//otherwise it's its own thing.
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
			}
			if ((gp.cor_prefs || gp.ind_pref) && pref_qtls.size() == 0)//sanity check!
			{
				for (j = 0; j < gp.num_chrom; j++)
				{
					pref_qtls.push_back(tracker());
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
				gp.ind_pref = true;
			}
		}
		if (gp.thresholds_evolve && !gp.thresholds_in_supergene)//put this here in case supergene is called but thresholds aren't there.
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				if (gp.parent_trait || gp.parent_conditional)
					parent_thresh_qtls.push_back(tracker());
				if (gp.court_trait || gp.courter_conditional)
					courter_thresh_qtls.push_back(tracker());
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (gp.parent_trait || gp.parent_conditional)
						parent_thresh_qtls[j].per_locus.push_back(randnum(gp.num_markers));
					if (gp.court_trait || gp.courter_conditional)
						courter_thresh_qtls[j].per_locus.push_back(randnum(gp.num_markers));
				}
			}
		}
		//Set up environmental effect QTL locations if necessary
		if (gp.gene_network)
		{
			//randomly choose gp.num_env_qtl from the set
			if (gp.court_trait)
				assign_env_qtls(courter_env_qtls, gp);
			if (gp.parent_trait)
				assign_env_qtls(parent_env_qtls, gp);
			if (gp.cor_prefs || gp.ind_pref)
				assign_env_qtls(pref_env_qtls, gp, gp.num_qtl/gp.num_chrom);
			if (gp.thresholds_evolve)
			{
				if (gp.court_trait || gp.courter_conditional)
					assign_env_qtls(cthresh_env_qtls, gp);
				if (gp.parent_conditional || gp.parent_trait)
					assign_env_qtls(pthresh_env_qtls, gp);
			}
		}

		//set up trackers
		for (j = 0; j < gp.num_chrom; j++)
		{
			maf.push_back(tracker());
			ref_allele.push_back(tracker());
			hs.push_back(tracker());
			for (jj = 0; jj < gp.num_markers; jj++)
			{
				ref_allele[j].per_locus.push_back(double());
				maf[j].per_locus.push_back(double());
				hs[j].per_locus.push_back(double());
			}
		}
		//set up adults
		initialize_adults(gp);//just set up the vectors
		//set up progeny
		initialize_progeny(gp);//just set up the vectors
		//assign loci and alleles
		vector<double> tempallele1, tempallele2, tempallele3;
		for (j = 0; j < gp.num_alleles; j++) //set up specific alleles.
		{
			tempallele1.push_back(randnorm(0, gp.allelic_std_dev/gp.num_qtl));
			tempallele2.push_back(randnorm(0, gp.allelic_std_dev / gp.num_qtl));
			tempallele3.push_back(randnorm(0, gp.allelic_std_dev / gp.num_qtl));
		}
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (j < gp.carrying_capacity)
				adults[j].alive = true;
			else
				adults[j].alive = false;//it's an empty slot for migrants
			adults[j].lifetime_rs = 0;
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				if (gp.court_trait)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
						adults[j].paternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].parent_ae[jjj] = tempallele2[j%gp.num_alleles];
						adults[j].paternal[jj].parent_ae[jjj] = tempallele2[j%gp.num_alleles];
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					for (jjj = 0; jjj < adults[j].maternal[jj].pref_ae.size(); jjj++)
					{
						adults[j].maternal[jj].pref_ae[jjj] = tempallele3[j%gp.num_alleles];
						adults[j].paternal[jj].pref_ae[jjj] = tempallele3[j%gp.num_alleles];
					}
				}
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci[jjj] = j%gp.num_alleles;
					adults[j].paternal[jj].loci[jjj] = j%gp.num_alleles;
				}
			}
			if (gp.court_trait) 
			{
				if (gp.gene_network)
					adults[j].calc_courter_trait(gp, courter_env_qtls, 0);//to initialize no cue but yes stability check
				else
					adults[j].calc_courter_trait(gp); 
			}
			if (gp.parent_trait)
			{
				if (gp.gene_network)
					adults[j].calc_parent_trait(gp, parent_env_qtls, 0);
				else
					adults[j].calc_parent_trait(gp);
			}
			if (genrand() < 0.5)
			{
				adults[j].female = true;
				adults[j].pot_rs = gp.max_fecund;
				num_fem++;
			}
			else
			{
				adults[j].female = false;
				num_mal++;
			}
		}//carrying_capacity
		//set up morphs
		set_thresholds(gp); //this also takes care of conditional (random) traits
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (gp.gene_network)
				adults[j].update_traits(gp, courter_thresh, parent_thresh, pref_thresh, 0, courter_env_qtls, parent_env_qtls, pref_env_qtls, pthresh_env_qtls, cthresh_env_qtls);
			else
				adults[j].update_traits(gp, courter_thresh, parent_thresh, pref_thresh);
		}
		population_size = gp.carrying_capacity;
		
	}
	bool sanity_checks(parameters& gp)
	{
		int p, pp;
		bool run_program, each_ok;
		run_program = true;
		if (gp.supergene)
		{
			each_ok = true;
			if (gp.no_genetics)
				each_ok = false;
			for (p = 0; p < gp.num_chrom; p++)
			{
				if (gp.qtl_per_chrom[p] != 0 && gp.qtl_per_chrom[p] != gp.num_qtl)
					each_ok = false;
			}
			if (gp.supergene_prop < (gp.num_qtl / gp.num_markers))
				each_ok = false;
			if (!gp.court_trait && !gp.parent_trait)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Supergene parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.cor_prefs)
		{
			each_ok = true;
			if (!gp.court_trait)
				each_ok = false;
			if (gp.ind_pref)
				each_ok = false;
			if (gp.random_mating)//just reset it
				gp.random_mating = false;
			if (pref_qtls.size() == 0)
				each_ok = false;
			for (p = 0; p < gp.num_chrom; p++)
			{
				for (pp = 0; pp < gp.qtl_per_chrom[p]; pp++)
				{
					if (pref_qtls[p].per_locus[pp] != courter_qtls[p].per_locus[pp])
						each_ok = false;
				}
			}
			if (!each_ok)
			{
				std::cout << "\nERROR: Correlated female preference parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.ind_pref)
		{
			each_ok = true;
			if (gp.cor_prefs)
				each_ok = false;
			if (gp.random_mating)
				gp.random_mating = false;//just reset it here
			if (pref_qtls.size() == 0)
				each_ok = false;
			if (!gp.court_trait && !gp.parent_trait && !gp.courter_conditional && !gp.parent_conditional)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Independent female preference parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.court_trait)
		{
			each_ok = true;
			if (courter_qtls.size() == 0)
				each_ok = false;
			if (courter_thresh == 0)
				each_ok = false;
			if (gp.rs_c < 0 || gp.rs_nc < 0)
				each_ok = false;
			if (gp.courter_conditional)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Courtship trait parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.parent_trait)
		{
			each_ok = true;
			if (parent_qtls.size() == 0)
				each_ok = false;
			if (parent_thresh == 0)
				each_ok = false;
			if (gp.rs_p < 0 || gp.rs_np < 0)
				each_ok = false;
			if (gp.parent_conditional)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Parent trait parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.gene_network)
		{
			each_ok = true;
			if (gp.no_genetics)
				each_ok = false;
			if (!gp.court_trait && !gp.parent_trait && !gp.ind_pref && !gp.cor_prefs)
				each_ok = false;
			if (gp.num_env_qtl <= 0 || gp.num_env_qtl >= gp.num_qtl)
				each_ok = false;
			if (gp.court_trait)
			{
				if (courter_env_qtls.size() == 0)
					each_ok = false;
				if (courter_qtls.size() == 0)
					each_ok = false;
			}
			if (gp.parent_trait)
			{
				if (parent_env_qtls.size() == 0)
					each_ok = false;
				if (parent_qtls.size() == 0)
					each_ok = false;
			}
			if (gp.ind_pref || gp.cor_prefs)
			{
				if (pref_env_qtls.size() == 0)
					each_ok = false;
				if (pref_qtls.size() == 0)
					each_ok = false;
			}
			if (gp.thresholds_evolve)
			{
				if (gp.courter_conditional || gp.court_trait)
				{
					if (cthresh_env_qtls.size() == 0)
						each_ok = false;
					if (courter_thresh_qtls.size() == 0)
						each_ok = false;
				}
				if (gp.parent_conditional || gp.parent_trait)
				{
					if (pthresh_env_qtls.size() == 0)
						each_ok = false;
					if (parent_thresh_qtls.size() == 0)
						each_ok = false;
				}
			}
			if (!each_ok)
			{
				std::cout << "\nERROR: Gene network parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.env_cue)
		{
			each_ok = true;
			if (!gp.gene_network)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Environmental cue parameters are insane.";
				run_program = each_ok;
			}
		}
		if (!gp.random_mating)
		{//if there's not random mating then there must at least be a courtship trait
			each_ok = true;
			if (!gp.court_trait && !gp.parent_trait && !gp.courter_conditional && !gp.parent_conditional)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Random mating parameters are insane.";
				run_program = each_ok;
			}
		}
		if (gp.FD_pref)
		{//frequency dependent mate choice; does not rely on female preference trait
			each_ok = true;
			if (gp.random_mating) //reset it here
				gp.random_mating = false;
			if (!gp.court_trait && !gp.parent_trait && !gp.courter_conditional && !gp.parent_conditional)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Frequency dependent mate choice is insane.";
				run_program = each_ok;
			}
		}
		if (gp.CD_pref)
		{//condition dependent mate choice, depdent on female's preference trait/fecundity
			each_ok = true;
			if (gp.random_mating) //reset it here
				gp.random_mating = false; 
			if (!gp.court_trait && !gp.parent_trait && !gp.courter_conditional && !gp.parent_conditional)
				each_ok = false;
			if (!gp.ind_pref && !gp.cor_prefs)
				each_ok = false;
			if (pref_qtls.size() == 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Condition dependent mate choice is insane.";
				run_program = each_ok;
			}
		}
		if (gp.FD_court)
		{
			each_ok = true;
			if (!gp.court_trait && !gp.courter_conditional)
				each_ok = false;
			if (courter_qtls.size() == 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Frequency dependent selection on courtship trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.FD_parent)
		{
			each_ok = true;
			if (!gp.parent_trait && !gp.parent_conditional)
				each_ok = false;
			if (parent_qtls.size() == 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Frequency dependent selection on parent trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.CD_court)
		{
			each_ok = true;
			if (!gp.polygyny)
				gp.polygyny = true;
			if (!gp.court_trait && !gp.courter_conditional)
				each_ok = false;
			if (courter_qtls.size() == 0)
				each_ok = false;
			if (gp.cond_adj <= 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Condition dependent selection on courtship trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.CD_parent)
		{
			each_ok = true;
			if (!gp.polygyny)
				gp.polygyny = true; 
			if (!gp.parent_trait && !gp.parent_conditional)
				each_ok = false;
			if (parent_qtls.size() == 0)
				each_ok = false;
			if (gp.cond_adj <= 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Condition dependent selection on parent trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.courter_conditional)
		{
			each_ok = true;
			if (courter_qtls.size() > 0 || courter_env_qtls.size() > 0)
				each_ok = false;
			if (gp.court_trait)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Conditional courter trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.parent_conditional)
		{
			each_ok = true;
			if (parent_qtls.size() > 0 || parent_env_qtls.size() > 0)
				each_ok = false;
			if (gp.parent_trait)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Conditional parent trait is insane.";
				run_program = each_ok;
			}
		}
		if (gp.thresholds_evolve)
		{
			each_ok = true;
			if (courter_thresh_qtls.size() == 0 && parent_thresh_qtls.size() == 0)
				each_ok = false;
			if (!each_ok)
			{
				std::cout << "\nERROR: Evolving thresholds are insane.";
				run_program = each_ok;
			}
		}
		return run_program;
	}
	
	//pop size
	void determine_pop_size(parameters gp)
	{
		int j;
		population_size = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive)
				population_size++;
		}
		if (population_size == 0)
			std::cout << "\nPopulation has crashed!";
	}
	void determine_sex_nums(parameters gp)
	{
		int j;
		num_mal = num_fem = 0;
		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].alive)
			{
				if (adults[j].female)
					num_fem++;
				else
					num_mal++;
			}
		}
		if (num_mal == 0)
			std::cout << "\nNo males in the population!";
		if (num_fem == 0)
			std::cout << "\nNo females in the population!";
	}

	//set up male traits
	double calc_mean_courter_ae(parameters gp)
	{
		int j,count;
		double mean;
		mean = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				if (gp.gene_network)
					adults[j].calc_courter_trait(gp, courter_env_qtls, 0);
				else
					adults[j].calc_courter_trait(gp);
				mean = mean + adults[j].courter_trait;
				count++;
			}
		}
		mean = mean / count;		
		return(mean);
	}
	double calc_courter_sd(parameters gp, double courter_mean)
	{//assumes you've already calculated the mean
		int j, count;
		double sd, var;
		sd = var = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				var = var + (adults[j].courter_trait - courter_mean)*(adults[j].courter_trait - courter_mean);
				count++;
			}
		}
		var = var / count;
		sd = sqrt(var);
		return sd;
	}
	double calc_mean_parent_ae(parameters gp)
	{
		int j, count;
		double mean;
		mean = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				if (gp.gene_network)
					adults[j].calc_parent_trait(gp, parent_env_qtls, 0);
				else
					adults[j].calc_parent_trait(gp);
				mean = mean + adults[j].parent_trait;
				count++;
			}
		}
		mean = mean / count;
		return(mean);
	}
	double calc_parent_sd(parameters gp, double parent_mean)
	{//assumes you've already calculated the mean
		int j, count;
		double sd, var;
		sd = var = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				var = var + (adults[j].parent_trait - parent_mean)*(adults[j].parent_trait - parent_mean);
				count++;
			}
		}
		var = var / count;
		sd = sqrt(var);
		return sd;
	}
	double calc_freq_courter(parameters gp)
	{
		int j;
		double freqT, count;
		freqT = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				if (adults[j].courter)
					freqT++;
				count++;
			}
		}
		freqT = freqT / count;
		return freqT;
	}
	double calc_freq_parent(parameters gp)
	{
		int j;
		double freqT, count;
		freqT = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				if (adults[j].parent)
					freqT++;
				count++;
			}
		}
		freqT = freqT / count;
		return freqT;
	}
	double calc_freq_pref(parameters gp)
	{
		int j;
		double freqT, count;
		freqT = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				if (adults[j].female_pref==1)
					freqT++;
				count++;
			}
		}
		freqT = freqT / count;
		return freqT;
	}
	vector<double> calc_freq_morphs(parameters gp)
	{
		int j,male_count;
		vector<double> morph_freqs;
		male_count = 0;
		if ((gp.parent_trait || gp.parent_conditional) && (gp.court_trait || gp.courter_conditional))//sanity check
		{
			for (j = 0; j < 4; j++) //initialize
				morph_freqs.push_back(0);
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && !adults[j].female)
				{
					if (!adults[j].courter && !adults[j].parent)
						morph_freqs[0]++;
					if (adults[j].courter && !adults[j].parent)
						morph_freqs[1]++;
					if (!adults[j].courter && adults[j].parent)
						morph_freqs[2]++;
					if (adults[j].courter && adults[j].parent)
						morph_freqs[3]++;
					male_count++;
				}
			}
			for (j = 0; j < 4; j++)
				morph_freqs[j] = morph_freqs[j] / male_count;
		}
		return morph_freqs;
	}
	void parent_fd_rs(parameters gp)
	{
		int j;
		double freqT;
		freqT = calc_freq_parent(gp);
		if (freqT < 0.5) //parents have higher RS
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && !adults[j].female)
				{
					if(adults[j].parent)
						adults[j].pot_rs = gp.rs_p;
					else
						adults[j].pot_rs = gp.rs_np;
				}
			}
		}
		else //non parents have higher RS
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && !adults[j].female)
				{
					if (adults[j].parent)
						adults[j].pot_rs = gp.rs_np;
					else
						adults[j].pot_rs = gp.rs_p;
				}
			}
		}
	}

	//all of the ways to calculate preferences
	void pref_better_morph(parameters gp)//prefer the court/parental trait
	{
		int j;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
				adults[j].female_pref = 1;
		}
	}
	void pref_femtrait(parameters gp)
	{//female's own trait is measured against the preferred threshold
		for (int j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				if (gp.gene_network)
					adults[j].calc_preference_trait(gp, pref_thresh, pref_env_qtls, 0);//what is the environmental cue??
				else
					adults[j].calc_preference_trait(gp, pref_thresh);
			}
		}
	}
	void pref_fd_courter(parameters gp)
	{
		int j;
		double freqT;
		freqT= 0;
		freqT = calc_freq_courter(gp);
		if (freqT < 0.5) // females prefer courter (courter is less frequent)
		{
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = 1;
			}
		}
		else //females prefer non-courter (courter is more frequent)
		{
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = 0;
			}
		}
	}
	void pref_fd_parent(parameters gp)
	{
		int j;
		double freqT;
		freqT= 0;
		freqT = calc_freq_parent(gp);
		if (freqT < 0.5) // parent is less frequent; prefer parent
		{
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = 1;
			}
		}
		else //parent is more frequent; prefer non-parent
		{
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = 0;
			}
		}
	}
	void pref_cd_courter(parameters gp, double meanT, double meanF, double xswitch)
	{
		int j;
		double mean;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				if(gp.gene_network)
					adults[j].calc_preference_trait(gp, pref_thresh,pref_env_qtls,0);//consider environmental cue...
				else
					adults[j].calc_preference_trait(gp, pref_thresh);
				mean = poissonrand(gp.max_fecund + adults[j].female_pref);
				adults[j].pot_rs = adults[j].female_pref + mean; //this is her fecundity
				if (xswitch < adults[j].pot_rs)//then she prefers rs_c
					adults[j].female_pref = meanT;
				else
					adults[j].female_pref = meanF;
			}
		}		
	}
	void assign_preference(parameters gp)
	{
		if (!gp.random_mating) //only waste the time if you're gonna use a preference
		{
			if (gp.ind_pref || gp.cor_prefs)//if female traits have a genetic architecture
			{
				if (!gp.FD_pref && !gp.CD_pref)//just plain old random mating
				{								//using the female's genetic architecture of preference to determine which morph
					pref_femtrait(gp);
				}
				if (gp.CD_pref)//then need to calculate female fecundity and compare it to S(rs_p - rs_np).
				{
					double xswitch, meanT, meanF;
					meanT = 1;
					meanF = 0;//based on morphs
					if (gp.court_trait)
					{
						xswitch = (num_mal / gp.max_encounters) / (gp.rs_c - gp.rs_nc);
					}
					else
					{
						xswitch = (num_mal / gp.max_encounters) / (gp.rs_p - gp.rs_np);
					}
					pref_cd_courter(gp, meanT, meanF, xswitch);
				}
				if (gp.FD_pref)//then need the frequency of morphs and female prefers less frequent one.
				{
					if (gp.court_trait || gp.courter_conditional)
						pref_fd_courter(gp);
					else //it's based on the parent trait
						pref_fd_parent(gp);

				}
			}
			else
			{
				if (gp.FD_pref)//then need the frequency of morphs and female prefers less frequent one.
				{
					if (gp.court_trait || gp.courter_conditional)
						pref_fd_courter(gp);
					else //it's based on the parent trait
						pref_fd_parent(gp);
				}
				else
					pref_better_morph(gp);//females just across the board prefer the better morph
			}								//regardless of genetic architecture
		}
	}

	//calculating cues
	void nesting_cue_eval(parameters gp, int mal_id,bool courter_trait, double fem_opt)
	{
		double c;
		if (gp.gene_network)//sanity check
		{
			if (courter_trait == true)
			{
				if (gp.env_cue)
					c = double(adults[mal_id].mate_found / gp.max_num_mates) - exp((adults[mal_id].courter - fem_opt)*(adults[mal_id].courter - fem_opt));
				else
					c = 0;
				adults[mal_id].calc_courter_trait(gp, courter_env_qtls, c,false);
				if (gp.thresholds_evolve)
					adults[mal_id].determine_threshold(gp, courter_thresh, parent_thresh, cthresh_env_qtls, pthresh_env_qtls, c, false);
				adults[mal_id].assign_court_morph(gp);
			}
			else
			{
				if (gp.env_cue)
					c = double(adults[mal_id].mate_found / gp.max_num_mates) - exp((adults[mal_id].parent - fem_opt)*(adults[mal_id].parent - fem_opt));
				else
					c = 0;
				adults[mal_id].calc_parent_trait(gp, parent_env_qtls, c, false);
				if (gp.thresholds_evolve)
					adults[mal_id].determine_threshold(gp, courter_thresh, parent_thresh, cthresh_env_qtls, pthresh_env_qtls, c, false);
				adults[mal_id].assign_parent_morph(gp);
			}
		}
	}
	void preference_cue_eval(parameters gp, int fem_id, int num_encounters, int male_trait)
	{
		double c;
		if (gp.gene_network)
		{
			if (gp.env_cue)
				c = double(num_encounters) / gp.max_encounters - exp((adults[fem_id].female_pref - male_trait)*(adults[fem_id].female_pref - male_trait));
			else
				c = 0;
			adults[fem_id].calc_preference_trait(gp, pref_thresh, pref_env_qtls, c, false);
			//I could add frequency dependence but I'm not sure how...
		}
	}
	void parent_cue_eval(parameters gp, int mal_id)
	{
		double c;
		if (gp.gene_network)
		{
			if (gp.parent_trait)
			{
				if (gp.env_cue)
					c = double(adults[mal_id].mate_found) / gp.max_num_mates - exp((adults[mal_id].pot_rs - gp.rs_p)*(adults[mal_id].pot_rs - gp.rs_p));
				else
					c = 0;
				adults[mal_id].calc_parent_trait(gp, parent_env_qtls, c, false);
				if (gp.thresholds_evolve)
					adults[mal_id].determine_threshold(gp, courter_thresh, parent_thresh, cthresh_env_qtls, pthresh_env_qtls, c, false);
				adults[mal_id].assign_parent_morph(gp);
			}
		}
	}
	//functions for mating
	void recombine_network(parameters gp, int progeny_index, int mom_index, int dad_index)
	{
		int jj, jjj;
		for (jj = 0; jj < gp.num_qtl; jj++)
		{
			for (jjj = 0; jjj < gp.num_qtl; jjj++)
			{
				if (gp.court_trait)
				{
					if (genrand() < 0.5)
					{
						progeny[progeny_index].courter_Y[jj].per_locus[jjj] = adults[mom_index].courter_Y[jj].per_locus[jjj];
						progeny[progeny_index].courter_int[jj].per_locus[jjj] = adults[mom_index].courter_int[jj].per_locus[jjj];
					}
					else
					{
						progeny[progeny_index].courter_Y[jj].per_locus[jjj] = adults[dad_index].courter_Y[jj].per_locus[jjj];
						progeny[progeny_index].courter_int[jj].per_locus[jjj] = adults[dad_index].courter_int[jj].per_locus[jjj];
					}
				}
				if (gp.parent_trait)
				{
					if (genrand() < 0.5)
					{
						progeny[progeny_index].parent_Y[jj].per_locus[jjj] = adults[mom_index].parent_Y[jj].per_locus[jjj];
						progeny[progeny_index].parent_int[jj].per_locus[jjj] = adults[mom_index].parent_int[jj].per_locus[jjj];
					}
					else
					{
						progeny[progeny_index].parent_Y[jj].per_locus[jjj] = adults[dad_index].parent_Y[jj].per_locus[jjj];
						progeny[progeny_index].parent_int[jj].per_locus[jjj] = adults[dad_index].parent_int[jj].per_locus[jjj];
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					if (genrand() < 0.5)
					{
						progeny[progeny_index].pref_Y[jj].per_locus[jjj] = adults[mom_index].pref_Y[jj].per_locus[jjj];
						progeny[progeny_index].pref_int[jj].per_locus[jjj] = adults[mom_index].pref_int[jj].per_locus[jjj];
					}
					else
					{
						progeny[progeny_index].pref_Y[jj].per_locus[jjj] = adults[dad_index].pref_Y[jj].per_locus[jjj];
						progeny[progeny_index].pref_int[jj].per_locus[jjj] = adults[dad_index].pref_int[jj].per_locus[jjj];
					}
				}
				if (gp.thresholds_evolve)
				{
					if(gp.court_trait || gp.courter_conditional)
					{
						if (genrand() < 0.5)
						{
							progeny[progeny_index].cthresh_Y[jj].per_locus[jjj] = adults[mom_index].cthresh_Y[jj].per_locus[jjj];
							progeny[progeny_index].cthresh_int[jj].per_locus[jjj] = adults[mom_index].cthresh_int[jj].per_locus[jjj];
						}
						else
						{
							progeny[progeny_index].cthresh_Y[jj].per_locus[jjj] = adults[dad_index].cthresh_Y[jj].per_locus[jjj];
							progeny[progeny_index].cthresh_int[jj].per_locus[jjj] = adults[dad_index].cthresh_int[jj].per_locus[jjj];
						}
					}
					if (gp.parent_trait || gp.parent_conditional)
					{
						if (genrand() < 0.5)
						{
							progeny[progeny_index].pthresh_Y[jj].per_locus[jjj] = adults[mom_index].pthresh_Y[jj].per_locus[jjj];
							progeny[progeny_index].pthresh_int[jj].per_locus[jjj] = adults[mom_index].pthresh_int[jj].per_locus[jjj];
						}
						else
						{
							progeny[progeny_index].pthresh_Y[jj].per_locus[jjj] = adults[dad_index].pthresh_Y[jj].per_locus[jjj];
							progeny[progeny_index].pthresh_int[jj].per_locus[jjj] = adults[dad_index].pthresh_int[jj].per_locus[jjj];
						}
					}
				}
			}
		}
	}
	void recombine_chromosome(chromosome& chrom, individual &parent, double expected_recomb_events, int which_chrom, parameters gp, bool alive, int prog_id)//tracker& court_qtls, tracker& parent_qtls, tracker& pref_qtls, 
	{
		int RCi, RCj;
		int NumberRecombEvents = 0;
		int segment_start[22], segment_end[20];
		int break_point[20];
		vector<bool> recombine;
		if (expected_recomb_events == 0)
			NumberRecombEvents = 0;
		if (expected_recomb_events < 6)
			NumberRecombEvents = poissonrand(expected_recomb_events);
		if (expected_recomb_events > 6)
			NumberRecombEvents = positiveroundnorm(expected_recomb_events, sqrt(expected_recomb_events));
		for (RCi = 0; RCi < 20; RCi++)
			break_point[RCi] = gp.num_markers + 1;
		bool segment_maternal[22];
		int num_segments;
		bool start_maternal;
		if (NumberRecombEvents > 20)
			NumberRecombEvents = 20;
		//actually recombine
		if (NumberRecombEvents > 0)
		{
			for (RCi = 0; RCi < NumberRecombEvents; RCi++) 
				break_point[RCi] = randnum(gp.num_markers);

			//sort breakpoints
			sort(begin(break_point), end(break_point));
			//first segment maternal or paternal?
			if (genrand() < 0.5)
				start_maternal = true;
			else
				start_maternal = false;

			num_segments = 1;
			segment_start[0] = 0;
			segment_maternal[0] = start_maternal;
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
			{
				segment_end[RCi] = break_point[RCi];
				segment_start[RCi + 1] = break_point[RCi];
				if (gp.supergene)
				{
					if (gp.court_trait && courter_qtls[which_chrom].per_locus.size() > 0)
					{
						double first_qtl, last_qtl;
						first_qtl = courter_qtls[which_chrom].first_qtl();
						last_qtl = courter_qtls[which_chrom].last_qtl();
						if (first_qtl >= segment_start[RCi] && last_qtl >= segment_end[RCi])
						{
							alive = false;
						}
					}
					if (gp.parent_trait && parent_qtls[which_chrom].per_locus.size() > 0)
					{
						int first_qtl, last_qtl;
						first_qtl = parent_qtls[which_chrom].first_qtl();
						last_qtl = parent_qtls[which_chrom].last_qtl();
						if(first_qtl >= segment_start[RCi] && last_qtl <= segment_end[RCi])
						{
							alive = false;
						}
					}
				}

				if (segment_maternal[RCi])
					segment_maternal[RCi + 1] = false;
				else
					segment_maternal[RCi + 1] = true;
				num_segments++;
			}
			segment_end[RCi] = gp.num_markers;
			//pass allelic info to recombined chromosome
			for (RCi = 0; RCi < num_segments; RCi++)
			{
				if(alive == true)
				{
					if (segment_maternal[RCi])
					{
						for (RCj = segment_start[RCi]; RCj < segment_end[RCi]; RCj++)
						{
							chrom.loci[RCj] = parent.maternal[which_chrom].loci[RCj];
						}
					}
					else
					{
						for (RCj = segment_start[RCi]; RCj < segment_end[RCi]; RCj++)
						{
							chrom.loci[RCj] = parent.paternal[which_chrom].loci[RCj];
						}
					}
					for (RCj = 0; RCj < gp.qtl_per_chrom[which_chrom]; RCj++)
					{
						if (gp.court_trait)
						{
							if (courter_qtls[which_chrom].per_locus[RCj] >= segment_start[RCi] && courter_qtls[which_chrom].per_locus[RCj] >= segment_end[RCi])
							{
								if (segment_maternal[RCi])
									chrom.courter_ae[RCj] = parent.maternal[which_chrom].courter_ae[RCj];
								else
									chrom.courter_ae[RCj] = parent.paternal[which_chrom].courter_ae[RCj];
							}
						}
						if (gp.parent_trait)
						{
							if (parent_qtls[which_chrom].per_locus[RCj] >= segment_start[RCi] && parent_qtls[which_chrom].per_locus[RCj] >= segment_end[RCi])
							{
								if (segment_maternal[RCi])
									chrom.parent_ae[RCj] = parent.maternal[which_chrom].parent_ae[RCj];
								else
									chrom.parent_ae[RCj] = parent.paternal[which_chrom].parent_ae[RCj];
								
							}
							
						}
						if (gp.thresholds_evolve)
						{
							if (gp.courter_conditional || gp.court_trait)
							{
								if (courter_thresh_qtls[which_chrom].per_locus[RCj] >= segment_start[RCi] && courter_thresh_qtls[which_chrom].per_locus[RCj] >= segment_end[RCi])
								{
									if (segment_maternal[RCi])
										chrom.courter_thresh[RCj] = parent.maternal[which_chrom].courter_thresh[RCj];
									else
										chrom.courter_thresh[RCj] = parent.paternal[which_chrom].courter_thresh[RCj];
								}
							}
							if (gp.parent_conditional || gp.parent_trait)
							{
								if (parent_thresh_qtls[which_chrom].per_locus[RCj] >= segment_start[RCi] && parent_thresh_qtls[which_chrom].per_locus[RCj] >= segment_end[RCi])
								{
									if (segment_maternal[RCi])
										chrom.parent_thresh[RCj] = parent.maternal[which_chrom].parent_thresh[RCj];
									else
										chrom.parent_thresh[RCj] = parent.paternal[which_chrom].parent_thresh[RCj];
								}
							}
						}
						if (gp.ind_pref || gp.cor_prefs)
						{
							for (RCj = 0; RCj < parent.maternal[which_chrom].pref_ae.size(); RCj++)
							{
								if (pref_qtls[which_chrom].per_locus[RCj] >= segment_start[RCi] && pref_qtls[which_chrom].per_locus[RCj] >= segment_end[RCi])
								{
									if (segment_maternal[RCi])
										chrom.pref_ae[RCj] = parent.maternal[which_chrom].pref_ae[RCj];
									else
										chrom.pref_ae[RCj] = parent.paternal[which_chrom].pref_ae[RCj];
								}
							}
						}	
					}
				}//alive
			}//num segments			
		}//numberrecombevents > 0
		else
		{
			//no recombination
			if (genrand() < 0.5)
			{//maternal
				for (RCi = 0; RCi < gp.num_markers; RCi++)
					chrom.loci[RCi] = parent.maternal[which_chrom].loci[RCi];
				for (RCi = 0; RCi < gp.qtl_per_chrom[which_chrom]; RCi++)
				{
					if (gp.court_trait)
					{
						chrom.courter_ae[RCi] = parent.maternal[which_chrom].courter_ae[RCi];
						chrom.courter_ae[RCi] = parent.maternal[which_chrom].courter_ae[RCi];
					}
					if (gp.parent_trait)
					{
						chrom.parent_ae[RCi] = parent.maternal[which_chrom].parent_ae[RCi];
						chrom.parent_ae[RCi] = parent.maternal[which_chrom].parent_ae[RCi];
					}
					if (gp.thresholds_evolve)
					{
						if (gp.courter_conditional || gp.court_trait)
						{
							chrom.courter_thresh[RCi] = parent.maternal[which_chrom].courter_thresh[RCi];
							chrom.courter_thresh[RCi] = parent.maternal[which_chrom].courter_thresh[RCi];
						}
						if (gp.parent_conditional || gp.parent_trait)
						{
							chrom.parent_thresh[RCi] = parent.maternal[which_chrom].parent_thresh[RCi];
							chrom.parent_thresh[RCi] = parent.maternal[which_chrom].parent_thresh[RCi];
						}
					}
				}
				
				if (gp.ind_pref || gp.cor_prefs)
				{
					for (RCi = 0; RCi < parent.maternal[which_chrom].pref_ae.size(); RCi++)
					{
						chrom.pref_ae[RCi] = parent.maternal[which_chrom].pref_ae[RCi];
						chrom.pref_ae[RCi] = parent.maternal[which_chrom].pref_ae[RCi];
					}
				}
			}
			else
			{//paternal
				for (RCi = 0; RCi < gp.num_markers; RCi++)
					chrom.loci[RCi] = parent.paternal[which_chrom].loci[RCi];
				for (RCi = 0; RCi < gp.qtl_per_chrom[which_chrom]; RCi++)
				{
					if (gp.court_trait)
					{
						chrom.courter_ae[RCi] = parent.paternal[which_chrom].courter_ae[RCi];
						chrom.courter_ae[RCi] = parent.paternal[which_chrom].courter_ae[RCi];	
					}
					if (gp.parent_trait)
					{
						chrom.parent_ae[RCi] = parent.paternal[which_chrom].parent_ae[RCi];
						chrom.parent_ae[RCi] = parent.paternal[which_chrom].parent_ae[RCi];
					}	
					if (gp.thresholds_evolve)
					{
						if (gp.courter_conditional || gp.court_trait)
						{
							chrom.courter_thresh[RCi] = parent.paternal[which_chrom].courter_thresh[RCi];
							chrom.courter_thresh[RCi] = parent.paternal[which_chrom].courter_thresh[RCi];
						}
						if (gp.parent_conditional || gp.parent_trait)
						{
							chrom.parent_thresh[RCi] = parent.paternal[which_chrom].parent_thresh[RCi];
							chrom.parent_thresh[RCi] = parent.paternal[which_chrom].parent_thresh[RCi];
						}
					}
				}
				if (gp.ind_pref || gp.cor_prefs)
				{
					for (RCi = 0; RCi < parent.paternal[which_chrom].pref_ae.size(); RCi++)
					{
						chrom.pref_ae[RCi] = parent.paternal[which_chrom].pref_ae[RCi];
						chrom.pref_ae[RCi] = parent.paternal[which_chrom].pref_ae[RCi];
					}
				}
			}
		}
	}
	
	void making_babies(parameters gp, int fecundity, int& num_progeny, int mom_index, int dad_index)
	{
		int jj, jjj;

		for (jj = 0; jj < fecundity; jj++)
		{
			if (num_progeny >= (gp.max_fecund*gp.carrying_capacity))
				num_progeny = (gp.max_fecund*gp.carrying_capacity) - 1;

			progeny[num_progeny].alive = true;
			progeny[num_progeny].mom = mom_index;
			progeny[num_progeny].dad = dad_index;
			if (genrand() < 0.5)
				progeny[num_progeny].female = true;
			else
				progeny[num_progeny].female = false;
			
			for (jjj = 0; jjj < gp.num_chrom; jjj++)
			{
				if (progeny[num_progeny].alive)
				{
					recombine_chromosome(progeny[num_progeny].maternal[jjj], adults[mom_index], gp.recombination_rate, jjj, gp, progeny[num_progeny].alive, num_progeny);
					recombine_chromosome(progeny[num_progeny].paternal[jjj], adults[dad_index], gp.recombination_rate, jjj, gp, progeny[num_progeny].alive, num_progeny);
				}
			}
			if (progeny[num_progeny].alive)
			{
				//also need to do mutation
				if (gp.gene_network)//need to evaluate if it's compatible
				{
					recombine_network(gp, num_progeny, mom_index, dad_index);
					progeny[num_progeny].mutation(gp, courter_qtls, parent_qtls, pref_qtls, courter_thresh_qtls, parent_thresh_qtls, courter_env_qtls, parent_env_qtls, cthresh_env_qtls, pthresh_env_qtls, pref_env_qtls);
					if (gp.courter_conditional || gp.parent_conditional)
						progeny[num_progeny].assign_conditional_traits(gp);//do this before updating traits, but after mutation
					progeny[num_progeny].update_traits(gp, courter_thresh, parent_thresh, pref_thresh,
						0, courter_env_qtls, parent_env_qtls, pref_env_qtls, pthresh_env_qtls,cthresh_env_qtls);
				}
				else
				{
					progeny[num_progeny].mutation(gp, courter_qtls, parent_qtls, pref_qtls,courter_thresh_qtls,parent_thresh_qtls,courter_env_qtls,parent_env_qtls,cthresh_env_qtls,pthresh_env_qtls,pref_env_qtls);
					if (gp.courter_conditional || gp.parent_conditional)
						progeny[num_progeny].assign_conditional_traits(gp);//do this before updating traits, but after mutation
					progeny[num_progeny].update_traits(gp, courter_thresh, parent_thresh, pref_thresh);
				}
				if (progeny[num_progeny].alive)
				{										
					num_progeny++;
				}
			}
		}
	}
	void per_female_mating(parameters gp, vector<int> & male_index)
	{
		int j, first, first_progeny, fem_ms;
		bool nest_found;
		//mating happens
		fem_ms = 0;
		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].female && adults[j].alive) //loop through the females.
			{
				nest_found = choose_nest(j, male_index, gp);
				if (nest_found)
				{
					if (gp.gene_network)
						parent_cue_eval(gp, adults[j].mate_id);
					first_progeny = fertilization(j, gp);
					nest_survival(gp, first_progeny, adults[j].mate_id);
					fem_ms++;
				}
			}
		}
	}
	void nest_and_fertilize(parameters gp, bool output, string out_name)
	{
		int j, first_progeny;
		bool nest_found;
		vector<int> male_index;
		ofstream out_file;
		int fem_ms;

		fem_ms = num_progeny = num_fem = num_mal = 0;

		//set up pop level metrics
		determine_pop_size(gp);
		assign_preference(gp);
		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].female && adults[j].alive)
				num_fem++;
			else
			{
				if (adults[j].alive)
				{
					num_mal++;
					male_index.push_back(j);
					if (gp.court_trait || gp.courter_conditional)//make sure morphs are assigned.
						adults[j].assign_court_morph(gp);
					if (gp.parent_trait || gp.parent_conditional)
						adults[j].assign_parent_morph(gp);
					//if parent traits don't exist then set them all to false so all males have equal chance.
					if (!gp.parent_conditional && !gp.parent_trait)
						adults[j].parent = false;
				}
			}
			adults[j].mate_found = 0;
		}
		if(gp.verbose)
			std::cout << ", " << num_mal << " males, " << num_fem << " females" << flush;
		if(gp.per_fem_mating)
			per_female_mating(gp, male_index);
		else
		{
			for (j = 0; j < adults.size(); j++)
			{
				if (adults[j].female && adults[j].alive) //loop through the females.
				{
					nest_found = choose_nest(j, male_index, gp);
					if (nest_found)
					{
						if (gp.gene_network)
							parent_cue_eval(gp, adults[j].mate_id);
						first_progeny = fertilization(j, gp);
						nest_survival(gp, first_progeny, adults[j].mate_id);
						fem_ms++;
					}
				}
			}
		}
		if (gp.verbose)
			std::cout <<", and " << fem_ms << " mated" << std::flush;
	}
	bool random_mating(parameters gp,int fem_index, vector<int> & male_index)
	{
		int male_id, encounters, irndnum;
		bool mate_found = false;
		encounters = 0;
		while (!mate_found && (encounters < gp.max_encounters))//female mates once
		{
			irndnum = randnum(num_mal);
			male_id = male_index[irndnum];
			if (adults[male_id].alive)
			{
				if (gp.polygyny || adults[male_id].mate_found == 0)//either polygyny is ok or if monogamy the male hasn't mated yet
				{
					mate_found = true;
					adults[fem_index].mate_found++;
					adults[male_id].mate_found++;
					adults[fem_index].mate_id = male_id;
				}
			}
			encounters++;
		}
		if (mate_found)
		{
			if (gp.CD_court)
			{
				if (adults[male_id].courter)//it decreases
					adults[male_id].courter_trait = adults[male_id].courter_trait - gp.cond_adj;
				else//it increases
					adults[male_id].courter_trait = adults[male_id].courter_trait + gp.cond_adj;
				adults[male_id].assign_court_morph(gp);
			}
			if (gp.CD_parent)
			{
				if (adults[male_id].parent)//it decreases
					adults[male_id].parent_trait = adults[male_id].parent_trait - gp.cond_adj;
				else//it increases
					adults[male_id].parent_trait = adults[male_id].parent_trait + gp.cond_adj;
				adults[male_id].assign_parent_morph(gp);
			}
		}
		return mate_found;
	}
	bool choose_nest(int fem_index, vector<int> & male_index, parameters gp)
	{//females choose one male to give her eggs to.
		int male_id, encounters, irndnum;
		bool mate_found = false;
		encounters = 0;
		if (gp.random_mating)//then they randomly find males
		{
			mate_found = random_mating(gp, fem_index, male_index);
		}//end random mating
		else //preference for one morph or the other 
		{ 
			vector <int> acceptable_males;
			while (encounters < gp.max_encounters)
			{
				irndnum = randnum(num_mal);
				male_id = male_index[irndnum];
				if (adults[male_id].alive)
				{
					if (adults[male_id].mate_found < gp.max_num_mates)
					{
						if (gp.gene_network)
						{//redetermine traits and thresholds, given the mating experience and optimality of trait
							nesting_cue_eval(gp, male_id, gp.court_trait, adults[fem_index].female_pref);
							if ((gp.ind_pref || gp.cor_prefs) && gp.court_trait)
								preference_cue_eval(gp, fem_index, encounters, adults[male_id].courter);
							if ((gp.ind_pref || gp.cor_prefs) && !gp.court_trait)
								preference_cue_eval(gp, fem_index, encounters, adults[male_id].parent);
						}
						if (gp.court_trait || gp.courter_conditional)
						{
							if (adults[fem_index].female_pref == adults[male_id].courter)	
								acceptable_males.push_back(male_id);
						}
						else//parent_trait
						{
							if (adults[fem_index].female_pref == adults[male_id].parent)
								acceptable_males.push_back(male_id);
						}
					}
				}
				encounters++;
			}
			if (acceptable_males.size() >= 1)
			{
				int randbest = randnum(acceptable_males.size());//if there are multiple males, she chooses one randomly
				male_id = acceptable_males[randbest];
				mate_found = true;
				adults[fem_index].mate_found++;
				adults[male_id].mate_found++;
				adults[fem_index].mate_id = male_id;
				if (gp.CD_court)
				{
					if (adults[male_id].courter)//it decreases
						adults[male_id].courter_trait = adults[male_id].courter_trait - gp.cond_adj;
					else//it increases
						adults[male_id].courter_trait = adults[male_id].courter_trait + gp.cond_adj;
					adults[male_id].assign_court_morph(gp);
				}
				if (gp.CD_parent)
				{
					if (adults[male_id].parent)//it decreases
						adults[male_id].parent_trait = adults[male_id].parent_trait - gp.cond_adj;
					else//it increases
						adults[male_id].parent_trait = adults[male_id].parent_trait + gp.cond_adj;
					adults[male_id].assign_parent_morph(gp);
				}
			}
			else
			{
				mate_found = random_mating(gp, fem_index, male_index);
			}
		}//end of finding the mates
		return mate_found;
	}
	int fertilization(int fem_id,parameters gp)
	{
		int k,kk,fecundity, randmale, max_sperm, first_progeny, num_mates, num_prog;
		int male_id = adults[fem_id].mate_id;
		vector<int> male_ids;
		vector<double> fecundity_share;	
		vector <int> sneakers;
		first_progeny = num_progeny;//set up this tracker
		
		max_sperm = adults[male_id].pot_rs;
		
		//identify sneakers/eligible bachelors
		for (k = 0; k < adults.size(); k++)
		{
			if (adults[k].alive && !adults[k].female && !adults[k].parent)
			{
				if(adults[k].pot_rs > 0)
					sneakers.push_back(k);
			}
		}
		if (sneakers.size() == 0)//if all males are parents, no one will sneak
		{
			if (adults[male_id].pot_rs > 0)//make sure the male has the sperm to fertilize the offspring
			{
				if (adults[male_id].pot_rs >= adults[fem_id].pot_rs)//they can only produce as many offspring as the male has sperm for
					fecundity = adults[fem_id].pot_rs;
				else
					fecundity = adults[male_id].pot_rs;
				making_babies(gp, fecundity, num_progeny, fem_id, male_id);
				adults[male_id].pot_rs = adults[male_id].pot_rs - fecundity; //reduce male's RS based on how many babies he's already made.
			}
		}
		else
		{
			//get the non-parental male mates
			male_ids.push_back(male_id);
			fecundity_share.push_back(0);
			if (sneakers.size() < gp.max_num_mates)
				num_mates = sneakers.size() + 1;
			else
				num_mates = gp.max_num_mates;
			int tries = 0;
			while (male_ids.size() < num_mates && tries < 10000)
			{
				randmale = randnum(sneakers.size());
				if (randmale != male_id)
				{
					if (!adults[sneakers[randmale]].parent && (adults[sneakers[randmale]].pot_rs > 0))//if it's a sneaker: sanity check!
					{
						fecundity_share.push_back(0);
						male_ids.push_back(sneakers[randmale]);
						int male_sperm = (gp.sperm_comp_r*adults[sneakers[randmale]].pot_rs);
						if (male_sperm > adults[sneakers[randmale]].pot_rs)
							male_sperm = adults[sneakers[randmale]].pot_rs;
						max_sperm = max_sperm + male_sperm;//if sperm_comp_r is 1, they're all equally weighted
						//adults[sneakers[randmale]].mate_found++;//keep track of individual reproductive success
					}
				}
				tries++;//just to break out of a while loop in case.
			}
			fecundity_share[0] = double(adults[male_id].pot_rs) / double(max_sperm); //it doesn't get weighted if r <= 1
			for (k = 1; k < fecundity_share.size(); k++)
				fecundity_share[k] = double(adults[male_ids[k]].pot_rs*gp.sperm_comp_r) / double(max_sperm);

			//now generate the offspring
			num_prog = 0;
			for (k = 0; k < fecundity_share.size(); k++)
			{
				fecundity = fecundity_share[k] * adults[fem_id].pot_rs;
				if ((num_prog + fecundity) > adults[fem_id].pot_rs)
					fecundity = adults[fem_id].pot_rs - num_prog;
				if (fecundity > adults[male_ids[k]].pot_rs)
					fecundity = min(adults[male_ids[k]].pot_rs, (adults[fem_id].pot_rs - num_prog));
				if (fecundity > 0)//make sure you can make this baby
				{
					making_babies(gp, fecundity, num_progeny, fem_id, male_ids[k]);
					adults[male_ids[k]].pot_rs = adults[male_ids[k]].pot_rs - fecundity; //reduce male's RS based on how many babies he's already made.
					num_prog = num_prog + fecundity;
				}
				if ((num_prog > gp.max_fecund) && gp.verbose)
					std::cout << "\n\tWARNING:" << num_prog << " were produced.";
			}
		}
		return first_progeny;
	}
	void nest_survival(parameters gp, int prog_start, int mal_id)
	{
		int k;
		double surv_rand;
		if (gp.parent_conditional || gp.parent_trait)//assign survival based on parent trait
		{											// otherwise all of them survive (don't need to do anything)
			for (k = prog_start; k < num_progeny; k++)
			{
				if (progeny[k].alive)//only deal with ones that are alive at this point
				{
					surv_rand = genrand();
					if (adults[mal_id].parent)
					{
						if (surv_rand < gp.egg_surv_parent)
							progeny[k].alive = true;
						else
							progeny[k].alive = false;
					}
					else
					{
						if (surv_rand < gp.egg_surv_noparent)
							progeny[k].alive = true;
						else
							progeny[k].alive = false;
					}
				}
			}
		}
	}

	//selection
	double via_against_courter(double selection_strength, int ind_id)
	{
		double surv_prob,optimum;
		optimum = 0;
		surv_prob = exp(-1 * (progeny[ind_id].courter - optimum)*(progeny[ind_id].courter - optimum)
			/ (2 * selection_strength));
		return surv_prob;
	}
	double via_against_parent(double selection_strength, int ind_id)
	{
		double surv_prob, optimum;
		optimum = 0;
		surv_prob = exp(-1 * (progeny[ind_id].parent - optimum)*(progeny[ind_id].parent - optimum)
			/ (2 * selection_strength));
		return surv_prob;
	}
	double via_fd_courter(parameters gp, int ind_id)
	{
		double surv_prob, optimum, courter_freq;
		courter_freq = calc_freq_courter(gp);
		if (courter_freq > 0.5)
			optimum = 0;
		else
			optimum = 1;
		surv_prob = exp(-1 * (progeny[ind_id].courter - optimum)*(progeny[ind_id].courter - optimum)
			/ (2 * gp.via_sel_strength));
		return surv_prob;
	}
	double via_fd_parent(parameters gp, int ind_id)
	{
		double surv_prob, optimum, parent_freq;
		parent_freq = calc_freq_parent(gp);
		if (parent_freq > 0.5)
			optimum = 1;
		else
			optimum = 0;
		surv_prob = exp(-1 * (progeny[ind_id].parent - optimum)*(progeny[ind_id].parent - optimum)
			/ (2 * gp.via_sel_strength));
		return surv_prob;
	}
	void viability_selection(parameters gp)
	{
		int j, ProgAlive;
		double dSurvProb;
		double drnum1;
		int malecount = 0;
		ProgAlive = 0;

		for (j = 0; j < num_progeny; j++)
		{
			if (progeny[j].alive)
			{
				if(gp.gene_network)
					progeny[num_progeny].update_traits(gp, courter_thresh, parent_thresh, pref_thresh,
						0, courter_env_qtls, parent_env_qtls, pref_env_qtls, pthresh_env_qtls, cthresh_env_qtls);
				else
					progeny[j].update_traits(gp, courter_thresh, parent_thresh, pref_thresh);
				if (progeny[j].female)
				{
					progeny[j].alive = true;
					ProgAlive++;
				}
				else//selection only on males
				{
					if (gp.via_sel_strength > 0)
					{
						dSurvProb = 1;
						if ((gp.court_trait || gp.courter_conditional) && !gp.FD_court)
							dSurvProb = dSurvProb*via_against_courter(gp.via_sel_strength, j);
						if (gp.FD_court)
							dSurvProb = dSurvProb*via_fd_courter(gp, j);
						if ((gp.parent_trait || gp.parent_conditional) && !gp.FD_court)
							dSurvProb = dSurvProb*via_against_parent(gp.via_sel_strength, j);
						if (gp.FD_parent)
							dSurvProb = dSurvProb*via_fd_parent(gp, j);
					}
					else
						dSurvProb = 1;
					//std::cout<<dSurvProb<<'\n';
					drnum1 = genrand();
					if (drnum1 < dSurvProb)
					{
						progeny[j].alive = true;
						ProgAlive++;
						malecount++;
					}
					else
						progeny[j].alive = false;
				}
			} // end of j
		}
		if (gp.verbose)
			std::cout << ", " << ProgAlive << " progeny" << std::flush;
	}

	//stochastic survival
	void pass_on_env_qtls(parameters gp, int adult_index, int progeny_index)
	{
		int pp, ppp;
		//pass on env effects
		for (pp = 0; pp < gp.num_qtl; pp++)
		{
			if (pp >= gp.num_env_qtl)
			{
				if (gp.court_trait)
				{
					adults[adult_index].courter_Z[pp - gp.num_env_qtl] = progeny[progeny_index].courter_Z[pp - gp.num_env_qtl];
					adults[adult_index].courter_x[pp - gp.num_env_qtl] = progeny[progeny_index].courter_x[pp - gp.num_env_qtl];
				}
				if (gp.parent_trait)
				{
					adults[adult_index].parent_Z[pp - gp.num_env_qtl] = progeny[progeny_index].parent_Z[pp - gp.num_env_qtl];
					adults[adult_index].parent_x[pp - gp.num_env_qtl] = progeny[progeny_index].parent_x[pp - gp.num_env_qtl];
				}
				if (gp.ind_pref || gp.cor_prefs)
				{
					adults[adult_index].pref_Z[pp - gp.num_env_qtl] = progeny[progeny_index].pref_Z[pp - gp.num_env_qtl];
					adults[adult_index].pref_x[pp - gp.num_env_qtl] = progeny[progeny_index].pref_x[pp - gp.num_env_qtl];
				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
					{
						adults[adult_index].cthresh_Z[pp - gp.num_env_qtl] = progeny[progeny_index].cthresh_Z[pp - gp.num_env_qtl];
						adults[adult_index].cthresh_x[pp - gp.num_env_qtl] = progeny[progeny_index].cthresh_x[pp - gp.num_env_qtl];
					}
					if (gp.parent_conditional || gp.parent_trait)
					{
						adults[adult_index].pthresh_Z[pp - gp.num_env_qtl] = progeny[progeny_index].pthresh_Z[pp - gp.num_env_qtl];
						adults[adult_index].pthresh_x[pp - gp.num_env_qtl] = progeny[progeny_index].pthresh_x[pp - gp.num_env_qtl];
					}
				}
			}
			for (ppp = 0; ppp < gp.num_qtl; ppp++)
			{
				if (gp.court_trait)
				{
					adults[adult_index].courter_Y[pp].per_locus[pp] = progeny[progeny_index].courter_Y[pp].per_locus[ppp];
					adults[adult_index].courter_int[pp].per_locus[ppp] = progeny[progeny_index].courter_int[pp].per_locus[ppp];
				}
				if (gp.parent_trait)
				{
					adults[adult_index].parent_Y[pp].per_locus[pp] = progeny[progeny_index].parent_Y[pp].per_locus[ppp];
					adults[adult_index].parent_int[pp].per_locus[ppp] = progeny[progeny_index].parent_int[pp].per_locus[ppp];
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					adults[adult_index].pref_Y[pp].per_locus[pp] = progeny[progeny_index].pref_Y[pp].per_locus[ppp];
					adults[adult_index].pref_int[pp].per_locus[ppp] = progeny[progeny_index].pref_int[pp].per_locus[ppp];

				}
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
					{
						adults[adult_index].cthresh_Y[pp].per_locus[pp] = progeny[progeny_index].cthresh_Y[pp].per_locus[ppp];
						adults[adult_index].cthresh_int[pp].per_locus[ppp] = progeny[progeny_index].cthresh_int[pp].per_locus[ppp];
					}
					if (gp.parent_conditional || gp.parent_trait)
					{
						adults[adult_index].pthresh_Y[pp].per_locus[pp] = progeny[progeny_index].pthresh_Y[pp].per_locus[ppp];
						adults[adult_index].pthresh_int[pp].per_locus[ppp] = progeny[progeny_index].pthresh_int[pp].per_locus[ppp];
					}
				}
			}
		}
	}
	void pass_on_loci(parameters gp, int adult_index, int progeny_index)
	{
		int pp, ppp;
		for (pp = 0; pp < gp.num_chrom; pp++)
		{
			for (ppp = 0; ppp < gp.num_markers; ppp++)
			{
				adults[adult_index].maternal[pp].loci[ppp] = progeny[progeny_index].maternal[pp].loci[ppp];
				adults[adult_index].paternal[pp].loci[ppp] = progeny[progeny_index].paternal[pp].loci[ppp];
			}
			for (ppp = 0; ppp < gp.qtl_per_chrom[pp]; ppp++)
			{
				if (gp.thresholds_evolve)
				{
					if (gp.courter_conditional || gp.court_trait)
					{
						adults[adult_index].maternal[pp].courter_thresh[ppp] = progeny[progeny_index].maternal[pp].courter_thresh[ppp];
						adults[adult_index].paternal[pp].courter_thresh[ppp] = progeny[progeny_index].paternal[pp].courter_thresh[ppp];
					}
					if (gp.parent_conditional || gp.parent_trait)
					{
						adults[adult_index].maternal[pp].parent_thresh[ppp] = progeny[progeny_index].maternal[pp].parent_thresh[ppp];
						adults[adult_index].paternal[pp].parent_thresh[ppp] = progeny[progeny_index].paternal[pp].parent_thresh[ppp];
					}
				}
				if (gp.court_trait)
				{
					adults[adult_index].maternal[pp].courter_ae[ppp] = progeny[progeny_index].maternal[pp].courter_ae[ppp];
					adults[adult_index].paternal[pp].courter_ae[ppp] = progeny[progeny_index].paternal[pp].courter_ae[ppp];
				}
				if (gp.parent_trait)
				{
					adults[adult_index].maternal[pp].parent_ae[ppp] = progeny[progeny_index].maternal[pp].parent_ae[ppp];
					adults[adult_index].paternal[pp].parent_ae[ppp] = progeny[progeny_index].paternal[pp].parent_ae[ppp];
				}
			}
			if (gp.ind_pref || gp.cor_prefs)
			{
				for (ppp = 0; ppp < progeny[progeny_index].maternal[pp].pref_ae.size(); ppp++)
				{
					adults[adult_index].maternal[pp].pref_ae[ppp] = progeny[progeny_index].maternal[pp].pref_ae[ppp];
					adults[adult_index].paternal[pp].pref_ae[ppp] = progeny[progeny_index].paternal[pp].pref_ae[ppp];
				}
			}
		}
	}
	void adult_from_prog(parameters gp, int adult_index, int progeny_index)
	{
		adults[adult_index].alive = true;
		adults[adult_index].mate_found = 0;
		adults[adult_index].lifetime_rs = 0;
		adults[adult_index].mate_id = int();
		pass_on_loci(gp, adult_index, progeny_index);
		if (progeny[progeny_index].female)
		{
			adults[adult_index].female = true;
			adults[adult_index].pot_rs = gp.max_fecund;
			num_fem++;
		}
		else
		{
			adults[adult_index].female = false;
			num_mal++;
		}
		if (gp.gene_network)
		{
			pass_on_env_qtls(gp, adult_index, progeny_index);
			adults[adult_index].update_traits(gp, courter_thresh, parent_thresh, pref_thresh,
				0, courter_env_qtls, parent_env_qtls, pref_env_qtls, pthresh_env_qtls, cthresh_env_qtls);
		}
		else
			adults[adult_index].update_traits(gp, courter_thresh, parent_thresh, pref_thresh);
	}
	void regulate_popsize(parameters gp)
	{
		int p;
		int num_adults_chosen, prog_alive;
		vector<int> itracker;
		num_mal = num_fem = num_adults_chosen = 0;
		if (num_progeny > gp.carrying_capacity)
		{//only some of the progeny survive
			for (p = 0; p < num_progeny; p++)
			{
				if(progeny[p].alive)
					itracker.push_back(p);				
			}
			random_shuffle(itracker.begin(),itracker.end());
			if (gp.carrying_capacity < itracker.size())
				prog_alive = gp.carrying_capacity;
			else
				prog_alive = itracker.size();
			for (p = 0; p < prog_alive; p++)
			{
				adult_from_prog(gp, p, itracker[p]);
				progeny[itracker[p]].reset(false);
				num_adults_chosen++;
			}
			for (p = prog_alive; p < adults.size(); p++) //designate any remainder as dead.
				adults[p].reset(false);
			for (p = prog_alive; p < progeny.size(); p++)
				progeny[p].reset(false);
		}
		else
		{//all of the progeny survive
			for (p = 0; p < num_progeny; p++)
			{
				if (progeny[p].alive)
				{
					adult_from_prog(gp, num_adults_chosen, p);
					progeny[p].reset(false);
					num_adults_chosen++;
				}
			}
			for (p = num_progeny; p < adults.size(); p++)
				adults[p].reset(false);
			for (p = num_progeny; p < progeny.size(); p++)
				progeny[p].reset(false);
		}
		population_size = num_adults_chosen;
		if (gp.verbose)
			std::cout << ", " << num_adults_chosen << " become adults" << std::flush;
	}
	
	void density_regulation(parameters gp)
	{
		int p, iNumAdultsChosen;
		double CarCapUnfilled, ProgLeft, KeepProb;
		double DRrandnum;
		num_mal = 0;
		num_fem = 0;
		ProgLeft = 0;
		//count the ones that are still alive
		for (p = 0; p < num_progeny; p++)
		{
			if (progeny[p].alive)
				ProgLeft++;
		}
		for (p = 0; p < gp.carrying_capacity; p++)//reset the adults
			adults[p].alive = false;
		CarCapUnfilled = gp.carrying_capacity;
		iNumAdultsChosen = 0;
		for (p = 0; p < num_progeny; p++)
		{
			if (progeny[p].alive)
			{
				if (ProgLeft == 0)
					KeepProb = 0;
				else
					KeepProb = CarCapUnfilled / ProgLeft;
				DRrandnum = genrand();
				if (DRrandnum<KeepProb)
				{//then turn it into an adult
					adults[iNumAdultsChosen].alive = true;
					adults[iNumAdultsChosen].mate_found = 0;
					pass_on_loci(gp, iNumAdultsChosen, p);
					if (gp.gene_network)
					{
						pass_on_env_qtls(gp, iNumAdultsChosen, p);
						adults[iNumAdultsChosen].update_traits(gp, courter_thresh, parent_thresh, pref_thresh,
							0, courter_env_qtls, parent_env_qtls, pref_env_qtls, pthresh_env_qtls, cthresh_env_qtls);
					}
					else
						adults[iNumAdultsChosen].update_traits(gp, courter_thresh, parent_thresh, pref_thresh);
					if (progeny[p].female) 
					{
						adults[iNumAdultsChosen].female = true;
						adults[iNumAdultsChosen].pot_rs = gp.max_fecund;
						num_fem++;
					}
					else 
					{
						adults[iNumAdultsChosen].female = false;
						num_mal++;
					}
					CarCapUnfilled = CarCapUnfilled - 1;
					iNumAdultsChosen++;
				}//end of if KeepProb
				else
					progeny[p].alive = false;
			}//end of if Alive
			ProgLeft = ProgLeft - 1;
		}//end of for p
		population_size = iNumAdultsChosen;
		if (population_size == 0)
			extinct = true;
	}//end Density Regulation

	//calc summary statistics
	void calc_allele_freqs(parameters gp)
	{
		int j, jj, jjj, k;
		for (j = 0; j < gp.num_chrom; j++)
		{
			for (jj = 0; jj < gp.num_markers; jj++)
			{
				double allele_freqs = 0;
				int maj_allele = -1;
				double max_af = 0;
				for (jjj = 0; jjj < gp.num_alleles; jjj++)
				{
					int alive_adults = 0;
					allele_freqs = 0;
					for (k = 0; k < gp.carrying_capacity; k++)
					{
						if (adults[k].alive)
						{
							if (adults[k].maternal[j].loci[jj] == jjj)
								allele_freqs++;
							if (adults[k].paternal[j].loci[jj] == jjj)
								allele_freqs++;
							alive_adults++;
						}
					}
					allele_freqs = allele_freqs / (2 * alive_adults);
					if (allele_freqs > max_af)
					{
						maj_allele = jjj;
						max_af = allele_freqs;
					}
				}
				ref_allele[j].per_locus[jj] = maj_allele;
				maf[j].per_locus[jj] = max_af;
			}
		}
	}
	void eval_rs(parameters gp) // calculate average reproductive success
	{
		int j;
		for (j = 0; j < adults.size(); j++)
		{
			adults[j].lifetime_rs = 0;
		}
		for (j = 0; j < progeny.size(); j++)
		{
			if (progeny[j].alive)
			{
				adults[progeny[j].mom].lifetime_rs++;
				adults[progeny[j].dad].lifetime_rs++;
			}
		}
		if (gp.verbose)
		{
			for (j = 0; j < adults.size(); j++)
			{
				if (adults[j].lifetime_rs > gp.rs_c)
				{
					if (adults[j].alive)
						std::cout << "\n\tLiving ind. " << j;
					else
						std::cout << "\n\tDead ind. " << j;
					if (adults[j].female)
						std::cout << ", female, had lifetime RS " << adults[j].lifetime_rs;
					else
					{
						if (adults[j].courter && adults[j].parent)
							std::cout << ", courter/parent, had lifetime RS " << adults[j].lifetime_rs;
						if (adults[j].courter && !adults[j].parent)
							std::cout << ", courter/non-parent, had lifetime RS " << adults[j].lifetime_rs;
						if (!adults[j].courter && adults[j].parent)
							std::cout << ", non-courter/parent, had lifetime RS " << adults[j].lifetime_rs;
						if (!adults[j].courter && !adults[j].parent)
							std::cout << ", non-courter/non-parent, had lifetime RS " << adults[j].lifetime_rs;
					}
				}
			}
		}
	}
	vector<double> avg_court_rs(parameters gp)
	{
		int k;
		int num_courter, num_parent, num_noncourter, num_nonparent;
		determine_sex_nums(gp); 
		num_courter = calc_freq_courter(gp)*num_mal;
		num_noncourter = num_mal - num_courter;
		num_parent = calc_freq_parent(gp)*num_mal;
		num_nonparent = num_mal - num_parent;
		eval_rs(gp);
		
		vector<double> rs;// courter_rs, noncourter_rs, parent_rs, nonparent_rs
		for (k = 0; k < 4; k++)
			rs.push_back(0);
		for (k = 0; k < adults.size(); k++)
		{
			if (adults[k].alive && !adults[k].female)
			{
				if (gp.courter_conditional || gp.court_trait)
				{
					if (adults[k].courter)
						rs[0] = rs[0] + adults[k].lifetime_rs;
					else
						rs[1] = rs[1] + adults[k].lifetime_rs;
				}
				if (gp.parent_conditional || gp.parent_trait)
				{
					if (adults[k].parent)
						rs[2] = rs[2] + adults[k].lifetime_rs;
					else
						rs[3] = rs[3] + adults[k].lifetime_rs;
				}
			}
		}
		if (gp.courter_conditional || gp.court_trait)
		{
			if (num_courter > 0)
				rs[0] = rs[0] / num_courter;
			else
				rs[0] = 0;
			if(num_noncourter > 0)
				rs[1] = rs[1] / num_noncourter;
			else
				rs[1] = 0;
		}
		if (gp.parent_conditional || gp.parent_trait)
		{
			if (num_parent > 0)
				rs[2] = rs[2] / num_parent;
			else
				rs[2] = 0;
			if (num_nonparent > 0)
				rs[3] = rs[3] / num_nonparent;
			else
				rs[3] = 0;
		}
		return rs;
	}
	
	//output
	void output_qtl_info(parameters gp, ofstream & qtlinfo_output, bool initial)//no environmental qtls
	{
		int j, jj;
		int count1, count2, count3, count4,count5;
		count1 = count2 = count3 =count4= count5 =0;
		for (j = 0; j < gp.num_chrom; j++)
		{
			if (gp.ind_pref || gp.cor_prefs)
			{
				for (jj = 0; jj < pref_qtls[j].per_locus.size(); jj++)
				{
					if (initial)
						qtlinfo_output << "\tPrefQTL" << count1;
					else
						qtlinfo_output << '\t' << j << "." << pref_qtls[j].per_locus[jj];
					count1++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tPrefQTL" << count1;
					else
						qtlinfo_output << "\tNA";
					count1++;
				}
			}
			if (gp.parent_trait)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentQTL" << count2;
					else
						qtlinfo_output << '\t' << j << "." << parent_qtls[j].per_locus[jj];
					count2++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentQTL" << count2;
					else
						qtlinfo_output << "\tNA";
					count2++;
				}
			}
			if (gp.court_trait)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tCourterQTL" << count3;
					else
						qtlinfo_output << '\t' << j << "." << courter_qtls[j].per_locus[jj];
					count3++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tCourterQTL" << count3;
					else
						qtlinfo_output << "\tNA";
					count3++;
				}
			}
			if (courter_thresh_qtls.size() > 0)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tCourterThresholdQTL" << count4;
					else
						qtlinfo_output << '\t' << j << "." << courter_thresh_qtls[j].per_locus[jj];
					count4++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tCourterThresholdQTL" << count4;
					else
						qtlinfo_output << "\tNA";
					count4++;
				}
			}
			if (parent_thresh_qtls.size() > 0)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentThresholdQTL" << count5;
					else
						qtlinfo_output << '\t' << j << "." << parent_thresh_qtls[j].per_locus[jj];
					count5++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentThresholdQTL" << count5;
					else
						qtlinfo_output << "\tNA";
					count5++;
				}
			}
		}
		qtlinfo_output << '\n';
	}
	void output_qtl_info(parameters gp, ofstream & qtlinfo_output, bool initial,
		vector<tracker>& court_env, vector<tracker>& parent_env, vector<tracker>&pref_env,
		vector<tracker>& cthresh_env, vector<tracker>&pthresh_env)
	{
		int j, jj;
		int count1, count2, count3, count4, count5;
		count1 = count2 = count3 = count4 = count5 = 0;
		for (j = 0; j < gp.num_chrom; j++)
		{
			if (gp.ind_pref || gp.cor_prefs)
			{
				for (jj = 0; jj < pref_qtls[j].per_locus.size(); jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tPrefQTL" << count1;
						if(pref_env[j].per_locus[jj]>=0)
							qtlinfo_output << "_env";						
					}
					else
						qtlinfo_output << '\t' << j << "." << pref_qtls[j].per_locus[jj];
					count1++;
				}
			}
			else
			{
				if (gp.ind_pref || gp.cor_prefs)
				{
					for (jj = 0; jj < pref_qtls[j].per_locus.size(); jj++)
					{
						if (initial)
							qtlinfo_output << "\tPrefQTL" << count1;
						else
							qtlinfo_output << "\tNA";
						count1++;
					}
				}
			}
			if (gp.parent_trait)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tParentQTL" << count2;
						if (parent_env[j].per_locus[jj] >= 0)
							qtlinfo_output << "_env";
					}
					else
						qtlinfo_output << '\t' << j << "." << parent_qtls[j].per_locus[jj];
					count2++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentQTL" << count2;
					else
						qtlinfo_output << "\tNA";
					count2++;
				}
			}
			if (gp.court_trait)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tCourterQTL" << count3;
						if (court_env[j].per_locus[jj] >= 0)
							qtlinfo_output << "_env";
					}
					else
						qtlinfo_output << '\t' << j << "." << courter_qtls[j].per_locus[jj];
					count3++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tCourterQTL" << count3;
					else
						qtlinfo_output << "\tNA";
					count3++;
				}
			}
			if (courter_thresh_qtls.size() > 0)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tCourterThresholdQTL" << count4;
						if (cthresh_env[j].per_locus[jj] >= 0)
							qtlinfo_output << "_env";
					}
					else
						qtlinfo_output << '\t' << j << "." << courter_thresh_qtls[j].per_locus[jj];
					count4++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tCourterThresholdQTL" << count4;
					}
					else
						qtlinfo_output << "\tNA";
					count4++;
				}
			}
			if (parent_thresh_qtls.size() > 0)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
					{
						qtlinfo_output << "\tParentThresholdQTL" << count5;
						if (pthresh_env[j].per_locus[jj] >= 0)
							qtlinfo_output << "_env";
					}
					else
						qtlinfo_output << '\t' << j << "." << parent_thresh_qtls[j].per_locus[jj];
					count5++;
				}
			}
			else
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tParentThresholdQTL" << count5;
					else
						qtlinfo_output << "\tNA";
					count5++;
				}
			}
		}
		qtlinfo_output << '\n';
	}
	
	void output_allele_freqs(parameters gp, ofstream & output_file)
	{
		int j, jj;
		calc_allele_freqs(gp);
		for (j = 0; j < gp.num_chrom; j++)
		{
			for (jj = 0; jj < gp.num_markers; jj++)
			{
				output_file << '\t' << maf[j].per_locus[jj];
			}
		}
		output_file << std::flush;
	}
	void output_summary_info(parameters gp, ofstream & summary_output)
	{
		double dtemp;
		vector<double> rs = avg_court_rs(gp);
		if (gp.parent_trait || gp.parent_conditional)
		{
			dtemp = calc_freq_parent(gp);
			summary_output << "\t" << parent_thresh << '\t' << dtemp << '\t' << rs[2] << '\t' << rs[3];
		}
		else
			summary_output << "\tNA\tNA\tNA\tNA";
		if (gp.court_trait || gp.courter_conditional)
		{
			dtemp = calc_freq_courter(gp);
			summary_output << "\t" << courter_thresh << '\t' << dtemp << '\t' << rs[0] << '\t' << rs[1];
		}
		else
			summary_output << "\tNA\tNA\tNA\tNA";
		if ((gp.parent_trait || gp.parent_conditional) && (gp.court_trait || gp.courter_conditional))
		{
			vector<double> morph_freqs = calc_freq_morphs(gp);
			for (int j = 0; j < 4; j++)
				summary_output << '\t' << morph_freqs[j];
		}
		else
			summary_output << "\tNA\tNA\tNA\tNA";
		if (gp.ind_pref || gp.cor_prefs)
		{
			dtemp = calc_freq_pref(gp);
			summary_output << '\t' << pref_thresh << '\t' << dtemp;
		}
		else
			summary_output << "\tNA\tNA";
		output_allele_freqs(gp, summary_output);
		summary_output << std::flush;
	}
	void output_genotypes_vcf(parameters gp, int pop_id)
	{
		int j, jj, jjj;
		if (gp.parent_trait || gp.court_trait || gp.courter_conditional || gp.parent_conditional || gp.ind_pref || gp.cor_prefs)
		{
			string vcf_name = gp.base_name + "_pop_" + to_string(pop_id) + ".vcf";
			ofstream vcf;
			vcf.open(vcf_name);
			calc_allele_freqs(gp);//also determines ref allele
			//output the header
			vcf << "##fileformat=VCFv4.0";
			string date = determine_date();
			vcf << "\n##fileDate=" << date;
			vcf << "\n##source=ARTsSimulation";
			vcf << "\n##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">";
			vcf << "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
			vcf << "\n##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
			for (j = 0; j < population_size; j++)
			{
				if (adults[j].alive)
				{
					if (adults[j].female)
						vcf << "\tFEM" << j << "_pref" << adults[j].female_pref;
					else
					{
						if (adults[j].courter && adults[j].parent)
							vcf << "\tMAL" << j << "_CP";
						if (adults[j].courter && !adults[j].parent)
							vcf << "\tMAL" << j << "_C";
						if (!adults[j].courter && adults[j].parent)
							vcf << "\tMAL" << j << "_P";
						if (!adults[j].courter && !adults[j].parent)
							vcf << "\tMAL" << j << "_NON";
					}
				}
			}
			for (j = 0; j < gp.num_chrom; j++)
			{
				for (jj = 0; jj < gp.num_markers; jj++)
				{
					vcf << '\n' << j << '\t' << jj << '\t' << j << "." << jj << '\t' << ref_allele[j].per_locus[jj];
					int cnt = 0;//output the reference alleles
					for (jjj = 0; jjj < gp.num_alleles; jjj++)
					{
						if (ref_allele[j].per_locus[jj] != jjj)
						{
							if (cnt == 0)
								vcf << "\t" << jjj;
							else
								vcf << "," << jjj;
							cnt++;
						}
					}
					vcf << "\t100\tPASS\tAF=" << maf[j].per_locus[jj] << "\tGT";
					for (jjj = 0; jjj < population_size; jjj++)
					{
						if (adults[jjj].alive)
						{
							vcf << '\t' << adults[jjj].maternal[j].loci[jj] << "/" << adults[jjj].paternal[j].loci[jj];
						}
					}
				}
			}
			vcf.close();
		}
	}
	void output_trait_info(parameters gp, int pop_id, ofstream & output)
	{
		int j;
		eval_rs(gp);
		for (j = 0; j < adults.size(); j++)
		{
			output << '\n' << pop_id << '\t' << j;// "Pop\tIndividual\tSex\tCourter\tCourtTrait\tParent\tParentTrait\tPreference\tPrefTrait\tMateFound\tPotRS\tLifetimeRS\tAlive";
			if (adults[j].female)
				output << "\tFEMALE";
			else
				output << "\tMALE";
			output << '\t' << adults[j].courter << '\t' << adults[j].courter_trait << '\t' << adults[j].parent << '\t' << adults[j].parent_trait 
				<< '\t' << adults[j].female_pref << '\t' << adults[j].pref_trait << '\t' << adults[j].mate_found << '\t' << adults[j].pot_rs 
				<< '\t' << adults[j].lifetime_rs <<'\t' <<adults[j].alive;
		}
	}
};
