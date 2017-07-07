#pragma once

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
	int population_size, num_mal, num_fem, num_progeny, sex_trait, max_num_migrants, migrant_index;
	double sex_theta, sex_omega, courter_thresh, parent_thresh;
	bool extinct, random_mating;
	vector<double> theta, mean_fem_traits, mean_mal_traits, d_parentfreq, d_courterfreq;
	vector<individual> adults;
	vector<individual> progeny;
	vector<tracker> courter_qtls,parent_qtls,pref_qtls,courter_env_qtls,parent_env_qtls,pref_env_qtls, maf, hs, maj_allele, ref_allele;


	population()
	{
		population_size = num_mal = num_fem = num_progeny = sex_trait = max_num_migrants = migrant_index = int();
		theta = mean_fem_traits = mean_mal_traits = vector<double>();
		sex_theta = sex_omega = courter_thresh = parent_thresh = double();
		adults = progeny = vector<individual>();
		courter_env_qtls = courter_qtls = parent_env_qtls = parent_qtls = pref_env_qtls = pref_qtls = maf = hs = vector<tracker>();
		random_mating = false;
		extinct = false;
		sgenrand(time(0));
	}

	void initialize_supergene(parameters gp)
	{
		int m, mm, mmm;
		double supergene_nmarkers = gp.supergene_prop*gp.num_markers;
		//chose a chromosome for the supergene
		int sg_chrom = randnum(gp.num_chrom);
		for (m = 0; m < gp.num_chrom; m++)
		{
			if (sg_chrom == m)
				gp.qtl_per_chrom[m] = gp.num_qtl;
			else
				gp.qtl_per_chrom[m] = 0;
		}
		//choose marker start
		int sg_start = randnum(gp.num_markers - supergene_nmarkers);
		if (gp.court_trait && gp.parent_trait)
		{
			supergene_nmarkers = 2 * supergene_nmarkers;
			if (supergene_nmarkers <= 2 * gp.num_qtl)//make sure there's enough space in the supergene region.
				supergene_nmarkers = 2 * gp.num_qtl;
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
		else
		{
			if (supergene_nmarkers <= gp.num_qtl)//make sure there's enough space in the supergene region.
				supergene_nmarkers = gp.num_qtl;
			for (m = 0; m < gp.num_chrom; m++)
			{
				if (gp.court_trait)
					courter_qtls.push_back(tracker());
				if (gp.parent_trait)
					parent_qtls.push_back(tracker());
				for (mm = 0; mm < gp.qtl_per_chrom[m]; mm++)
				{
					if (gp.court_trait)
						courter_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
					if (gp.parent_trait)
						parent_qtls[m].per_locus.push_back((sg_start + randnum(supergene_nmarkers)));
				}
			}
		}
	}
	void initialize_adults(parameters gp)
	{
		int j, jj, jjj;
		for (j = 0; j < (gp.carrying_capacity + max_num_migrants); j++)
		{
			adults.push_back(individual());
			if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
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
				if (gp.cor_prefs || gp.ind_pref || !random_mating)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].pref_ae.push_back(double());
						adults[j].paternal[jj].pref_ae.push_back(double());
					}
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
		for (j = 0; j < (max_potential_fecund*(gp.carrying_capacity + max_num_migrants)); j++)
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
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						progeny[j].maternal[jj].parent_ae.push_back(double());
						progeny[j].paternal[jj].parent_ae.push_back(double());
					}
				}
				if (gp.ind_pref || gp.cor_prefs || !random_mating)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						progeny[j].maternal[jj].pref_ae.push_back(double());
						progeny[j].paternal[jj].pref_ae.push_back(double());
					}
				}
			}
		}
	}
	void assign_env_qtls(vector<tracker>&env_qtls, parameters gp)
	{
		int j, jj, jjj;
		//assign qtls to be environmental.
		for (j = 0; j < gp.num_chrom; j++)
		{
			env_qtls.push_back(tracker());
			for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
			{
				env_qtls[j].per_locus.push_back(-1);
			}
		}
		int index = 0;
		int this_qtl = randnum(gp.num_qtl);
		for (jj = 0; jj < gp.num_chrom; jj++)
		{
			for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
			{
				index++;
				if (index == this_qtl)
					env_qtls[jj].per_locus[jjj] = this_qtl;
			}
		}
	}
	void initialize(parameters gp)
	{
		cout << "Initializing a population.\n";
		population_size = gp.carrying_capacity;
		if (!gp.cor_prefs && !gp.ind_pref)
			random_mating = true;
		int j, jj, jjj, k;
		num_fem = num_mal = 0;
		migrant_index = gp.carrying_capacity;
		sex_theta = -100;
		if (gp.supergene)
		{
			initialize_supergene(gp);
		}
		else
		{
			//set up the qtls for the traits
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
		}
		//set up females
		if (gp.cor_prefs == true)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				pref_qtls.push_back(tracker());
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
					pref_qtls[j].per_locus.push_back(courter_qtls[j].per_locus[jj]);
			}
		}
		if (gp.ind_pref == true)//female preference is not a part of a supergene if it's independent
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				pref_qtls.push_back(tracker());
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
					pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
			}
		}
		if (pref_qtls.size() == 0 && !random_mating)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				pref_qtls.push_back(tracker());
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
					pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
			}
			gp.ind_pref = true;
		}
		//if environmental effects
		if (gp.env_effects)
		{
			//randomly choose gp.num_env_qtl from the set
			for (j = 0; j < gp.num_env_qtl; j++)
			{
				if (gp.court_trait)
					assign_env_qtls(courter_env_qtls, gp);
				if (gp.parent_trait)
					assign_env_qtls(parent_env_qtls, gp);
				if (!random_mating)
					assign_env_qtls(pref_env_qtls, gp);
			}
		}

		//set up trackers
		for (j = 0; j < gp.num_chrom; j++)
		{
			maf.push_back(tracker());
			hs.push_back(tracker());
			for (jj = 0; jj < gp.num_markers; jj++)
			{
				maf[j].per_locus.push_back(double());
				hs[j].per_locus.push_back(double());
			}
		}
		//set up adults
		initialize_adults(gp);
		//set up progeny
		initialize_progeny(gp);
		//assign loci and alleles
		vector<double> tempallele1, tempallele2, tempallele3;
		for (j = 0; j < gp.num_alleles; j++) //set up specific alleles.
		{
			tempallele1.push_back(randnorm(0, gp.allelic_std_dev));
			tempallele2.push_back(randnorm(0, gp.allelic_std_dev));
			tempallele3.push_back(randnorm(0, gp.allelic_std_dev));
		}
		for (j = 0; j < (gp.carrying_capacity + max_num_migrants); j++)
		{
			if (j < gp.carrying_capacity)
				adults[j].alive = true;
			else
				adults[j].alive = false;//it's an empty slot for migrants
			if (gp.env_effects)//set up interactions
			{
				if (gp.court_trait == true)
				{
					for (jj = 0; jj < gp.num_qtl; jj++)
					{//the first ones are the environmentally regulated ones.
						adults[j].courter_int.push_back(tracker());//SxS matrix of gene interactions
						adults[j].courter_Y.push_back(tracker());	//SxS matrix of whether genes regulate each other					
						if (jj > gp.num_env_qtl)
						{
							adults[j].courter_Z.push_back(1);		//(S-E) vector of whether gene j contributes to trait i. 
							adults[j].courter_x.push_back(double());//(S-E) vector of weigths and contributions to trait i.
						}
						for (jjj = 0; jjj < gp.num_qtl; jjj++)
						{
							adults[j].courter_Y[jj].per_locus.push_back(1);
							adults[j].courter_int[jj].per_locus.push_back(randnorm(0, 0.5));
						}
					}
				}
				if (gp.parent_trait == true)
				{
					for (jj = 0; jj < gp.num_qtl; jj++)
					{
						adults[j].parent_int.push_back(tracker());
						adults[j].parent_Y.push_back(tracker());
						if (jj > gp.num_env_qtl)
						{
							adults[j].parent_Z.push_back(1);
							adults[j].parent_x.push_back(double());
						}
						for (jjj = 0; jjj < gp.num_qtl; jjj++)
						{
							adults[j].parent_Y[jj].per_locus.push_back(1);//they all start out with an interaction
							adults[j].parent_int[jj].per_locus.push_back(randnorm(0, 0.5));
						}
					}
				}
				if (!random_mating)
				{
					for (jj = 0; jj < gp.num_qtl; jj++)
					{
						adults[j].pref_int.push_back(tracker());
						adults[j].pref_Y.push_back(tracker());						
						if (jj > gp.num_env_qtl)
						{
							adults[j].pref_Z.push_back(1);
							adults[j].pref_x.push_back(double());
						}
						for (jjj = 0; jjj < gp.num_qtl; jjj++)
						{
							adults[j].pref_Y[jj].per_locus.push_back(1);
							adults[j].pref_int[jj].per_locus.push_back(randnorm(0, 0.5));
						}
					}
				}
			}
			
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				if (gp.court_trait || !random_mating)//because if preferences exist then there needs to be a male trait.
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
				if (!random_mating)
				{
					for (jjj = 0; jjj < gp.qtl_per_chrom[jj]; jjj++)
					{
						adults[j].maternal[jj].pref_ae[jjj] = tempallele3[j%gp.num_alleles];
						adults[j].paternal[jj].pref_ae[jjj] = tempallele3[j%gp.num_alleles];
					}
				}
				int index = 0;
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci[jjj] = j%gp.num_alleles;
					adults[j].paternal[jj].loci[jjj] = j%gp.num_alleles;
				}
			}
			if (gp.court_trait || !random_mating) //because if preferences exist then there needs to be a male trait.
			{
				adults[j].calc_courter_trait(gp); //to initialize env cue = 0 (the default)
			}
			if (gp.parent_trait)
			{
				adults[j].calc_parent_trait(gp);
			}
			if (genrand() < 0.5)
			{
				adults[j].female = true;
				num_fem++;
			}
			else
			{
				adults[j].female = false;
				num_mal++;
			}
		}//carrying_capacity
		//set up morphs
		if (gp.court_trait)
		{
			courter_thresh = calc_mean_courter_ae(gp);
			//assign morph to each ind
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (!adults[j].female)
				{
					if (adults[j].courter_trait < courter_thresh)
						adults[j].courter = false;
					else
						adults[j].courter = true;
				}
			}
		}
		else
			courter_thresh = 0;
		if (gp.parent_trait)
		{
			parent_thresh = calc_mean_parent_ae(gp);
			//assign morph to each ind
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (!adults[j].female)
				{
					if (adults[j].courter_trait < parent_thresh)
						adults[j].parent = false;
					else
						adults[j].parent = true;
				}
			}
		}
		else
			parent_thresh = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			adults[j].update_traits(gp, courter_thresh, parent_thresh);
		}
		population_size = gp.carrying_capacity;
		
	}

	void determine_pop_size(parameters gp)
	{
		int j;
		population_size = 0;
		for (j = 0; j < (gp.carrying_capacity + max_num_migrants); j++)
		{
			if (adults[j].alive)
				population_size++;
		}
		if (population_size == 0)
			cout << "\nPopulation has crashed!";
	}

	//set up male traits
	double calc_mean_courter_ae(parameters gp)
	{
		int j, jj, jjj, count;
		double mean;
		mean = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				adults[j].calc_courter_trait(gp);
				mean = mean + adults[j].courter_trait;
				count++;
			}
		}
		mean = mean / count;		
		return(mean);
	}
	double calc_mean_parent_ae(parameters gp)
	{
		int j, jj, jjj, count;
		double mean;
		mean = count = 0;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && !adults[j].female)
			{
				adults[j].calc_parent_trait(gp);
				mean = mean + adults[j].parent_trait;
				count++;
			}
		}
		mean = mean / count;
		return(mean);
	}
	double calc_freq_courter(parameters gp)
	{
		int j, jj;
		double mean, freqT, count;
		mean = freqT = count = 0;
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
		int j, jj;
		double mean, freqT, count;
		mean = freqT = count = 0;
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
	void parent_fd_rs(parameters gp)
	{
		int j, jj;
		double mean, freqT, count;
		freqT = calc_freq_parent(gp);
		mean = 0;
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
	void pref_femtrait(parameters gp)//if it's correlated, then preference should be based on female's own trait
	{
		for (int j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				if(gp.court_trait)
					adults[j].calc_preference_trait(gp, courter_thresh);
				else
					adults[j].calc_preference_trait(gp, parent_thresh);
			}
		}
	}
	void pref_fd_courter(parameters gp)
	{
		int j, jj;
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
		int j, jj;
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
	void pref_cd_courter_noenv(parameters gp, double meanT, double meanF, double xswitch)
	{
		int j, jj, count;
		double mean;
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				adults[j].calc_preference_trait(gp, courter_thresh);//consider environmental cue...
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
		int j, jj, count;
		double mean;
		if (!random_mating) //only waste the time if you're gonna use a preference
		{
			if (!gp.FD_pref && !gp.CD_pref) //then it's gaussian
			{
				if (gp.ind_pref)//if it's independent, not freq dependent, and not condition dependent
								//then it's based on mean male trait
					pref_better_morph(gp);
				else //if it's correlated, then preference should be based on female's own trait
					pref_femtrait(gp);
			}
			if (gp.FD_pref)//then need the frequency of morphs and female prefers less frequent one.
			{
				double freqT = 0;
				count = 0;
				if (gp.court_trait)
					pref_fd_courter(gp);
				else //it's based on the parent trait
					pref_fd_parent(gp);
				
			}
			if (gp.CD_pref)//then need to calculate female fecundity and compare it to S(rs_p - rs_np).
			{
				int j, jj, count;
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
				pref_cd_courter_noenv(gp, meanT, meanF, xswitch);
			}
		}
	}

	//functions for mating
	bool recombine_chromosome(chromosome& chrom, individual &parent, double expected_recomb_events, int which_chrom, parameters gp)//tracker& court_qtls, tracker& parent_qtls, tracker& pref_qtls, 
	{
		int RCi, RCj, RCk;
		int NumberRecombEvents = 0;
		int segment_start[22], segment_end[20];
		int break_point[20];
		bool recombined = false;

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
			recombined = true;
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
				if (segment_maternal[RCi])
				{
					for (RCj = segment_start[RCi]; RCj < segment_end[RCi]; RCj++)
					{
						int new_loc = parent.maternal[which_chrom].loci[RCj];
						chrom.loci[RCj] = parent.maternal[which_chrom].loci[RCj];
					}
				}
				else
				{
					for (RCj = segment_start[RCi]; RCj < segment_end[RCi]; RCj++)
					{
						int new_loc = parent.paternal[which_chrom].loci[RCj];
						chrom.loci[RCj] = parent.paternal[which_chrom].loci[RCj];
					}
				}
			}
			//now the QTLs
			if (gp.court_trait)
			{
				for (RCi = 0; RCi < gp.qtl_per_chrom[which_chrom]; RCi++)
				{
					for (RCj = 0; RCj < num_segments; RCj++)
					{
						if (courter_qtls[which_chrom].per_locus[RCi] >= segment_start[RCj] && courter_qtls[which_chrom].per_locus[RCi] >= segment_end[RCj])
						{
							if (segment_maternal[RCj])
								chrom.courter_ae[RCi] = parent.maternal[which_chrom].courter_ae[RCi];
							else
								chrom.courter_ae[RCi] = parent.paternal[which_chrom].courter_ae[RCi];
						}
					}
				}
			}
			if (gp.parent_trait)
			{
				for (RCi = 0; RCi < gp.qtl_per_chrom[which_chrom]; RCi++)
				{
					for (RCj = 0; RCj < num_segments; RCj++)
					{
						if (parent_qtls[which_chrom].per_locus[RCi] >= segment_start[RCj] && parent_qtls[which_chrom].per_locus[RCi] >= segment_end[RCj])
						{
							if (segment_maternal[RCj])
								chrom.parent_ae[RCi] = parent.maternal[which_chrom].parent_ae[RCi];
							else
								chrom.parent_ae[RCi] = parent.paternal[which_chrom].parent_ae[RCi];
						}
					}
				}
			}
			if (!random_mating)
			{
				for (RCi = 0; RCi < gp.qtl_per_chrom[which_chrom]; RCi++)
				{
					for (RCj = 0; RCj < num_segments; RCj++)
					{
						if (pref_qtls[which_chrom].per_locus[RCi] >= segment_start[RCj] && pref_qtls[which_chrom].per_locus[RCi] >= segment_end[RCj])
						{
							if (segment_maternal[RCj])
								chrom.pref_ae[RCi] = parent.maternal[which_chrom].pref_ae[RCi];
							else
								chrom.pref_ae[RCi] = parent.paternal[which_chrom].pref_ae[RCi];
						}
					}
				}
			}
		}//numberrecombevents > 0
		else
		{
			recombined = false;
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
					if (!random_mating)
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
					if (!random_mating)
					{
						chrom.pref_ae[RCi] = parent.paternal[which_chrom].pref_ae[RCi];
						chrom.pref_ae[RCi] = parent.paternal[which_chrom].pref_ae[RCi];
					}
				}

			}
		}
		return recombined;
	}
	double making_babies(parameters gp, int fecundity, int& num_progeny, int mom_index, int dad_index)
	{
		int j, jj, jjj, rec_val;
		double rel_rs, diff;
		bool recomb1, recomb2;
		rel_rs = 0;
		if (num_progeny >= (gp.max_fecund*gp.carrying_capacity))
			num_progeny = (gp.max_fecund*gp.carrying_capacity) - 1;
		for (jj = 0; jj < fecundity; jj++)
		{

			for (jjj = 0; jjj < gp.num_chrom; jjj++)
			{
				recomb1 = recombine_chromosome(progeny[num_progeny].maternal[jjj], adults[mom_index], gp.recombination_rate, jjj, gp);
				recomb2 = recombine_chromosome(progeny[num_progeny].paternal[jjj], adults[dad_index], gp.recombination_rate, jjj, gp);
			}
			//also need to do mutation
			if (gp.env_effects)//need to evaluate if it's compatible
			{
				progeny[num_progeny].alive = true;
				progeny[num_progeny].mutation_env(gp);
				if(gp.court_trait)
					progeny[num_progeny].calc_courter_trait(gp, 0);
				if (gp.parent_trait)
					progeny[num_progeny].calc_parent_trait(gp, 0);
				if (!random_mating)
					progeny[num_progeny].calc_preference_trait(gp, 0);
			}
			else
			{
				progeny[num_progeny].mutation(gp, courter_qtls, parent_qtls, pref_qtls);
				progeny[num_progeny].update_traits(gp, courter_thresh, parent_thresh);
				progeny[num_progeny].alive = true;
			}
			if (progeny[num_progeny].alive)
			{
				if (genrand() < 0.5)
					progeny[num_progeny].female = true;
				else
					progeny[num_progeny].female = false;
				rel_rs++;
				num_progeny++;
			}
		}
		rel_rs = rel_rs / fecundity;
		return rel_rs;
	}
	void gaussian_mating(bool output, string out_name, parameters gp)
	{
		//monogamous mating without choice
		int j, jj, jjj, male_id, encounters;
		double mate_prob, rel_rs;
		num_progeny = num_fem = num_mal = 0;
		bool mate_found;
		vector<int> male_index;
		vector<bool> mated;
		ofstream out_file;
		if (output)
		{
			cout << "\nWriting to file " << out_name;
			out_file.open(out_name);
		}
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
					mated.push_back(false);
				}
			}
			adults[j].mate_found = 0;
		}
		

		for (j = 0; j < adults.size(); j++)
		{
			mate_found = false;
			if (adults[j].female && adults[j].alive)
			{
				int count = 0;
				encounters = 0;
				if (random_mating)
				{
					while (!mate_found && encounters <= gp.max_encounters)//female mates once
					{
						int irndnum = randnum(num_mal);
						male_id = male_index[irndnum];
						if (adults[male_id].alive)
						{
							if (gp.polygyny || !mated[irndnum])
							{
								mate_found = true;
								mated[irndnum] = true;
								adults[j].mate_found++;
								adults[male_id].mate_found++;
							}
							encounters++;
						}
					}
				}//random mating
				else
				{ //gaussian: already calculated male traits and female prefs 
					while (!mate_found && encounters <= gp.max_encounters)
					{
						int irndnum = randnum(num_mal);
						male_id = male_index[irndnum];
						if (adults[male_id].alive)
						{
							if (gp.polygyny || !mated[irndnum])
							{
								if (gp.court_trait)
									mate_prob = exp(-0.5 * (adults[male_id].courter_trait - gp.gaussian_pref_mean)*
									(adults[male_id].courter_trait - gp.gaussian_pref_mean) / adults[j].female_pref);
								else
									mate_prob = exp(-0.5 * (adults[male_id].parent_trait - gp.gaussian_pref_mean)*
									(adults[male_id].parent_trait - gp.gaussian_pref_mean) / adults[j].female_pref);
								if (genrand() < mate_prob)
								{
									mate_found = true;
									mated[irndnum] = true;
									adults[j].mate_found++;
									adults[male_id].mate_found++;
									if (gp.CD_court)
									{
										if (adults[male_id].courter)//it decreases
											adults[male_id].courter_trait = adults[male_id].courter_trait - gp.cond_adj;
										else//it increases
											adults[male_id].courter_trait = adults[male_id].courter_trait + gp.cond_adj;
										adults[male_id].assign_court_morph(gp, courter_thresh);
									}
									if (gp.CD_parent)
									{
										if (adults[male_id].parent)//it decreases
											adults[male_id].parent_trait = adults[male_id].parent_trait - gp.cond_adj;
										else//it increases
											adults[male_id].parent_trait = adults[male_id].parent_trait + gp.cond_adj;
										adults[male_id].assign_parent_morph(gp, parent_thresh);
									}
								}
							}
							encounters++;
						}
					}//while
				}				
				if (mate_found)
				{
					int fecundity = min(adults[j].pot_rs, adults[male_id].pot_rs);
					rel_rs = making_babies(gp, fecundity, num_progeny, j, male_id);
					if (output)
					{
						out_file << '\n' << rel_rs << "\tMother" << '\t' << j;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++) {
								out_file << '\t' << adults[j].maternal[jjj].loci[x] << "/" << adults[j].paternal[jjj].loci[x];

							}
						}
						out_file << "\nFather" << '\t' << male_id;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++)
								out_file << '\t' << adults[male_id].maternal[jjj].loci[x] << "/" << adults[male_id].paternal[jjj].loci[x];
						}
						out_file << "\nOffspring" << '\t' << num_progeny;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++) {
								out_file << '\t' << progeny[num_progeny].maternal[jjj].loci[x] << "/" << progeny[num_progeny].paternal[jjj].loci[x];

							}
						}
					}
				}
			}
		}
		if (output)
			out_file.close();
	}//mating
	void bestofN_mating(bool output, string out_name, parameters gp)
	{
		//monogamous mating without choice
		int j, jj, jjj, male_id, encounters, irndnum;
		double mate_prob, rel_rs;
		num_progeny = num_fem = num_mal = 0;
		bool mate_found;
		vector<int> male_index;
		vector<bool> mated;
		ofstream out_file;
		if (output)
		{
			cout << "\nWriting to file " << out_name;
			out_file.open(out_name);
		}
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
					mated.push_back(false);
				}
			}
			adults[j].mate_found = 0;
		}


		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].female && adults[j].alive)
			{
				int count = 0;
				encounters = 0;
				mate_found = false;
				if (random_mating)
				{
					while (!mate_found && encounters < gp.max_encounters)//female mates once
					{
						irndnum = randnum(num_mal);
						male_id = male_index[irndnum];
						if (adults[male_id].alive)
						{
							if (gp.polygyny || !mated[irndnum])//either polygyny is ok or if monogamy the male hasn't mated yet
							{
								mate_found = true;
								mated[irndnum] = true;
								adults[j].mate_found++;
								adults[male_id].mate_found++;
							}
						}
						encounters++;
					}
				}//random mating
				else
				{ //preference for one morph or the other 
					vector <int> acceptable_males;
					while (encounters < gp.max_encounters)
					{
						irndnum = randnum(num_mal);
						male_id = male_index[irndnum];
						if (adults[male_id].alive)
						{
							if (gp.polygyny || !mated[irndnum])
							{
								if (gp.court_trait)
								{
									if (adults[j].female_pref == adults[male_id].courter);
										acceptable_males.push_back(male_id);
								}
								else//parent_trait
								{
									if (adults[j].female_pref == adults[male_id].parent);
										acceptable_males.push_back(male_id);
								}
							}
						}
						encounters++;
					}
					if (acceptable_males.size() >= 1)//if there are multiple males, she chooses one randomly
					{
						int randbest = randnum(acceptable_males.size());
						male_id = acceptable_males[randbest];
						mate_found = true;
						mated[irndnum] = true;
						adults[j].mate_found++;
						adults[male_id].mate_found++;
						if (gp.CD_court)
						{
							if (adults[male_id].courter)//it decreases
								adults[male_id].courter_trait = adults[male_id].courter_trait - gp.cond_adj;
							else//it increases
								adults[male_id].courter_trait = adults[male_id].courter_trait + gp.cond_adj;
							adults[male_id].assign_court_morph(gp, courter_thresh);
						}
						if (gp.CD_parent)
						{
							if (adults[male_id].parent)//it decreases
								adults[male_id].parent_trait = adults[male_id].parent_trait - gp.cond_adj;
							else//it increases
								adults[male_id].parent_trait = adults[male_id].parent_trait + gp.cond_adj;
							adults[male_id].assign_parent_morph(gp, parent_thresh);
						}
					}
				}//end of finding the mates
				if (mate_found)
				{
					int fecundity = min(adults[j].pot_rs, adults[male_id].pot_rs);
					rel_rs = making_babies(gp, fecundity, num_progeny, j, male_id);
					if (output)
					{
						out_file << '\n' << rel_rs << "\tMother" << '\t' << j;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++) {
								out_file << '\t' << adults[j].maternal[jjj].loci[x] << "/" << adults[j].paternal[jjj].loci[x];

							}
						}
						out_file << "\nFather" << '\t' << male_id;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++)
								out_file << '\t' << adults[male_id].maternal[jjj].loci[x] << "/" << adults[male_id].paternal[jjj].loci[x];
						}
						out_file << "\nOffspring" << '\t' << num_progeny;
						for (jjj = 0; jjj < gp.num_chrom; jjj++)
						{
							for (int x = 0; x < gp.num_markers; x++) {
								out_file << '\t' << progeny[num_progeny].maternal[jjj].loci[x] << "/" << progeny[num_progeny].paternal[jjj].loci[x];

							}
						}
					}
				}
			}
		}
		if (output)
			out_file.close();
	}
	void nest_and_fertilize(parameters gp, bool output, string out_name)
	{
		int j, jj, jjj;
		num_progeny = num_fem = num_mal = 0;
		bool nest_found;
		vector<int> male_index;
		ofstream out_file;

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
				}
			}
			adults[j].mate_found = 0;
		}
		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].female && adults[j].alive) //loop through the females.
			{
				nest_found = choose_nest(j, male_index, gp);
				if (nest_found)
				{
					random_fertilization(j, male_index, gp);
				}
			}
		}
	}
	bool choose_nest(int fem_index, vector<int> & male_index, parameters gp)
	{//females choose one male to give her eggs to.
		int male_id, encounters, irndnum;
		bool mate_found;
		
		int count = 0;
		encounters = 0;
		mate_found = false;
		if (random_mating)
		{
			while (!mate_found && encounters < gp.max_encounters)//female mates once
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
		}//random mating
		else
		{ //preference for one morph or the other 
			vector <int> acceptable_males;
			while (encounters < gp.max_encounters)
			{
				irndnum = randnum(num_mal);
				male_id = male_index[irndnum];
				if (adults[male_id].alive)
				{
					if (gp.polygyny || adults[male_id].mate_found == 0)
					{
						if (gp.court_trait)
						{
							if (adults[fem_index].female_pref == adults[male_id].courter);
							acceptable_males.push_back(male_id);
						}
						else//parent_trait
						{
							if (adults[fem_index].female_pref == adults[male_id].parent);
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
					adults[male_id].assign_court_morph(gp, courter_thresh);
				}
				if (gp.CD_parent)
				{
					if (adults[male_id].parent)//it decreases
						adults[male_id].parent_trait = adults[male_id].parent_trait - gp.cond_adj;
					else//it increases
						adults[male_id].parent_trait = adults[male_id].parent_trait + gp.cond_adj;
					adults[male_id].assign_parent_morph(gp, parent_thresh);
				}
			}
		}//end of finding the mates
		return mate_found;
	}
	void random_fertilization(int fem_id, vector<int> & male_index, parameters gp)
	{
		int k,fecundity, randmale, sneaker_id;
		double rel_rs;
		int male_id = adults[fem_id].mate_id;
		int num_eggs_left = adults[fem_id].pot_rs;
		
		//nesting male gets the first shot
		fecundity = randnum(adults[male_id].pot_rs);
		if (fecundity > num_eggs_left)
			fecundity = num_eggs_left;
		rel_rs = making_babies(gp, fecundity, num_progeny, fem_id, male_id);
		num_eggs_left = num_eggs_left - fecundity;
		//then if there are any eggs left, sneakers can fertilize
		while (num_eggs_left > 0)
		{
			randmale = randnum(male_index.size());
			if (!adults[male_index[randmale]].parent)//if it's a sneaker
			{
				sneaker_id = male_index[randmale];
				fecundity = randnum(adults[sneaker_id].pot_rs);
				if(fecundity > num_eggs_left)
					fecundity = num_eggs_left;
				rel_rs = making_babies(gp, fecundity, num_progeny, fem_id, sneaker_id);
				num_eggs_left = num_eggs_left - fecundity;
				adults[sneaker_id].mate_found++;
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
		double dOptimum;
		double phenSD = 0;
		double phenMean = 0;
		double num;
		int malecount = 0;
		ProgAlive = 0;

		for (j = 0; j < num_progeny; j++) 
		{
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
					if (gp.court_trait && !gp.FD_court)
						dSurvProb = dSurvProb*via_against_courter(gp.via_sel_strength, j);
					if(gp.FD_court)
						dSurvProb = dSurvProb*via_fd_courter(gp, j);
					if (gp.parent_trait && !gp.FD_court)
						dSurvProb = dSurvProb*via_against_parent(gp.via_sel_strength, j);
					if (gp.FD_parent)
						dSurvProb = dSurvProb*via_fd_parent(gp, j);					
				}
				else
					dSurvProb = 1;
				//cout<<dSurvProb<<'\n';
				drnum1 = genrand();
				if (drnum1 < dSurvProb)
				{
					progeny[j].alive = true;
					ProgAlive++;
				}
				else
					progeny[j].alive = false;
			}
		} // end of j
	}

	//stochastic survival
	void density_regulation(parameters gp)
	{
		int p, pp, ppp;
		int iNumAdultsChosen;
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
					for (pp = 0; pp < gp.num_chrom; pp++)
					{
						for (ppp = 0; ppp < gp.num_markers; ppp++)
						{
							adults[iNumAdultsChosen].maternal[pp].loci[ppp] = progeny[p].maternal[pp].loci[ppp];
							adults[iNumAdultsChosen].paternal[pp].loci[ppp] = progeny[p].paternal[pp].loci[ppp];
						}
						for (ppp = 0; ppp < gp.num_qtl; ppp++)
						{
							if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
							{
								adults[iNumAdultsChosen].maternal[pp].courter_ae[ppp] = progeny[p].maternal[pp].courter_ae[ppp];
								adults[iNumAdultsChosen].paternal[pp].courter_ae[ppp] = progeny[p].paternal[pp].courter_ae[ppp];
							}
							if (gp.parent_trait)
							{
								adults[iNumAdultsChosen].maternal[pp].parent_ae[ppp] = progeny[p].maternal[pp].parent_ae[ppp];
								adults[iNumAdultsChosen].paternal[pp].parent_ae[ppp] = progeny[p].paternal[pp].parent_ae[ppp];
							}
							if (!random_mating)
							{
								adults[iNumAdultsChosen].maternal[pp].pref_ae[ppp] = progeny[p].maternal[pp].pref_ae[ppp];
								adults[iNumAdultsChosen].paternal[pp].pref_ae[ppp] = progeny[p].paternal[pp].pref_ae[ppp];
							}
						}
					}
					adults[iNumAdultsChosen].update_traits(gp, courter_thresh, parent_thresh);
					if (progeny[p].female) {
						adults[iNumAdultsChosen].female = true;
						num_fem++;
					}
					else {
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

	//output
	void output_qtl_info(parameters gp, ofstream & qtlinfo_output, bool initial)
	{
		int j, jj;
		int count1, count2, count3;
		count1 = count2 = count3 = 0;
		for (j = 0; j < gp.num_chrom; j++)
		{
			if (!random_mating)
			{
				for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
				{
					if (initial)
						qtlinfo_output << "\tPrefQTL" << count1;
					else
						qtlinfo_output << '\t' << j << "." << pref_qtls[j].per_locus[jj];
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
		}
		qtlinfo_output << '\n';
	}
	void output_allele_freqs(parameters gp, ofstream & output_file)
	{
		int j, jj, jjj;
		double allele_freq;
		for (j = 0; j < gp.num_chrom; j++)
		{
			for (jj = 0; jj < gp.num_markers; jj++)
			{
				allele_freq = 0;
				for (jjj = 0; jjj < population_size; jjj++)
				{
					if (adults[jjj].alive)
					{
						if (adults[jjj].maternal[j].loci[jj] == 0)
							allele_freq++;
						if (adults[jjj].paternal[j].loci[jj] == 0)
							allele_freq++;
					}
				}
				allele_freq = allele_freq / (2*population_size);
				output_file << '\t' << allele_freq;
				maf[j].per_locus[jj] = allele_freq;
			}
		}
	}
	void output_summary_info(parameters gp, ofstream & summary_output)
	{
		double dtemp;
		if (gp.parent_trait)
		{
			dtemp = calc_freq_parent(gp);
			summary_output << "\t" << parent_thresh << '\t' << dtemp;
		}
		else
			summary_output << "\tNA\tNA";
		if (gp.court_trait)
		{
			dtemp = calc_freq_courter(gp);
			summary_output << "\t" << courter_thresh << '\t' << dtemp;
		}
		else
			summary_output << "\tNA\tNA";
		output_allele_freqs(gp, summary_output);
	}
	void output_genotypes_vcf(parameters gp, int pop_id)
	{
		int j, jj, jjj;
		string vcf_name = gp.base_name + "_pop_" + to_string(pop_id) + ".vcf";
		ofstream vcf;
		vcf.open(vcf_name);
		//output the header
		vcf << "##fileformat=VCFv4.0";
		string date = determine_date();
		vcf<< "\n##fileDate=" << date;
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
					if(adults[j].courter && adults[j].parent)
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
				vcf << '\n' << j << '\t' << jj << '\t' << j << "." << jj << "\t0\t1\t100\tPASS\tAF=" << maf[j].per_locus[jj] << "\tGT";
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
	void output_trait_info(parameters gp, int pop_id, ofstream & output)
	{
		int j, jj;
		for (j = 0; j < population_size; j++)
		{
			output << '\n' << pop_id << '\t' << j;// "Pop\tIndividual\tSex\tCourter\tParent\tPreference\tMateFound\tPotRS";
			if (adults[j].female)
				output << "\tFEMALE";
			else
				output << "\tMALE";
			output << '\t' << adults[j].courter << '\t' << adults[j].parent << '\t' << adults[j].female_pref << '\t' << adults[j].mate_found << '\t' << adults[j].pot_rs;
		}
	}
};