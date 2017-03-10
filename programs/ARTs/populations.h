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
	vector<double> theta, mean_fem_traits, mean_mal_traits;
	vector<individual> adults;
	vector<individual> progeny;
	vector<tracker> courter_qtls,parent_qtls,pref_qtls,courter_env_qtls,parent_env_qtls,pref_env_qtls, maf, hs;


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

	void initialize(parameters gp)
	{
		cout << "Initializing a population.\n";
		population_size = gp.carrying_capacity;
		if (gp.cor_prefs || gp.ind_pref)
			random_mating = true;
		int j, jj, jjj, k;
		num_fem = num_mal = 0;
		migrant_index = gp.carrying_capacity;
		sex_theta = -100;
		//set up the qtls for the traits
		if (gp.court_trait == true)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				courter_qtls.push_back(tracker());
				for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					courter_qtls[j].per_locus.push_back(randnum(gp.num_markers));
			}
		}
		if (gp.parent_trait == true)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				parent_qtls.push_back(tracker());
				for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					parent_qtls[j].per_locus.push_back(randnum(gp.num_markers));
			}
		}
		if (gp.cor_prefs == true)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				pref_qtls.push_back(tracker());
				for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					pref_qtls[j].per_locus.push_back(courter_qtls[j].per_locus[jj]);
			}
		}
		if (gp.ind_pref == true)
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				pref_qtls.push_back(tracker());
				for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					pref_qtls[j].per_locus.push_back(randnum(gp.num_markers));
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
		for (j = 0; j < (gp.carrying_capacity + max_num_migrants); j++)
		{
			adults.push_back(individual());
			if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
				adults[j].courter_trait = 0;
			if (gp.parent_trait)
				adults[j].parent_trait = 0;
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
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						adults[j].maternal[jj].courter_ae.push_back(double());
						adults[j].paternal[jj].courter_ae.push_back(double());
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						adults[j].maternal[jj].parent_ae.push_back(double());
						adults[j].paternal[jj].parent_ae.push_back(double());
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						adults[j].maternal[jj].pref_ae.push_back(double());
						adults[j].paternal[jj].pref_ae.push_back(double());
					}
				}
			}
		}
		//set up progeny
		for (j = 0; j < (gp.max_fecund*(gp.carrying_capacity + max_num_migrants)); j++)
		{
			progeny.push_back(individual());
			if (gp.court_trait || gp.ind_pref || gp.cor_prefs)
				progeny[j].courter_trait = 0;
			if (gp.parent_trait)
				progeny[j].parent_trait = 0;
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
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						progeny[j].maternal[jj].courter_ae.push_back(double());
						progeny[j].paternal[jj].courter_ae.push_back(double());
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						progeny[j].maternal[jj].parent_ae.push_back(double());
						progeny[j].paternal[jj].parent_ae.push_back(double());
					}
				}
				if (gp.ind_pref || gp.cor_prefs)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						progeny[j].maternal[jj].pref_ae.push_back(double());
						progeny[j].paternal[jj].pref_ae.push_back(double());
					}
				}
			}
		}
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
					for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					{
						adults[j].courter_int.push_back(tracker());
						adults[j].courter_Y.push_back(tracker());						
						if (jj > gp.num_env_qtl)
						{
							adults[j].courter_Z.push_back(1);
							adults[j].courter_x.push_back(double());
						}
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
						{
							adults[j].courter_Y[jj].per_locus.push_back(1);
							adults[j].courter_int[jj].per_locus.push_back(randnorm(0, 0.5));
						}
					}
				}
				if (gp.parent_trait == true)
				{
					for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					{
						adults[j].parent_int.push_back(tracker());
						adults[j].parent_Y.push_back(tracker());
						if (jj > gp.num_env_qtl)
						{
							adults[j].parent_Z.push_back(1);
							adults[j].parent_x.push_back(double());
						}
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
						{
							adults[j].parent_Y[jj].per_locus.push_back(1);//they all start out with an interaction
							adults[j].parent_int[jj].per_locus.push_back(randnorm(0, 0.5));
						}
					}
				}
				if (!random_mating)
				{
					for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					{
						adults[j].pref_int.push_back(tracker());
						adults[j].pref_Y.push_back(tracker());						
						if (jj > gp.num_env_qtl)
						{
							adults[j].pref_Z.push_back(1);
							adults[j].pref_x.push_back(double());
						}
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
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
					for (jjj = 0; jjj < (gp.num_env_qtl+gp.num_qtl); jjj++)
					{
						adults[j].maternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
						adults[j].paternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
					{
						adults[j].maternal[jj].parent_ae[jjj] = tempallele2[j%gp.num_alleles];
						adults[j].paternal[jj].parent_ae[jjj] = tempallele2[j%gp.num_alleles];
					}
				}
				if (!random_mating)
				{
					for (jjj = 0; jjj < (gp.num_env_qtl + gp.num_qtl); jjj++)
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
	void pref_gauss_malemu(parameters gp)
	{
		int j;
		double mean;
		if (gp.court_trait)
			mean = calc_mean_courter_ae(gp);
		else
			mean = calc_mean_parent_ae(gp);
		for (j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
				adults[j].female_pref = mean;
		}
	}
	void pref_gauss_femtrait_noenv(parameters gp)
	{
		for (int j = 0; j < gp.carrying_capacity; j++)
		{
			if (adults[j].alive && adults[j].female)
			{
				adults[j].calc_preference_trait(gp);
			}
		}
	}
	void pref_fd_courter(parameters gp)
	{
		int j, jj;
		double mean, freqT, count;
		mean = freqT = count = 0;
		freqT = calc_freq_courter(gp);
		mean = 0;
		if (freqT < 0.5) // calculate the mean
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].courter && adults[j].alive && !adults[j].female)
				{
					mean = mean + adults[j].courter_trait;
					count++;
				}
			}
			mean = mean / count;
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = mean;
			}
		}
		else
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (!adults[j].courter && adults[j].alive && !adults[j].female)
				{
					mean = mean + adults[j].courter_trait;
					count++;
				}
			}
			mean = mean / count;
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = mean;
			}
		}
	}
	void pref_fd_parent(parameters gp)
	{
		int j, jj;
		double mean, freqT, count;
		mean = freqT = count = 0;
		freqT = calc_freq_parent(gp);
		mean = 0;
		if (freqT < 0.5) // calculate the mean
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].parent && adults[j].alive && !adults[j].female)
				{
					mean = mean + adults[j].parent_trait;
					count++;
				}
			}
			mean = mean / count;
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = mean;
			}
		}
		else
		{
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (!adults[j].parent && adults[j].alive && !adults[j].female)
				{
					mean = adults[j].parent_trait;
					count++;
				}
			}
			mean = mean / count;
			//assign female preferences
			for (j = 0; j < gp.carrying_capacity; j++)
			{
				if (adults[j].alive && adults[j].female)
					adults[j].female_pref = mean;
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
				adults[j].calc_preference_trait(gp);//consider environmental cue...
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
					pref_gauss_malemu(gp);
				else //if it's correlated, then preference should be based on female's own trait
					pref_gauss_femtrait_noenv(gp);
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
				if (gp.court_trait)
				{
					xswitch = gp.max_encounters / (gp.rs_c - gp.rs_nc);
					meanT = meanF = 0;
					count = 0;
					for (jj = 0; jj < gp.carrying_capacity; jj++)
					{
						if (adults[jj].alive && !adults[jj].female)
						{
							if (adults[jj].courter)
								meanT = meanT + adults[jj].courter_trait;
							else
								meanF = meanF + adults[jj].courter_trait;
							count++;
						}
					}
					meanT = meanT / count;
					meanF = meanF / count;
				}
				else
				{
					xswitch = gp.max_encounters / (gp.rs_p - gp.rs_np);
					meanT = meanF = 0;
					count = 0;
					for (jj = 0; jj < gp.carrying_capacity; jj++)
					{
						if (adults[jj].alive && !adults[jj].female)
						{
							if (adults[jj].parent)
								meanT = meanT + adults[jj].parent_trait;
							else
								meanF = meanF + adults[jj].parent_trait;
							count++;
						}
					}
					meanT = meanT / count;
					meanF = meanF / count;
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
				for (RCi = 0; RCi < (gp.num_env_qtl + gp.num_qtl); RCi++)
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
				for (RCi = 0; RCi < (gp.num_env_qtl + gp.num_qtl); RCi++)
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
				for (RCi = 0; RCi < (gp.num_env_qtl + gp.num_qtl); RCi++)
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
				for (RCi = 0; RCi < (gp.num_env_qtl + gp.num_qtl); RCi++)
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
				for (RCi = 0; RCi < (gp.num_env_qtl + gp.num_qtl); RCi++)
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
	void mating(bool output, string out_name, parameters gp)
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

	//selection

};