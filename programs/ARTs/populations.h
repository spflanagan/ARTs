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
	double sex_theta, sex_omega;
	bool extinct, sexual_selection;
	vector<double> theta, mean_fem_traits, mean_mal_traits;
	vector<individual> adults;
	vector<individual> progeny;
	vector<tracker> courter_qtls,parent_qtls,pref_qtls,courter_env_qtls,parent_env_qtls,pref_env_qtls, maf, hs;


	population()
	{
		population_size = num_mal = num_fem = num_progeny = sex_trait = max_num_migrants = migrant_index = int();
		theta = mean_fem_traits = mean_mal_traits = vector<double>();
		sex_theta = sex_omega = double();
		adults = progeny = vector<individual>();
		courter_env_qtls = courter_qtls = parent_env_qtls = parent_qtls = pref_env_qtls = pref_qtls = maf = hs = vector<tracker>();
		sexual_selection = false;
		extinct = false;
		sgenrand(time(0));
	}

	void initialize(parameters gp)
	{
		cout << "Initializing a population.\n";
		population_size = gp.carrying_capacity;
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
				if (gp.court_trait || gp.ind_pref || gp.cor_prefs)//same container for males and females
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						adults[j].maternal[jj].courter_ae.push_back(double());
						adults[j].paternal[jj].courter_ae.push_back(double());
					}
					if (gp.env_effects)
					{
						for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
						{
							adults[j].maternal[jj].cenv_ae.push_back(double());
							adults[j].paternal[jj].cenv_ae.push_back(double());
						}
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						adults[j].maternal[jj].parent_ae.push_back(double());
						adults[j].paternal[jj].parent_ae.push_back(double());
					}
					if (gp.env_effects)
					{
						for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
						{
							adults[j].maternal[jj].penv_ae.push_back(double());
							adults[j].paternal[jj].penv_ae.push_back(double());
						}
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
				if (gp.court_trait || gp.ind_pref || gp.cor_prefs)//same container for males and females
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						progeny[j].maternal[jj].courter_ae.push_back(double());
						progeny[j].paternal[jj].courter_ae.push_back(double());
					}
					if (gp.env_effects)
					{
						for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
						{
							progeny[j].maternal[jj].cenv_ae.push_back(double());
							progeny[j].paternal[jj].cenv_ae.push_back(double());
						}
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						progeny[j].maternal[jj].parent_ae.push_back(double());
						progeny[j].paternal[jj].parent_ae.push_back(double());
					}
					if (gp.env_effects)
					{
						for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
						{
							progeny[j].maternal[jj].penv_ae.push_back(double());
							progeny[j].paternal[jj].penv_ae.push_back(double());
						}
					}
				}
			}
		}
		//assign loci and alleles
		vector<double> tempallele1, tempallele2;
		for (j = 0; j < gp.num_alleles; j++) //set up specific alleles.
		{
			tempallele1.push_back(randnorm(0, gp.allelic_std_dev));
			tempallele2.push_back(randnorm(0, gp.allelic_std_dev));
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
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
						{
							if (rand() <= 0.2)//it doesn't have an interaction
							{
								adults[j].courter_int[jj].per_locus.push_back(0);
							}
							else//it has a weight from a normal distribution with mean 0 and sd 0.5
							{
								adults[j].courter_int[jj].per_locus.push_back(randnorm(0, 0.5));
							}
						}
					}
				}
				if (gp.parent_trait == true)
				{
					for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					{
						adults[j].parent_int.push_back(tracker());
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
						{
							if (rand() <= 0.2)//it doesn't have an interaction
							{
								adults[j].parent_int[jj].per_locus.push_back(0);
							}
							else//it has a weight from a normal distribution with mean 0 and sd 0.5
							{
								adults[j].parent_int[jj].per_locus.push_back(randnorm(0, 0.5));
							}
						}
					}
				}
				if (gp.cor_prefs || gp.ind_pref)
				{
					for (jj = 0; jj < (gp.num_qtl + gp.num_env_qtl); jj++)
					{
						adults[j].pref_int.push_back(tracker());
						for (jjj = 0; jjj < (gp.num_qtl + gp.num_env_qtl); jjj++)
						{
							if (rand() <= 0.2)//it doesn't have an interaction
							{
								adults[j].pref_int[jj].per_locus.push_back(0);
							}
							else//it has a weight from a normal distribution with mean 0 and sd 0.5
							{
								adults[j].pref_int[jj].per_locus.push_back(randnorm(0, 0.5));
							}
						}
					}
				}
			}
			
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				if (gp.court_trait || gp.cor_prefs || gp.ind_pref)
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						adults[j].maternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
						adults[j].paternal[jj].courter_ae[jjj] = tempallele1[j%gp.num_alleles];
					}
					for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
					{
						adults[j].maternal[jj].cenv_ae[jjj] = tempallele2[j%gp.num_alleles];
						adults[j].paternal[jj].cenv_ae[jjj] = tempallele2[j%gp.num_alleles];
					}
				}
				if (gp.parent_trait)
				{
					for (jjj = 0; jjj < gp.num_qtl; jjj++)
					{
						adults[j].maternal[jj].parent_ae[jjj] = tempallele1[j%gp.num_alleles];
						adults[j].paternal[jj].parent_ae[jjj] = tempallele1[j%gp.num_alleles];
					}
					for (jjj = 0; jjj < gp.num_env_qtl; jjj++)
					{
						adults[j].maternal[jj].penv_ae[jjj] = tempallele2[j%gp.num_alleles];
						adults[j].paternal[jj].penv_ae[jjj] = tempallele2[j%gp.num_alleles];
					}
				}
				int index = 0;
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci[jjj] = j%gp.num_alleles;
					adults[j].paternal[jj].loci[jjj] = j%gp.num_alleles;
				}
			}
			if (gp.court_trait || gp.cor_prefs || gp.ind_pref)
			{
				//to initialize env cue = 0
				adults[j].calc_courter_trait(0, gp);
			}
			if (gp.parent_trait)
			{
				adults[j].calc_parent_trait(0, gp);
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

	bool recombine_chromosome(chromosome& chrom, individual &parent, double expected_recomb_events, int which_chrom, tracker& qtls, parameters gp)
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
				for (RCi = 0; RCi < gp.num_qtl; RCi++)
				{
					for (RCj = 0; RCj < num_segments; RCj++)
					{

					}
				}
			}
			for (RCk = 0; RCk < num_traits; RCk++)
			{
				for (RCi = 0; RCi < num_qtl; RCi++)
				{
					for (RCj = 0; RCj < num_segments; RCj++)
					{
						if (qtls.per_locus[RCi] >= segment_start[RCj] && qtls.per_locus[RCi] >= segment_end[RCj])
						{
							if (segment_maternal[RCj])
							{
								chrom.allelic_effects.mat_rows[RCk].row[RCi] = parent.maternal[which_chrom].allelic_effects.mat_rows[RCk].row[RCi];
							}
							else
							{
								chrom.allelic_effects.mat_rows[RCk].row[RCi] = parent.paternal[which_chrom].allelic_effects.mat_rows[RCk].row[RCi];
							}
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
			{
				for (RCi = 0; RCi < num_markers; RCi++)
					chrom.loci[RCi] = parent.maternal[which_chrom].loci[RCi];
				for (RCi = 0; RCi < num_traits; RCi++)
				{
					for (RCk = 0; RCk < num_qtl; RCk++)
					{
						chrom.allelic_effects.mat_rows[RCi].row[RCk] = parent.maternal[which_chrom].allelic_effects.mat_rows[RCi].row[RCk];
						chrom.allelic_effects.mat_rows[RCi].row[RCk] = parent.maternal[which_chrom].allelic_effects.mat_rows[RCi].row[RCk];
					}
				}
			}
			else
			{
				for (RCi = 0; RCi < num_markers; RCi++)
					chrom.loci[RCi] = parent.paternal[which_chrom].loci[RCi];
				for (RCi = 0; RCi < num_traits; RCi++)
				{
					for (RCk = 0; RCk < num_qtl; RCk++)
					{
						chrom.allelic_effects.mat_rows[RCi].row[RCk] = parent.paternal[which_chrom].allelic_effects.mat_rows[RCi].row[RCk];
						chrom.allelic_effects.mat_rows[RCi].row[RCk] = parent.paternal[which_chrom].allelic_effects.mat_rows[RCi].row[RCk];
					}
				}

			}
			chrom.allelic_effects.ncol = num_qtl;
			chrom.allelic_effects.nrow = num_traits;
		}
		return recombined;
	}
};