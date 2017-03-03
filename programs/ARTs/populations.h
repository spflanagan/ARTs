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
		//set up traits and G matrix
		
		//figure out qtls
		if (qtl_list.size() == 0)//if you need to assign qtl
		{
			for (j = 0; j < gp.num_chrom; j++)
			{
				qtl_list.push_back(tracker());
				for (jj = 0; jj < gp.num_qtl; jj++)
					qtl_list[j].per_locus.push_back(randnum(gp.num_markers));//need to index it by chromosome
			}
		}
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
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				adults[j].maternal.push_back(chromosome());
				adults[j].paternal.push_back(chromosome());
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci.push_back(int());
					adults[j].paternal[jj].loci.push_back(int());
				}
				for (jjj = 0; jjj < gp.num_traits; jjj++)
				{
					adults[j].maternal[jj].allelic_effects.mat_rows.push_back(rows());
					adults[j].paternal[jj].allelic_effects.mat_rows.push_back(rows());
					for (k = 0; k < gp.num_qtl; k++)
					{
						adults[j].maternal[jj].allelic_effects.mat_rows[jjj].row.push_back(double());
						adults[j].paternal[jj].allelic_effects.mat_rows[jjj].row.push_back(double());
					}
				}
				adults[j].maternal[jj].allelic_effects.ncol = adults[j].paternal[jj].allelic_effects.ncol = gp.num_qtl;
				adults[j].maternal[jj].allelic_effects.nrow = adults[j].paternal[jj].allelic_effects.nrow = gp.num_traits;
			}
			for (jj = 0; jj < gp.num_traits; jj++)
			{
				adults[j].genotypes.push_back(double());
				adults[j].traits.push_back(double());
			}
		}
		//set up progeny
		for (j = 0; j < (gp.max_fecund*(gp.carrying_capacity + max_num_migrants)); j++)
		{
			progeny.push_back(individual());
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				progeny[j].maternal.push_back(chromosome());
				progeny[j].paternal.push_back(chromosome());
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					progeny[j].maternal[jj].loci.push_back(int());
					progeny[j].paternal[jj].loci.push_back(int());
				}
				for (jjj = 0; jjj < gp.num_traits; jjj++)
				{
					progeny[j].maternal[jj].allelic_effects.mat_rows.push_back(rows());
					progeny[j].paternal[jj].allelic_effects.mat_rows.push_back(rows());
					for (k = 0; k < gp.num_qtl; k++)
					{
						progeny[j].maternal[jj].allelic_effects.mat_rows[jjj].row.push_back(double());
						progeny[j].paternal[jj].allelic_effects.mat_rows[jjj].row.push_back(double());
					}
				}
				progeny[j].maternal[jj].allelic_effects.ncol = progeny[j].paternal[jj].allelic_effects.ncol = gp.num_qtl;
				progeny[j].maternal[jj].allelic_effects.nrow = progeny[j].paternal[jj].allelic_effects.nrow = gp.num_traits;
			}
			for (jj = 0; jj < gp.num_traits; jj++)
			{
				progeny[j].genotypes.push_back(double());
				progeny[j].traits.push_back(double());
			}
		}
		//assign loci and alleles
		vector<double> tempallele1, tempallele2;
		if (!gp.random_g_alleles)
		{
			for (j = 0; j < gp.num_alleles; j++)
			{
				tempallele1.push_back(randnorm(0, gp.allelic_std_dev));
				tempallele2.push_back(randnorm(0, gp.allelic_std_dev));
			}
		}
		for (j = 0; j < (gp.carrying_capacity + max_num_migrants); j++)
		{
			if (j < gp.carrying_capacity)
				adults[j].alive = true;
			else
				adults[j].alive = false;//it's an empty slot for migrants
			for (jj = 0; jj < gp.num_chrom; jj++)
			{
				for (jjj = 0; jjj < gp.num_traits; jjj++)
				{
					for (k = 0; k < gp.num_qtl; k++)
					{
						if (!gp.random_g_alleles)
						{
							adults[j].maternal[jj].allelic_effects.mat_rows[jjj].row[k] = tempallele1[j%gp.num_alleles];
							adults[j].paternal[jj].allelic_effects.mat_rows[jjj].row[k] = tempallele1[j%gp.num_alleles];
						}
						else
						{
							adults[j].maternal[jj].allelic_effects.mat_rows[jjj].row[k] = randnorm(0, gp.allelic_std_dev);
							adults[j].paternal[jj].allelic_effects.mat_rows[jjj].row[k] = randnorm(0, gp.allelic_std_dev);
						}
						adults[j].genotypes[jjj] = adults[j].genotypes[jjj] + adults[j].maternal[jj].allelic_effects.mat_rows[jjj].row[k] + adults[j].paternal[jj].allelic_effects.mat_rows[jjj].row[k];
					}
				}
				int index = 0;
				for (jjj = 0; jjj < gp.num_markers; jjj++)
				{
					adults[j].maternal[jj].loci[jjj] = j%gp.num_alleles;
					adults[j].paternal[jj].loci[jjj] = j%gp.num_alleles;
				}

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
			for (jj = 0; jj < gp.num_traits; jj++)
				adults[j].traits[jj] = adults[j].genotypes[jj] + randnorm(0, gp.env_sd);
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
};