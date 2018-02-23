//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Last Updated: 2 March 2017
//Date Started: 2 March 2017
//Purpose: model alternative reproductive tactics with explicit genetic architectures underlying traits and preferences.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include "random_numbers.h"
#include "populations.h"
#include "classes.h"

using namespace std;

int main(int argc, char*argv[])
{
	//useful variables
	int i, ii, iii, num_eq_tries;
	parameters global_params;
	vector<population> pops;
	bool command_line, run;
	vector<double> courter_freqs, parent_freqs;
	vector<bool> eq_reached;

	//parse parameters
	if (argc == 1)
		command_line = false; //this is also true if running command line with no params
	else
		command_line = true;
	if (command_line)
	{
		run = global_params.parse_parameters(argc, argv);
		if (!run)
		{
			return 0;
		}
		else
			cout << "\nRunning the ARTs model with output to base name " << global_params.base_name << '\n';
	}
	else
	{
		cout << "\nRunning the ARTs model with default parameters.\n";
		global_params.set_defaults();
		//OPTIONAL SET PARAMETERS HERE FOR TESTING
		global_params.court_trait= true;
		global_params.no_genetics = true;
		global_params.carrying_capacity = 1000;
		global_params.num_init_gen = 5;
		global_params.num_exp_gen = 2;
		global_params.base_name = "../../results/testing_court-nogenetics";
		global_params.dependent_params();
		global_params.verbose = true;
	}
	
	//output

	string summary_output_name, trait_output_name, qtlinfo_output_name,popdyn_output_name;
	ofstream summary_output, trait_output, qtlinfo_output, popdyn_output;

	trait_output_name = global_params.base_name + "_traits.txt";

	popdyn_output_name = global_params.base_name + "_popdyn.txt";
	popdyn_output.open(popdyn_output_name);
	popdyn_output << "Generation\tPop\tPopSize\tNumMal\tNumFem\tNumProgeny";

	qtlinfo_output_name = global_params.base_name + "_qtlinfo.txt";
	qtlinfo_output.open(qtlinfo_output_name);
	qtlinfo_output << "Pop";

	summary_output_name = global_params.base_name + "_summary.txt";
	summary_output.open(summary_output_name);
	summary_output << "Generation\tPop\tParentThresh\tParentFreq\tParentW\tNonParentW\tCourterThresh\tCourterFreq\tCourterW\tNonCourterW"
		<< "\tFreqNcNp\tFreqCNp\tFreqNcP\tFreqCP";
	for (i = 0; i < global_params.num_chrom; i++)
	{
		for (ii = 0; ii < global_params.num_markers; ii++)
			summary_output << "\tMarker" << i << "." << ii;
	}		

	//Initialize
	cout << "\nInitializing " << global_params.num_pops << " populations.\n";
	for (i = 0; i < global_params.num_pops; i++)
	{
		pops.push_back(population());
		pops[i].initialize(global_params);
		run = pops[i].sanity_checks(global_params);
		if (!run)
		{
			if(command_line)
				return 0;
			else
			{
				cout << "\nInput integer to close dialogue\n";
				cin >> ii;
				return ii;
			}
		}
		courter_freqs.push_back(0);
		parent_freqs.push_back(0);
		eq_reached.push_back(false);
		//output QTL info
		if (i == 0)//if it's the first/only pop, write the header
		{
			if (global_params.gene_network)
				pops[i].output_qtl_info(global_params, qtlinfo_output, true,
					pops[i].courter_env_qtls,pops[i].parent_env_qtls,pops[i].pref_env_qtls,
					pops[i].cthresh_env_qtls, pops[i].pthresh_env_qtls);
			else
				pops[i].output_qtl_info(global_params, qtlinfo_output, true);
		}
		qtlinfo_output << "Pop" << i;
		pops[i].output_qtl_info(global_params, qtlinfo_output,false);
	}
	global_params.output_parameters();
	qtlinfo_output.close();
	//Run the program
	cout << "\nRunning " << global_params.num_init_gen << " initial generations\n";
	for (i = 0; i < global_params.num_init_gen; i++)
	{
		for (ii = 0; ii < global_params.num_pops; ii++)
		{
			if(pops[ii].population_size > 0)
			{
				if (global_params.verbose)
					cout << i << std::flush;
				else
				{
					if (i % 1000 == 0)
						cout << "\nInitial generation " << i + 1 << " beginning.";
				}
				pops[ii].determine_pop_size(global_params);
				if (global_params.verbose)
					cout << ", " << pops[ii].population_size << " adults" << flush;
				//mating (includes assigning preferences, recombination, and mutation)
				bool write_to_file = false;
				string temp_file_name;
				/*if (i % 1000 == 0)
				{
					write_to_file = true;
					stringstream temp_name;
					temp_name << global_params.base_name << "genotypes" << i << ".txt";
					temp_file_name = temp_name.str();
				}
				else
					write_to_file = false;*/
				pops[ii].nest_and_fertilize(global_params, write_to_file, temp_file_name);

				//viability selection
				pops[ii].viability_selection(global_params);

				//output summary stats
				summary_output << "\n" << i << "\tPop" << ii;
				pops[ii].output_summary_info(global_params, summary_output);//includes allele freqs and RS

				//stochastic survival
				//pops[ii].density_regulation(global_params);
				pops[ii].regulate_popsize(global_params);
				//track frequencies
				double new_parent, new_courter;
				if (global_params.parent_trait)
				{
					new_parent = pops[ii].calc_freq_parent(global_params);
					if (i > 0)
						pops[ii].d_parentfreq.push_back((new_parent - parent_freqs[ii]));
					parent_freqs[ii] = new_parent;
					if (global_params.verbose)
						cout << ", " << new_parent << " parents" << std::flush;
				}
				if (global_params.court_trait)
				{
					new_courter = pops[ii].calc_freq_courter(global_params);
					if (i > 0)
						pops[ii].d_courterfreq.push_back((new_courter - courter_freqs[ii]));
					courter_freqs[ii] = new_courter;

					if (global_params.verbose)
						cout << ", " << new_courter << " courters" << std::flush;
				}
			}
			else
			{
				cout << "\nPopulation has crashed at initial generation " << i << '\n';
				if (command_line)
					return 0;
				else
				{
					cout << "\nInput integer to close dialogue.\n";
					cin >> iii;
					return 0;
				}
			}
			//output population info
			popdyn_output << '\n' << i << '\t' << ii << '\t' << pops[ii].population_size << '\t' << pops[ii].num_mal << '\t' << pops[ii].num_fem << '\t' << pops[ii].num_progeny;
			if(global_params.verbose)
				cout << endl;
		}	
	}
	//calc variance in change in frequencies
	if (global_params.parent_trait)
	{
		for (i = 0; i < global_params.num_pops; i++)
		{
			double mean = 0;
			for (ii = 0; ii < pops[i].d_parentfreq.size(); ii++)
				mean = mean + pops[i].d_parentfreq[ii];
			mean = mean / pops[i].d_parentfreq.size();
			parent_freqs[i] = 0;//use this for the variance
			for (ii = 0; ii < pops[i].d_parentfreq.size(); ii++)
				parent_freqs[i] = parent_freqs[i] + (pops[i].d_parentfreq[ii] - mean)*(pops[i].d_parentfreq[ii] - mean);
			parent_freqs[i] = parent_freqs[i] / pops[i].d_parentfreq.size();
		}
	}
	if (global_params.court_trait)
	{
		for (i = 0; i < global_params.num_pops; i++)
		{
			double mean = 0;
			for (ii = 0; ii < pops[i].d_courterfreq.size(); ii++)
				mean = mean + pops[i].d_courterfreq[ii];
			mean = mean / pops[i].d_courterfreq.size();
			courter_freqs[i] = 0;//use this for the variance
			for (ii = 0; ii < pops[i].d_courterfreq.size(); ii++)
				courter_freqs[i] = courter_freqs[i] + (pops[i].d_courterfreq[ii] - mean)*(pops[i].d_courterfreq[ii] - mean);
			courter_freqs[i] = courter_freqs[i] / pops[i].d_courterfreq.size();
		}
	}

	cout << "\nEvaluating equilibrium";
	num_eq_tries = 0;
	trait_output.open(trait_output_name);
	trait_output << "Pop\tIndividual\tSex\tCourter\tCourtTrait\tParent\tParentTrait\tPreference\tPrefTrait\tMateFound\tPotRS\tLifetimeRS\tAlive";
	while (num_eq_tries < global_params.num_exp_gen)
	{
		for (i = 0; i < global_params.num_pops; i++)
		{
			if (pops[i].population_size > 0)
			{
				//is it at equilibrium?
				if (!eq_reached[i])
				{
					if (global_params.parent_trait && global_params.court_trait)
					{
						if (pops[i].d_parentfreq.back() <= parent_freqs[i] && pops[i].d_courterfreq.back() <= courter_freqs[i])
							eq_reached[i] = true;
					}
					if (global_params.parent_trait && !global_params.court_trait)
					{
						if (pops[i].d_parentfreq.back() <= parent_freqs[i])
							eq_reached[i] = true;
					}
					if (!global_params.parent_trait && global_params.court_trait)
					{
						if (pops[i].d_courterfreq.back() <= courter_freqs[i])
							eq_reached[i] = true;
					}
					if (eq_reached[i] == false)//then we try again
					{
						if (global_params.verbose)
							cout << '\n' << i << std::flush;
						else
						{
							if (num_eq_tries % 1000 == 0)
								cout << "\nEquilibrium generation " << num_eq_tries + 1 << " beginning.";
						}
						pops[i].determine_pop_size(global_params);						
						//mating (includes assiging preferences, recombination, and mutation)
						pops[i].nest_and_fertilize(global_params, false, "temp");
						//selection
						pops[i].viability_selection(global_params);
						//output summary stats
						summary_output << "\n" << global_params.num_init_gen + num_eq_tries << "\tPop" << i;
						pops[i].output_summary_info(global_params, summary_output);//includes allele freqs
						//stochastic survival
						//pops[i].density_regulation(global_params);
						pops[i].regulate_popsize(global_params);
						//track frequencies
						double new_parent, new_courter;
						if (global_params.parent_trait)
						{
							new_parent = pops[i].calc_freq_parent(global_params);
							if (i > 0)
								pops[i].d_parentfreq.push_back((new_parent - parent_freqs[i]));
						}
						if (global_params.court_trait)
						{
							new_courter = pops[i].calc_freq_courter(global_params);
							if (i > 0)
								pops[i].d_courterfreq.push_back((new_courter - courter_freqs[i]));
						}
					}
					else//output some stuff
					{
						cout << "\nEquilibrium reached for population " << i << " at generation " << global_params.num_init_gen + num_eq_tries;
						pops[i].output_genotypes_vcf(global_params, i);
						pops[i].output_trait_info(global_params, i, trait_output);
						num_eq_tries = global_params.num_exp_gen;
					}
				}
				else
				{
					cout << "\nPopulation has crashed at experimental generation " << num_eq_tries << '\n';
					if (command_line)
						return 0;
					else
					{
						cout << "\nInput integer to close dialogue.\n";
						cin >> ii;
						return 0;
					}
				}
				//output population info
				popdyn_output << '\n' << global_params.num_init_gen + num_eq_tries << '\t' << i << '\t' 
					<< pops[i].population_size << '\t' << pops[i].num_mal << '\t' << pops[i].num_fem << '\t' << pops[i].num_progeny;
				num_eq_tries++;
			}
		}
	}
	//check to see if they reached equilibrium
	for (i = 0; i < global_params.num_pops; i++)
	{
		if (!eq_reached[i])
		{
			cout << "\nNo equilibrium could be reached for population " << i << " with population size " << pops[i].population_size;
			pops[i].output_genotypes_vcf(global_params, i);
			pops[i].output_trait_info(global_params, i, trait_output);
		}
	}
	
	
	//output
	summary_output.close();
	trait_output.close();
	popdyn_output.close();
	cout << "\nDone!\n";
	if (command_line)
		return 0;
	else
	{
		cout << "\nInput integer to quit.\n";
		cin >> i;
		return 0;
	}
}