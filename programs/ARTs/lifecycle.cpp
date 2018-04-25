//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Last Updated: 12 March 2018
//Date Started: 2 March 2017
//Purpose: model alternative reproductive tactics with explicit genetic architectures underlying traits and preferences.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <chrono>
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
			std::cout << "\nRunning the ARTs model with output to base name " << global_params.base_name << '\n';
	}
	else
	{
		std::cout << "\nRunning the ARTs model with default parameters.\n";
		global_params.set_defaults();
		//OPTIONAL SET PARAMETERS HERE FOR TESTING
		global_params.parent_trait= false;
		global_params.court_trait = true;
		global_params.no_genetics = true;
		global_params.num_init_gen = 5;
		global_params.num_exp_gen = 2;
		global_params.base_name = "../../results/testing";
		global_params.dependent_params();
		global_params.verbose = false;
		global_params.optimize = true;
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
		<< "\tFreqNcNp\tFreqCNp\tFreqNcP\tFreqCP\tPrefThresh\tPrefFreq";
	for (i = 0; i < global_params.num_chrom; i++)
	{
		for (ii = 0; ii < global_params.num_markers; ii++)
			summary_output << "\tMarker" << i << "." << ii;
	}		

	//Initialize
	std::cout << "\nInitializing " << global_params.num_pops << " populations.\n" << std::flush;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	//start with the same base population
	if (global_params.same_base)
	{
		pops.push_back(population());
		pops[0].initialize(global_params);
		run = pops[0].sanity_checks(global_params);
		if (!run)
		{
			summary_output.close();
			popdyn_output.close();
			if (command_line)
				return 0;
			else
			{
				std::cout << "\nInput integer to close dialogue\n";
				cin >> ii;
				return ii;
			}
		}
		for (i = 1; i < global_params.num_pops; i++)
		{
			pops.push_back(population());
			pops[i] = pops[0];
			run = pops[i].sanity_checks(global_params);
			if (!run)
			{
				summary_output.close();
				popdyn_output.close();
				if (command_line)
					return 0;
				else
				{
					std::cout << "\nInput integer to close dialogue\n";
					cin >> ii;
					return ii;
				}
			}
		}
	}
	else
	{
		for (i = 0; i < global_params.num_pops; i++)
		{
			pops.push_back(population());
			pops[i].initialize(global_params);
			run = pops[i].sanity_checks(global_params);
			if (!run)
			{
				summary_output.close();
				popdyn_output.close();
				if (command_line)
					return 0;
				else
				{
					std::cout << "\nInput integer to close dialogue\n";
					cin >> ii;
					return ii;
				}
			}
			
		}
	}
	for (i = 0; i < global_params.num_pops; i++)
	{
		//set up equilibrium and frequency trackers
		courter_freqs.push_back(0);
		parent_freqs.push_back(0);
		if (global_params.parent_trait)
		{
			pops[i].d_parentfreq.push_back(pops[i].calc_freq_parent(global_params));
			if (global_params.verbose)
				std::cout << "\n\t" << pops[i].d_parentfreq[0] << " parents" << std::flush;
		}
		if (global_params.court_trait)
		{
			pops[i].d_courterfreq.push_back(pops[i].calc_freq_courter(global_params));
			if (global_params.verbose)
				std::cout << "\n\t" << pops[i].d_courterfreq[0] << " courters" << std::flush;
		}
		eq_reached.push_back(false);
		//output QTL info
		if (i == 0)//if it's the first/only pop, write the header
		{
			if (global_params.gene_network)
				pops[i].output_qtl_info(global_params, qtlinfo_output, true,
					pops[i].courter_env_qtls, pops[i].parent_env_qtls, pops[i].pref_env_qtls,
					pops[i].cthresh_env_qtls, pops[i].pthresh_env_qtls);
			else
				pops[i].output_qtl_info(global_params, qtlinfo_output, true);
		}
		qtlinfo_output << "Pop" << i;
		pops[i].output_qtl_info(global_params, qtlinfo_output, false);
	}
		
	
	global_params.output_parameters();
	qtlinfo_output.close();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	if (global_params.optimize)
		std::cout << "\n   Initialize took " << duration << " seconds.";
	
	//Run the program
	std::cout << "\nRunning " << global_params.num_init_gen << " initial generations\n" << std::flush;
	for (i = 0; i < global_params.num_init_gen; i++)
	{
		for (ii = 0; ii < global_params.num_pops; ii++)
		{
			if(pops[ii].population_size > 0)
			{
				if (global_params.verbose)
					std::cout << i << std::flush;
				else
				{
					if (i % 1000 == 0)
						std::cout << "\nInitial generation " << i + 1 << " beginning." << std::flush;
				}
				pops[ii].determine_pop_size(global_params);
				if (global_params.verbose)
					std::cout << ", " << pops[ii].population_size << " adults" << flush;
				//mating (includes assigning preferences, recombination, and mutation)
				bool write_to_file = false;
				string temp_file_name;
				t1 = std::chrono::high_resolution_clock::now();
				pops[ii].nest_and_fertilize(global_params, write_to_file, temp_file_name);
				t2 = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				if (global_params.optimize)
					std::cout << "\n   nest_and_fertilize took " << duration << " seconds.";
				//viability selection
				t1 = std::chrono::high_resolution_clock::now();
				pops[ii].viability_selection(global_params);
				t2 = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				if (global_params.optimize)
					std::cout << "\n   viability_selection took " << duration << " seconds.";
				//output summary stats
				summary_output << "\n" << i << "\tPop" << ii;
				pops[ii].output_summary_info(global_params, summary_output);//includes allele freqs and RS

				//stochastic survival
				//pops[ii].density_regulation(global_params);
				t1 = std::chrono::high_resolution_clock::now();
				pops[ii].regulate_popsize(global_params);
				t2 = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				if (global_params.optimize)
					std::cout << "\n   regulate_popsize took " << duration << " seconds.";
				//output summary stats
				summary_output << "\n" << i << "\tPop" << ii;
			}
			else
			{
				std::cout << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
				summary_output.close();
				popdyn_output.close();
				if (command_line)
					return 0;
				else
				{
					std::cout << "\nInput integer to close dialogue.\n" << std::flush;
					cin >> iii;
					return 0;
				}
			}
			//output population info
			popdyn_output << '\n' << i << '\t' << ii << '\t' << pops[ii].population_size << '\t' << pops[ii].num_mal << '\t' << pops[ii].num_fem << '\t' << pops[ii].num_progeny << std::flush;
			if(global_params.verbose)
				std::cout << endl;
		}	
	}
	
	trait_output.open(trait_output_name);
	trait_output << "Pop\tIndividual\tSex\tCourter\tCourtTrait\tParent\tParentTrait\tPreference\tPrefTrait\tMateFound\tPotRS\tLifetimeRS\tAlive";
	//run the last 2000 generations
	for(i = 0; i < global_params.num_exp_gen; i++)
	{
		t1 = std::chrono::high_resolution_clock::now();
		for (ii = 0; ii < global_params.num_pops; ii++)
		{
			if (pops[ii].population_size > 0)
			{
				if (global_params.verbose)
					std::cout << '\n' << global_params.num_init_gen + i << std::flush;
					
				pops[ii].determine_pop_size(global_params);						
				//mating (includes assiging preferences, recombination, and mutation)
				pops[ii].nest_and_fertilize(global_params, false, "temp");
				//selection
				pops[ii].viability_selection(global_params);
				//output summary stats
				summary_output << "\n" << global_params.num_init_gen + i << "\tPop" << ii;
				pops[ii].output_summary_info(global_params, summary_output);//includes allele freqs
				//stochastic survival
				//pops[i].density_regulation(global_params);
				pops[ii].regulate_popsize(global_params);
				//track frequencies
				double new_parent, new_courter;
				if (global_params.parent_trait)
				{
					new_parent = pops[ii].calc_freq_parent(global_params);
					pops[ii].d_parentfreq.push_back((new_parent - parent_freqs[ii]));
					parent_freqs[ii] = new_parent;
					if (global_params.verbose)
						std::cout << ", " << new_parent << " parents" << std::flush;
				}
				if (global_params.court_trait)
				{
					new_courter = pops[ii].calc_freq_courter(global_params);
					pops[ii].d_courterfreq.push_back((new_courter - courter_freqs[ii]));
					courter_freqs[ii] = new_courter;
					if (global_params.verbose)
						std::cout << ", " << new_courter << " courters" << std::flush;
				}
				//output population info
				popdyn_output << '\n' << global_params.num_init_gen + i << '\t' << ii << '\t' 
					<< pops[ii].population_size << '\t' << pops[ii].num_mal << '\t' << pops[ii].num_fem << '\t' << pops[ii].num_progeny << std::flush;
				
			}//if popsize > 0
			else
			{
				std::cout << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
				summary_output.close();
				trait_output.close();
				popdyn_output.close();
				if (command_line)
					return 0;
				else
				{
					std::cout << "\nInput integer to close dialogue.\n" << std::flush;
					cin >> ii;
					return 0;
				}
			}
		}
		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
		if (global_params.optimize)
			std::cout << "\n   experimental generation " << i << " took " << duration << " seconds.";
	}
	
	//Evaluate stasis/equilibrium
	std::cout << "\nEvaluating equilibrium" << std::flush;
	for (i = 0; i < global_params.num_pops; i++)
	{
		//calc variance in change in frequencies
		if (global_params.parent_trait)
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
		if (global_params.court_trait)
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
		//How do the changes over time compare to the variance?
		if (global_params.parent_trait && global_params.court_trait)
		{
			if  (pops[i].d_parentfreq.back() <= parent_freqs[i] && pops[i].d_courterfreq.back() <= courter_freqs[i])
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
		//check to see if they reached equilibrium
		if (!eq_reached[i])
		{
			std::cout << "\nNo equilibrium could be reached for population " << i << " with population size " << pops[i].population_size << std::flush;
			pops[i].output_genotypes_vcf(global_params, i);
			pops[i].output_trait_info(global_params, i, trait_output);
		}
	}
	
	//close output files
	summary_output.close();
	trait_output.close();
	popdyn_output.close();
	std::cout << "\nDone!\n" << std::flush;
	if (command_line)
		return 0;
	else
	{
		std::cout << "\nInput integer to quit.\n" << std::flush;
		cin >> i;
		return 0;
	}
}