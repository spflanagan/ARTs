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
	int i, ii, iii, num_eq_tries, crash_counter;
	parameters global_params;
	vector<population> pops;
	bool command_line, run;
	vector<double> courter_freqs, parent_freqs;
	vector<bool> eq_reached;

	//set up log file
	string log_name;
	ofstream log_out;
	
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
		{
			if (global_params.log_file)
            {
                log_name = global_params.base_name + ".log";
                log_out.open(log_name);
                log_out << "\nRunning the ARTs model (from the command line) with output to base name " << global_params.base_name << '\n' <<std::flush;
            }
			std::cout << "\nRunning the ARTs model (from the command line) with output to base name " << global_params.base_name << '\n'<<std::flush;
		}
	}
	else
	{
		global_params.log_file = false;
		std::cout << "\nRunning the ARTs model with default parameters.\n";
		global_params.set_defaults();
		//OPTIONAL SET PARAMETERS HERE FOR TESTING
		global_params.parent_trait= false;
		global_params.court_trait = true;
		global_params.supergene = true;
		global_params.supergene_prop = 0.05;
		global_params.num_init_gen = 5;
		global_params.num_exp_gen = 2;
		global_params.base_name = "../../results/courter_supergene_propTest";
		global_params.dependent_params();
		global_params.verbose = false;
		global_params.optimize = true;
	}
	
	//output

	string summary_output_name, trait_output_name, qtlinfo_output_name,markers_output_name, ae_output_name;
	ofstream summary_output, trait_output, qtlinfo_output, markers_output, ae_output;
	

	trait_output_name = global_params.base_name + "_traits.txt";
    trait_output.open(trait_output_name);
    trait_output << "Gen\tPop\tIndividual\tSex\tCourter\tCourtTrait\tParent\tParentTrait\tPreference\tPrefTrait\tMateFound\tPotRS\tLifetimeRS\tAlive";
    
	markers_output_name = global_params.base_name + "_markers.txt";
	markers_output.open(markers_output_name);
	markers_output << "Generation\tPop";
    for (i = 0; i < global_params.num_chrom; i++)
    {
        for (ii = 0; ii < global_params.num_markers; ii++)
            markers_output << "\tMarker" << i << "." << ii;
    }

	qtlinfo_output_name = global_params.base_name + "_qtlinfo.txt";
	qtlinfo_output.open(qtlinfo_output_name);
	qtlinfo_output << "Pop";

	ae_output_name = global_params.base_name + "_allelic-effects.txt";
	ae_output.open(ae_output_name);
	ae_output << "Gen\tPop";

	summary_output_name = global_params.base_name + "_summary.txt";
	summary_output.open(summary_output_name);
	summary_output << "Generation\tPop\tPopSize\tNumMal\tNumFem\tNumProgeny\tParentThresh\tParentFreq\tParentAEmean\tParentAEsd\tParentW\tNonParentW\tCourterThresh\tCourterFreq\tCourterAEmean\tCourterAEsd\tCourterW\tNonCourterW"
		<< "\tFreqNcNp\tFreqCNp\tFreqNcP\tFreqCP\tPrefThresh\tPrefFreq";
	

	//Initialize
	if(global_params.log_file)
		log_out << "\nInitializing " << global_params.num_pops << " populations.\n";
	else
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
			markers_output.close();
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
				markers_output.close();
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
				markers_output.close();
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
			{
				if (global_params.log_file)
					log_out<< "\n\t" << pops[i].d_parentfreq[0] << " parents";
				else
					std::cout << "\n\t" << pops[i].d_parentfreq[0] << " parents" << std::flush;
			}
		}
		if (global_params.court_trait)
		{
			pops[i].d_courterfreq.push_back(pops[i].calc_freq_courter(global_params));
			if (global_params.verbose)
			{
				if (global_params.log_file)
					log_out << "\n\t" << pops[i].d_courterfreq[0] << " courters" ;
				else
					std::cout << "\n\t" << pops[i].d_courterfreq[0] << " courters" << std::flush;
			}				
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
			pops[i].output_allelic_effects(global_params,ae_output, true);
		}
		qtlinfo_output << "Pop" << i;
		pops[i].output_qtl_info(global_params, qtlinfo_output, false);
		ae_output << '\n' << 0 << "\tPop" << i;
		pops[i].output_allelic_effects(global_params,ae_output, false);
	}
		
	
	global_params.output_parameters();
	qtlinfo_output.close();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	if (global_params.optimize)
	{
		if (global_params.log_file)
			log_out << "\n   Initialize took " << duration << " seconds.";
		else
			std::cout << "\n   Initialize took " << duration << " seconds." << std::flush;

	}
	
	
	//Run the program
	if (global_params.log_file)
		log_out <<  "\nRunning " << global_params.num_init_gen << " initial generations\n" ;
	else
		std::cout << "\nRunning " << global_params.num_init_gen << " initial generations\n" << std::flush;
	for (i = 0; i < global_params.num_init_gen; i++)
	{
		crash_counter = 0;
		for (ii = 0; ii < global_params.num_pops; ii++)
		{
			pops[ii].determine_pop_size(global_params);
			if(pops[ii].population_size > 0)
			{
				if(i == 0 && global_params.ae_vcf)
				{
					pops[ii].output_ae_vcf(global_params, ii);
				}
				
				if (global_params.verbose)
				{
					if (global_params.log_file)
						log_out << i;
					else
						std::cout << i << std::flush;
				}					
				else
				{
					if (i % 1000 == 0)
					{
						if(global_params.log_file)
							log_out << "\nInitial generation " << i + 1 << " beginning.";
						std::cout << "\nInitial generation " << i + 1 << " beginning." << std::flush;
					}
				}
				if (global_params.verbose)
				{
					if (global_params.log_file)
                        log_out<< ", " << pops[ii].population_size << " adults; ";
					else
                        std::cout << ", " << pops[ii].population_size << " adults; " << std::flush;
				}
				//mating (includes assigning preferences, recombination, and mutation)
				bool write_to_file = false;
				string temp_file_name;
				t1 = std::chrono::high_resolution_clock::now();
				pops[ii].nest_and_fertilize(global_params, write_to_file, temp_file_name);
				//output trait values for gen 0 to include mating success
				if(i == 0) pops[ii].output_trait_info(global_params, i, ii, trait_output, log_out);
				t2 = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				if (global_params.optimize)
				{
					if (global_params.log_file)
						log_out << "\n   nest_and_fertilize took " << duration << " seconds.";
					else
						std::cout << "\n   nest_and_fertilize took " << duration << " seconds." << std::flush;
				}
				//viability selection
                if(global_params.viability_selection)
                {
                    t1 = std::chrono::high_resolution_clock::now();
                    pops[ii].viability_selection(global_params);
                    t2 = std::chrono::high_resolution_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
                    if (global_params.optimize)
                    {
                        if (global_params.log_file)
                            log_out<< "\n   viability_selection took " << duration << " seconds.";
                        else
                            std::cout << "\n   viability_selection took " << duration << " seconds." << std::flush;
                    }
                }
				//output summary stats
				summary_output << "\n" << i << "\tPop" << ii;
				pops[ii].output_summary_info(global_params, summary_output, log_out);//includes RS
                //output allele frequencies
                markers_output << "\n" << i << "\tPop" << ii;
                pops[ii].output_allele_freqs(global_params, markers_output);
				ae_output << '\n' << i+1 << "\tPop" << ii;
				pops[ii].output_allelic_effects(global_params,ae_output, false);
				//stochastic survival
				//pops[ii].density_regulation(global_params);
				t1 = std::chrono::high_resolution_clock::now();
				pops[ii].regulate_popsize(global_params);
				t2 = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				if (global_params.optimize)
				{
					if (global_params.log_file)
						log_out << "\n   regulate_popsize took " << duration << " seconds.";
					else
						std::cout << "\n   regulate_popsize took " << duration << " seconds." << std::flush;
				}
			}
			else //if the population size is 0 the population has crashed
			{
				if(!pops[ii].extinct)
				{				
					if (global_params.log_file)
						log_out << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
					else
						std::cout << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
					crash_counter++;
				}
				pops[ii].extinct = true;
			}
			
			
			if(global_params.verbose && !global_params.log_file)
				std::cout << std::flush;
		}
		// if all of the populations have crashed then we can stop.	
		if (crash_counter == global_params.num_pops+1)
		{
			if (global_params.log_file)
				log_out << "\n" << global_params.base_name << ": All populations have crashed at experimental generation " << i << '\n' << std::flush;
			else
				std::cout << "\n" << global_params.base_name << ": All populations have crashed at experimental generation " << i << '\n' << std::flush;
			summary_output.close();
			markers_output.close();
			ae_output.close();
			if (!command_line)
				return 0;
			else
			{
				std::cout << "\nInput integer to close dialogue.\n" << std::flush;
				cin >> iii;
				return 0;
			}
		}
	}
	
	//run the last 2000 generations
	for(i = 0; i < global_params.num_exp_gen; i++)
	{
		t1 = std::chrono::high_resolution_clock::now();
		for (ii = 0; ii < global_params.num_pops; ii++)
		{
            
            pops[ii].determine_pop_size(global_params);    
			if (pops[ii].population_size > 0)
			{
				if (global_params.verbose)
				{
					if (global_params.log_file)
						log_out << '\n' << global_params.num_init_gen + i;
					else
						std::cout << '\n' << global_params.num_init_gen + i << std::flush;
				}
                else
                {
                    if (i % 1000 == 0)
                    {
                        if(global_params.log_file)
                            log_out << "\nExperimental generation " << i + 1 << " beginning.";
                        std::cout << "\nexperimental generation " << i + 1 << " beginning." << std::flush;
                    }
                }
										
				//mating (includes assiging preferences, recombination, and mutation)
				pops[ii].nest_and_fertilize(global_params, false, "temp");
				//selection
                if(global_params.viability_selection)
                    pops[ii].viability_selection(global_params);
				//output summary stats
				summary_output << "\n" << global_params.num_init_gen + i << "\tPop" << ii;
				pops[ii].output_summary_info(global_params, summary_output, log_out);
                markers_output << '\n' <<global_params.num_init_gen + i << "\tPop" << ii;
                pops[ii].output_allele_freqs(global_params, markers_output);
				ae_output << '\n' << global_params.num_init_gen +i+1 << "\tPop" << ii;
				pops[ii].output_allelic_effects(global_params,ae_output, false);
				if(i == (global_params.num_exp_gen - 1))
				{
					//output the trait values for the final generation of each population
					pops[ii].output_trait_info(global_params,global_params.num_init_gen + global_params.num_exp_gen, ii, trait_output, log_out);
				}
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
					{
						if (global_params.log_file)
							log_out << ", " << new_parent << " parents";
						else
							std::cout << ", " << new_parent << " parents" << std::flush;
					}
						
				}
				if (global_params.court_trait)
				{
					new_courter = pops[ii].calc_freq_courter(global_params);
					pops[ii].d_courterfreq.push_back((new_courter - courter_freqs[ii]));
					courter_freqs[ii] = new_courter;
					if (global_params.verbose)
					{
						if (global_params.log_file)
							log_out << ", " << new_courter << " courters";
						else
							std::cout << ", " << new_courter << " courters" << std::flush;
					}
				}
				
			}//if popsize > 0
			else
			{
				if (!pops[ii].extinct)
				{
					if (global_params.log_file)
                        log_out << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
					else
						std::cout << "\n" << global_params.base_name << ": Population" << ii << " has crashed at experimental generation " << i << '\n' << std::flush;
					pops[ii].output_trait_info(global_params, i, ii, trait_output, log_out);
					crash_counter++;
				}
					
				pops[ii].extinct = true;
				
			}
		}
		// stop if the populations have all crashed
		if (crash_counter == global_params.num_pops +1)
		{
			if (global_params.log_file)
				log_out<< "\n" << global_params.base_name << ": All populations have crashed at experimental generation " << i << '\n' << std::flush;
			else
				std::cout << "\n" << global_params.base_name << ": All populations have crashed at experimental generation " << i << '\n' << std::flush;
			trait_output.close();
			summary_output.close();
			markers_output.close();
			ae_output.close();
			if (command_line)
				return 0;
			else
			{
				std::cout << "\nInput integer to close dialogue.\n" << std::flush;
				cin >> iii;
				return 0;
			}
		}
		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
		if (global_params.optimize)
		{
			if (global_params.log_file)
				log_out  << "\n   experimental generation " << i << " took " << duration << " seconds.";
			else
				std::cout << "\n   experimental generation " << i << " took " << duration << " seconds.";
		}
	}
	
	//Evaluate stasis/equilibrium
	if (global_params.log_file)
		log_out << "\nEvaluating equilibrium" << std::flush;
	else
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
			if (global_params.log_file)
				log_out	<< "\nNo equilibrium could be reached for population " << i << " with population size " << pops[i].population_size << " by generation " << global_params.num_init_gen + global_params.num_exp_gen;
			else
				std::cout << "\nNo equilibrium could be reached for population " << i << " with population size " << pops[i].population_size << " by generation " << global_params.num_init_gen + global_params.num_exp_gen << std::flush;
            if(global_params.output_vcf)
                pops[i].output_genotypes_vcf(global_params, i);	
		}
	}
	
	//close output files
	summary_output.close();
	trait_output.close();
	markers_output.close();
	ae_output.close();
	if (global_params.log_file)
		log_out<< "\nDone!\n";
	else
		std::cout << "\nDone!\n" << std::flush;
    log_out.close();
	if (command_line)
		return 0;
	else
	{
		std::cout << "\nInput integer to quit.\n" << std::flush;
		cin >> i;
		return 0;
	}
}
