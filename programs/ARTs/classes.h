#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include "random_numbers.h"

using namespace std;


class chromosome
{
public:
	vector<int> loci;
	
	chromosome()
	{
		loci = vector<int>();
	}


};

class individual
{
public:
	double courter_trait, parent_trait, courter_gt, parent_gt; //use courter_trait as pref if it's a female
	vector<chromosome> maternal, paternal;
	bool female, alive, courter, parent;
	int mate_found;

	individual()
	{
		courter_trait = courter_gt = parent_trait = parent_gt = double();
		maternal = paternal = vector<chromosome>();
		mate_found = int();
		female = alive = parent = courter = bool();
	}



};

class parameters
{
public:
	int carrying_capacity, num_sampled, num_chrom, num_markers, num_qtl, max_fecund, max_encounters, num_alleles;
	int num_pops, num_init_gen, num_exp_gen, num_mal_traits, num_ld_comparisons;
	double mutation_rate, env_sd, rmu, recombination_rate, allelic_std_dev;
	string base_name;
	bool polygyny;

	parameters()
	{
		carrying_capacity = num_sampled = num_chrom = num_markers = num_qtl = max_fecund = max_encounters = num_alleles = int();
		num_pops = num_init_gen = num_exp_gen = int();
		mutation_rate = env_sd = rmu = recombination_rate = allelic_std_dev = double();
		polygyny = bool();
	}

	void set_defaults()
	{
		num_init_gen = 10000;//1000 for testing
		num_exp_gen = 2000;//200 for testing
		num_pops = 1;
		carrying_capacity = 1000;//1000
		num_sampled = 50;
		num_chrom = 4;
		num_markers = 1000;//1000
		num_qtl = 50 / num_chrom;//per chrom
		max_fecund = 4;
		max_encounters = 50;
		num_alleles = 2; //biallelic to start
		mutation_rate = 0.0002;
		num_mal_traits = 1;//either one or two
		env_sd = 0;
		rmu = 0;//this value can be altered
		recombination_rate = 0.2;//0.2
		allelic_std_dev = 0.5;
		polygyny = false;
		base_name = "../../results/arts_";
		num_ld_comparisons = 100;
	}

	void help_message()
	{
		cout << "\n\t\tHELP MENU\n";
		cout << "\nSimulation model of G-matrix stability and fsts. Below are the parameters to input to the model. (Defaults in parentheses)\n";
		cout << "-b:\tBase for output file names (gsim_)\n";
		cout << "-K:\tcarrying capacity (1000)\n";
		cout << "-s:\tnumber of individuals Sampled. (50)\n";
		cout << "-t:\tnumber of Traits (2)\n";
		cout << "-c:\tnumber of Chromosomes (4)\n";
		cout << "-x:\tnumber of markers per chromosome (1000)\n";
		cout << "-q:\ttotal number of Quantitative trait loci. (50)\n";
		cout << "-f:\tmaximum Fecundity. (4)\n";
		cout << "-e:\tmaximum number of Encounters during mating. (50)\n";
		cout << "-a:\tnumber of alleles (2).\n";
		cout << "-p:\tnumber of populations. (2)\n";
		cout << "-i:\tnumber of initial generations. (1000)\n";
		cout << "-g:\tnumber of experimental generations (200).\n";
		cout << "-v:\tenVironmental standard deviation (0).\n";
		cout << "-mu:\tmaximum mutation rate (0.0002).\n";
		cout << "-rmu:\trmu value (0).\n";
		cout << "-mv:\tMutational Variances, enter as many numbers as there are traits separated by commas (0.05,0.05)\n";
		cout << "-r:\tRecombination rate. (0.2) \n";
		cout << "-asd:\tAllelic Standard Deviation (0.5)\n";
		cout << "-m:\tMigration rate. (0.1)\n";
		cout << "--selection:\tSelection type. Choices are: none, Gaussian (Gaussian)\n";
		cout << "--selection-file:\tpath to file containing information on selection parameters.(selection.txt)\n\t\tFor more info use flat --selection-help.\n";
		cout << "--migration:\tMigration type. Choices are: all, ibd, cline, mainland-island, infinite-island (all)\n";
		cout << "--selection-help:\tOutputs details for selection file specification.\n";
		cout << "--polygyny:\tIf true, males will mate multiple times. If false, monogamy occurse (false)\n";
		cout << "-h:\tPrint this help message.\n";
	}

	void selection_help_message()
	{
		cout << "\n\t\tINSTRUCTIONS FOR SELECTION FILE\n";
		cout << "The first line should begin with 'SELECTION' followed by a space and one of the following keywords : gaussian, sexual_gaussian\n";
		cout << "The rest of the lines should be split into two sections: OMEGAS and THETAS.\n";
		cout << "Place the phrase 'OMEGA' on the line preceding the omega values and the phrase 'THETAS' on the line preceding the theta values.\n";
		cout << "Omega is a symmetric matrix with num_traits rows and columns. Each population should have its own matrix specified, ";
		cout << "so each population will have num_traits rows with num_traits tab-delimited columns.\n";
		cout << "The diagonal values are w11,w22, etc, but instead of w12 etc on the off-diagonal please provide the correlation coefficient (rw).\n";
		cout << "It is acceptable to have a blank line between poulation matrices.\n";
		cout << "Theta values should be represented in tab-delimited columns, one column for each trait. ";
		cout << "Each row following THETA represents a population.\n";
		cout << "The omega values describe the strength of selection on each trait, and smaller values represent stronger selection.\n";
		cout << "Theta values represent trait optima. ";
		cout << "If you would like to have the optima centered around the mean, input 'mean' instead of a number in the appropriate column/row.\n";
		cout << "If multiple episodes of selection are included, just repeat the same format after antoher SELECTION line.\n";
		cout << "Here's an example:\n\n";
		cout << "\tgaussian\n\tOMEGAS\n\t9\t0\n\t0\t9\n\t49\t0.5\n\t0.5\t49\n\tTHETAS\n\tmean\t0\n\tmean\t0";
		cout << "\n\nIf you are specifying sexual selection (sexual_gaussian), add the index for the trait experiencing sexual selection (e.g., 0 or 1) ";
		cout << "after the line with SELECTION sexual_gaussian and before OMEGAS or THETAS.\n";
	}

	bool parse_parameters(int argc, char*argv[])
	{
		int j;
		bool run_program = true;
		set_defaults();
		string tempstring1, tempstring2;
		if (argc == 1)
		{
			help_message();
			run_program = false;
		}
		if (argc > 1)
		{
			tempstring1 = argv[1];
			if (tempstring1 == "-h" || tempstring1 == "--selection-help")
			{
				if (tempstring1 == "-h")
					help_message();
				if (tempstring2 == "--selection-help")
					selection_help_message();
				run_program = false;
			}
			else
			{
				run_program = true;
				for (j = 1; j < argc - 1; j++)
				{
					tempstring1 = argv[j];
					tempstring2 = argv[j + 1];
					if (tempstring1 == "-b")
						base_name = tempstring2;
					if (tempstring1 == "-K")
						carrying_capacity = atoi(tempstring2.c_str());
					if (tempstring1 == "-s")
						num_sampled = atoi(tempstring2.c_str());
					if (tempstring1 == "-t")
						num_traits = atoi(tempstring2.c_str());
					if (tempstring1 == "-c")
						num_chrom = atoi(tempstring2.c_str());
					if (tempstring1 == "-x")
						num_markers = atoi(tempstring2.c_str());
					if (tempstring1 == "-q")
						num_qtl = atoi(tempstring2.c_str());
					if (tempstring1 == "-f")
						max_fecund = atoi(tempstring2.c_str());
					if (tempstring1 == "-e")
						max_encounters = atoi(tempstring2.c_str());
					if (tempstring1 == "-a")
						num_alleles = atoi(tempstring2.c_str());
					if (tempstring1 == "-p")
						num_pops = atoi(tempstring2.c_str());
					if (tempstring1 == "-i")
						num_init_gen = atoi(tempstring2.c_str());
					if (tempstring1 == "-g")
						num_exp_gen = atoi(tempstring2.c_str());
					if (tempstring1 == "-v")
						env_sd = atof(tempstring2.c_str());
					if (tempstring1 == "-mu")
						max_mutation_rate = atof(tempstring2.c_str());
					if (tempstring1 == "-rmu")
						rmu = atof(tempstring2.c_str());
					if (tempstring1 == "-r")
						recombination_rate = atof(tempstring2.c_str());
					if (tempstring1 == "-asd")
						allelic_std_dev = atof(tempstring2.c_str());
					if (tempstring1 == "-mv")
					{
						stringstream ss;
						ss << tempstring2;
						while (getline(ss, tempstring1, ','))
							mutational_variance.push_back(atof(tempstring1.c_str()));
					}
					if (tempstring1 == "--selection")
						selection_type = tempstring2;
					if (tempstring1 == "--selection-file")
						selection_file_name = tempstring2;
					if (tempstring1 == "--migration")
						migration_type = tempstring2;
					if (tempstring1 == "-m")
						migration_rate = atof(tempstring2.c_str());
					if (tempstring1 == "--polygyny")
					{
						if (tempstring2 == "T" | tempstring2 == "true" || tempstring2 == "TRUE" || tempstring2 == "True" || tempstring2 == "t")
							polygyny = true;
						else
							polygyny = false;
					}
				}
			}
		}
		return run_program;
	}

	void output_parameters()
	{
		ofstream param_out;
		string param_out_name = base_name + "parameters.txt";
		param_out.open(param_out_name);
		param_out << "Carrying capacity:\t" << carrying_capacity;
		param_out << "\nNum Sampled:\t" << num_sampled;
		param_out << "\nNum Chrom:\t" << num_chrom;
		param_out << "\nNum markers:\t" << num_markers;
		param_out << "\nNum QTL:\t" << num_qtl;
		param_out << "\nMaximum Fecundity:\t" << max_fecund;
		param_out << "\nMax encounters:\t" << max_encounters;
		param_out << "\nNum alleles:\t" << num_alleles;
		param_out << "\nNum Pops:\t" << num_pops;
		param_out << "\nNum Initial generations:\t" << num_init_gen;
		param_out << "\nNum Experimental generations:\t" << num_exp_gen;
		param_out << "\nNum Traits:\t" << num_traits;
		param_out << "\nNum LD comparisons:\t" << num_ld_comparisons;
		param_out << "\nMax mutation rate:\t" << max_mutation_rate;
		param_out << "\nEnvironmental SD:\t" << env_sd;
		param_out << "\nrmu:\t" << rmu;
		param_out << "\nRecombination rate:\t" << recombination_rate;
		param_out << "\nAllelic standard deviation:\t" << allelic_std_dev;
		param_out << "\nMigration_rate:\t" << migration_rate;
		param_out << "\nSelection type:\t" << selection_type;
		param_out << "\nMigration type:\t" << migration_type;
		param_out << "\nSelection file name:\t" << selection_file_name;
		param_out << "\nPolygyny:\t" << polygyny;
		for (int j = 0; j < mutational_variance.size(); j++)
			param_out << "\nTrait" << j << " mutational variance:\t" << mutational_variance[j];
		param_out.close();
	}
};
