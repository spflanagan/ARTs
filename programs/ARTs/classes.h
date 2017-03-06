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

class tracker
{
public:
	vector<double> per_locus;

	tracker()
	{
		per_locus = vector<double>();
	}


};

class chromosome
{
public:
	vector<int> loci;
	vector<double> courter_ae, parent_ae, cenv_ae, penv_ae;


	chromosome()
	{
		loci = vector<int>();
		courter_ae = parent_ae = cenv_ae = penv_ae = vector<double>();
	}


};

class individual
{
public:
	double courter_trait, parent_trait, courter_gt, parent_gt; //use courter_trait as pref if it's a female
	vector<chromosome> maternal, paternal;
	vector<tracker> courter_int, parent_int, pref_int;
	bool female, alive, courter, parent;
	int mate_found;

	individual()
	{
		courter_trait = courter_gt = parent_trait = parent_gt = double();
		maternal = paternal = vector<chromosome>();
		courter_int = parent_int = pref_int = vector<tracker>();
		mate_found = int();
		female = alive = parent = courter = bool();
	}

	void calc_courter_trait(double env_cue, parameters gp)
	{
		int k, kk, kkk;
		courter_trait = 0;
		if (!gp.env_effects)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.num_alleles; kk++)
					courter_trait = courter_trait + maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk];
			}
		}
		else
		{
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.num_env_qtl; kk++)
				{
					for (kkk = 0; kkk < (gp.num_env_qtl + gp.num_qtl); kkk++)
					{
						courter_trait = courter_trait +
							((maternal[k].cenv_ae[kk] + paternal[k].cenv_ae[kk]) * courter_int[qtl_index].per_locus[kkk])
							+ env_cue;
					}
					qtl_index++;
				}
				for (kk = 0; kk < gp.num_qtl; kk++)
				{
					for (kkk = 0; kkk < (gp.num_env_qtl + gp.num_qtl); kkk++)
					{
						courter_trait = courter_trait +
							((maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk]) * courter_int[qtl_index].per_locus[kkk]);
					}
					qtl_index++;
				}
			}
			courter_trait = 1 / (1 + exp(-1 * courter_trait));
		}
	}

	void calc_parent_trait(double env_cue, parameters gp)
	{
		int k, kk, kkk;
		parent_trait = 0;
		if (!gp.env_effects)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.num_alleles; kk++)
					parent_trait = parent_trait + maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk];
			}
		}
		else
		{
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.num_env_qtl; kk++)
				{
					for (kkk = 0; kkk < (gp.num_env_qtl + gp.num_qtl); kkk++)
					{
						parent_trait = parent_trait +
							((maternal[k].penv_ae[kk] + paternal[k].penv_ae[kk]) * parent_int[qtl_index].per_locus[kkk])
							+ env_cue;
					}
					qtl_index++;
				}
				for (kk = 0; kk < gp.num_qtl; kk++)
				{
					for (kkk = 0; kkk < (gp.num_env_qtl + gp.num_qtl); kkk++)
					{
						parent_trait = parent_trait +
							((maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk]) * parent_int[qtl_index].per_locus[kkk]);
					}
					qtl_index++;
				}
			}
			parent_trait = 1 / (1 + exp(-1 * parent_trait));
		}
	}
};

class parameters
{
public:
	int carrying_capacity, num_sampled, num_chrom, num_markers, num_qtl, num_env_qtl, max_fecund, max_encounters, num_alleles;
	int num_pops, num_init_gen, num_exp_gen, num_ld_comparisons;
	double mutation_rate, recombination_rate, allelic_std_dev;
	string base_name;
	bool court_trait, parent_trait, env_effects, cor_prefs, ind_pref, FD_pref, CD_pref, FD_court, FD_parent,CD_court, CD_parent;

	parameters()
	{
		carrying_capacity = num_sampled = num_chrom = num_markers = num_qtl = max_fecund = max_encounters = num_alleles = int();
		num_pops = num_init_gen = num_exp_gen = int();
		mutation_rate =recombination_rate = allelic_std_dev = double();
		env_effects = court_trait = parent_trait = cor_prefs = ind_pref = FD_pref = CD_pref = FD_court = FD_parent = CD_court = CD_parent = bool();
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
		num_env_qtl = 0;
		max_fecund = 4;
		max_encounters = 50;
		num_alleles = 2; //biallelic to start
		mutation_rate = 0.0002;
		recombination_rate = 0.2;//0.2
		allelic_std_dev = 0.5;
		court_trait = true;
		parent_trait = false;
		env_effects = false;
		cor_prefs= ind_pref= FD_pref= CD_pref= FD_court= FD_parent= CD_court= CD_parent = false;//no selection
		base_name = "../../results/arts_";
		num_ld_comparisons = 100;
	}

	void help_message()
	{
		cout << "\n\t\tHELP MENU\n";
		cout << "\nSimulation model of G-matrix stability and fsts. Below are the parameters to input to the model. (Defaults in parentheses)\n";
		cout << "-b:\tBase for output file names (arts_)\n";
		cout << "-K:\tcarrying capacity (1000)\n";
		cout << "-s:\tnumber of individuals Sampled. (50)\n";
		cout << "-c:\tnumber of Chromosomes (4)\n";
		cout << "-x:\tnumber of markers per chromosome (1000)\n";
		cout << "-q:\ttotal number of Quantitative trait loci. (50)\n";
		cout << "-eq:\ttotal number of QTL responding the the Environment (half of -q if --plasticity flag included)\n";
		cout << "-f:\tmaximum Fecundity. (4)\n";
		cout << "-e:\tmaximum number of Encounters during mating. (50)\n";
		cout << "-a:\tnumber of alleles (2).\n";
		cout << "-p:\tnumber of populations. (2)\n";
		cout << "-i:\tnumber of initial generations. (1000)\n";
		cout << "-g:\tnumber of experimental generations (200).\n";
		cout << "-mu:\tmaximum mutation rate (0.0002).\n";
		cout << "-r:\tRecombination rate. (0.2) \n";
		cout << "-asd:\tAllelic Standard Deviation (0.5)\n";
		cout << "--plasticity:\tModel a plastic morph, where genotype is determined by interactions between genes.\n";
		cout << "--freq-dependent-preference:\tInclude this flag if preferences should be based on the frequency of male morphs (defaults to independent).\n";
		cout << "--condition-dependent-preference:\tInclude this flag if female preferences are determined by female condition (defaults to independent).\n";
		cout << "--courter:\tInclude this flag if males should have the courter trait (a trait affecting mating probabilities). Defaults to Gaussian preference for randomly chosen morph unless other flags included. \n";
		cout << "--parent:\tInclude this flag if males should have the parental trait (a trait affecting offspring survival). Defaults to Gaussian preference for randomly chosen morph unless other flags included. \n";
		cout << "--freq-dependent-courter:\tIf the courter trait is experiencing frequency dependent selection.\n";
		cout << "--freq-dependent-parent:\tIf the parent trait is experiencing frequency dependent selection.\n";
		cout << "--condition-dependent-courter:\tIf the courter trait is determined by male condition.\n";
		cout << "--condition-dependent-parent:\tIf the parent trait is determined by male condition.\n";
		cout << "--independent-pref:\tSpecifies an independent female preference (defaults to Gaussian preference for randomly chosen morph unless other flags included). \n";
		cout << "--correlated-pref:\tSpecifies a female preference correlated with the male courter trait (defaults to Gaussian preference for randomly chosen morph unless other flags included).\n";
		cout << "-h:\tPrint this help message.\n";
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
				run_program = false;
			}
			else
			{
				run_program = true;
				for (j = 1; j < argc - 1; j++)
				{
					tempstring1 = argv[j];
					if(tempstring1.substr(0,2) != "--")
						tempstring2 = argv[j + 1];
					if (tempstring1 == "-b")
						base_name = tempstring2;
					if (tempstring1 == "-K")
						carrying_capacity = atoi(tempstring2.c_str());
					if (tempstring1 == "-s")
						num_sampled = atoi(tempstring2.c_str());
					if (tempstring1 == "-c")
						num_chrom = atoi(tempstring2.c_str());
					if (tempstring1 == "-x")
						num_markers = atoi(tempstring2.c_str());
					if (tempstring1 == "-q")
						num_qtl = atoi(tempstring2.c_str());
					if (tempstring1 == "-eq")
						num_env_qtl = atoi(tempstring2.c_str());
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
					if (tempstring1 == "-mu")
						mutation_rate = atof(tempstring2.c_str());
					if (tempstring1 == "-r")
						recombination_rate = atof(tempstring2.c_str());
					if (tempstring1 == "-asd")
						allelic_std_dev = atof(tempstring2.c_str());
					if (tempstring1 == "--plasticity")
						env_effects = true;
					if (tempstring1 == "--freq-dependent-preference")
						FD_pref = ind_pref = true;
					if (tempstring1 == "--condition-dependent-preference")
						CD_pref = ind_pref = true;
					if (tempstring1 == "--courter")
						court_trait = ind_pref = true;
					if (tempstring1 == "--parent")
						parent_trait = ind_pref = true;
					if (tempstring1 == "--freq-dependent-courter")
						FD_court = court_trait= true;
					if (tempstring1 == "--freq-dependent-parent")
						FD_parent = parent_trait = true;
					if (tempstring1 == "--condition-dependent-courter")
						CD_court = court_trait= true;
					if (tempstring1 == "--condition-dependent-parent")
						CD_parent =parent_trait= true;
					if (tempstring1 == "--independent-pref")
						ind_pref = true;
					if (tempstring1 == "--correlated-pref")
						cor_prefs = court_trait = true;
				}
				if (env_effects)
				{
					if (num_env_qtl == 0)
						num_env_qtl = num_qtl / 2;
					else
						num_env_qtl = num_env_qtl / num_chrom;
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
		param_out << "\nNum LD comparisons:\t" << num_ld_comparisons;
		param_out << "\nMutation rate:\t" << mutation_rate;
		param_out << "\nRecombination rate:\t" << recombination_rate;
		param_out << "\nAllelic standard deviation:\t" << allelic_std_dev;
		if (env_effects)
			param_out << "\nPlasticity";
		if (FD_pref)
			param_out << "\n--freq-dependent-preference";	
		if (CD_pref)
			param_out << "\n--condition-dependent-preference";
		if (court_trait)
			param_out << "\n--courter";
		if (parent_trait)
			param_out << "\n--parent";
		if (FD_court)
			param_out << "\n--freq-dependent-courter";	 
		if (FD_parent)
			param_out << "\n--freq-dependent-parent";
		if (CD_court)
			param_out << "\n--condition-dependent-courter";
		if (CD_parent)
			param_out << "\n--condition-dependent-parent";
		if (ind_pref)
			param_out << "\n--independent-pref";
		if (cor_prefs )
			param_out<< "\n--correlated-pref";
			
		param_out.close();
	}
};
