#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include "random_numbers.h"
#include <chrono>

using namespace std;

//This is adapted from random_shuffle
//and the Knuth-Fisher-Yates shuffle algorithm
//adapted from: https://blog.codinghorror.com/the-danger-of-naivete/
void shuffle_vec(vector<int> & vec)
{
	int temp;
	for (int i = vec.size() - 1; i > 0; i--)
	{
		int n = randnum(i + 1);
		temp = vec[i];
		vec[i] = vec[n];
		vec[n] = temp;
	}
}

class tracker
{
public:
	vector<double> per_locus;

	tracker()
	{
		per_locus = vector<double>();
	}

	int first_qtl()
	{
		int q;
		double min = per_locus[0];
		for (q = 0; q < per_locus.size(); q++)
		{
			if (per_locus[q] < min)
				min = per_locus[q];
		}
		return min;
	}

	int last_qtl()
	{
		int q;
		double max = per_locus[0];
		for (q = 0; q < per_locus.size(); q++)
		{
			if (per_locus[q] > max)
				max = per_locus[q];
		}
		return max;
	}
};

class chromosome
{
public:
	vector<int> loci;
	vector<double> courter_ae, parent_ae, pref_ae;
	vector<double> courter_thresh, parent_thresh;


	chromosome()
	{
		loci = vector<int>();
		courter_ae = parent_ae = pref_ae = courter_thresh = parent_thresh = vector<double>();
	}


};

class nest
{
public:
	int mom, nest_dad;
	vector<int> all_dads;
	vector<double> off_props;

	nest()
	{
		mom = nest_dad = int();
		all_dads = vector<int>();
		off_props = vector<double>();
	}
};

class parameters
{
public:
	int carrying_capacity, num_sampled, num_chrom, num_markers, num_qtl, num_env_qtl, max_fecund, max_encounters, num_alleles;
	int num_pops, num_init_gen, num_exp_gen, num_ld_comparisons, rs_c, rs_nc, rs_p, rs_np, max_num_mates;
	double mutation_rate, mutational_var, recombination_rate, allelic_std_dev, gaussian_pref_mean, cond_adj, via_sel_strength, supergene_prop;
	double sperm_comp_r, egg_surv_parent, egg_surv_noparent;
	string base_name;
	bool same_base, court_trait, parent_trait, gene_network, env_cue, cor_prefs, ind_pref, FD_pref, CD_pref, FD_court, FD_parent,CD_court, CD_parent, polygyny, cor_mal_traits;
    bool density_dependent,all_sneak, per_fem_mating, supergene, random_mating, courter_conditional, parent_conditional, thresholds_evolve, thresholds_in_supergene, verbose, no_genetics, linked_additive,optimize, output_vcf, viability_selection;
	vector <int> qtl_per_chrom;
	bool log_file, ae_vcf;

	parameters()
	{
		carrying_capacity = num_sampled = num_chrom = num_markers = num_qtl = max_fecund = max_encounters = num_alleles = int();
		num_pops = num_init_gen = num_exp_gen = max_num_mates = int();
		mutation_rate =recombination_rate = allelic_std_dev = egg_surv_parent = egg_surv_noparent = sperm_comp_r = double();
		same_base = gene_network = env_cue = court_trait = parent_trait = cor_prefs = ind_pref = FD_pref = CD_pref = FD_court = FD_parent = CD_court = CD_parent = polygyny = cor_mal_traits = supergene =  bool();
		density_dependent = all_sneak = per_fem_mating = random_mating = courter_conditional = parent_conditional = thresholds_evolve = thresholds_in_supergene = no_genetics=linked_additive = optimize= output_vcf=viability_selection = bool();
		qtl_per_chrom = vector<int>();
		log_file = ae_vcf = bool();
	}

	void set_defaults()
	{
		same_base = true; //does it always start with the same base populations for all of the populations?
		num_pops = 1; //default 1 
		num_init_gen = 10000;//10000 normally
		num_exp_gen = 2000;//2000 normally		
		carrying_capacity = 1000;//1000
		num_sampled = 50;//default 50
		num_chrom = 4;//default 4
		num_markers = 1000;//1000
		num_qtl = 50;//default: 50
		num_env_qtl = 0;//default 0
		max_fecund = 4;//default 4
		max_encounters = 50;//default 50
		num_alleles = 2; //biallelic to start
		mutation_rate = 0.0002;//default 0.0002
		mutational_var = 0;//default 0
		recombination_rate = 0.2;//0.2
		allelic_std_dev = 0.5;//default 0.5
		supergene_prop = 0.1; //default 0.1 (10% of num_markers = 100)
		density_dependent = true; //default: true (sets density dependent mate choice)
		random_mating = true;//default: true
		all_sneak = false; //default: false
		per_fem_mating = false; //the original implementation of mating
		supergene =  false; //default: false
		court_trait = courter_conditional = false; //default: false
		parent_trait = parent_conditional = false;//default false
		linked_additive = true; //default: true, traits are determined by additive genetic variance.
		thresholds_evolve = thresholds_in_supergene = false; //default: false
		env_cue = gene_network = false;//default false
		cor_prefs= ind_pref= FD_pref= CD_pref= FD_court= FD_parent= CD_court= CD_parent = false;//no selection
		cor_mal_traits = false;//default false
		polygyny = false;//default false
		base_name = "../../results/arts";//default: "../../results/arts"
		num_ld_comparisons = 100;//default 100
		rs_c = 8;//default 8
		rs_nc = 4;//default 4
		rs_p = 8;//default 8
		rs_np = 4;//default 4
		max_num_mates = 3;//default 3 - used to cap number of sneakers allowed to fertilize a nest
		egg_surv_noparent = 0.1;//default 0.1
		egg_surv_parent = 0.9;//default 0.9
		sperm_comp_r = 0.5; //default 0.5
		gaussian_pref_mean = 0;//default 0
		via_sel_strength = 50;//unsure what value to put here (currently 50)
		cond_adj = 0.1;//amount to add/subtract to condition dependent traits (default 0.1)
		verbose = false; //outputs info about every step during every initial generation -- good for debugging
		no_genetics = false; //removes genetic architecture, is just unlinked additive genetic variance
		optimize = false; //if true, it outputs time taken for each step. (default false)
        output_vcf = false; //if true, a vcf will be output
        viability_selection = false; // if true, viability selection will act on offspring
		log_file = true; //saves std out to log instead of outputting it to the console
		ae_vcf = false; //creates a vcf of allelic effects
	}

	void help_message()
	{
		std::cout << "\n\t\tHELP MENU\n";
		std::cout << "\nSimulation model of the evolution of alternative reproductive tactics with complex genetic architectures.\n";
		std::cout << "Below are the parameters to input to the model. (Defaults in parentheses)\n";
		std::cout << "-b:\tBase for output file names (arts)\n";
		std::cout << "-K:\tcarrying capacity (1000)\n";
		std::cout << "-s:\tnumber of individuals Sampled. (50)\n";
		std::cout << "-c:\tnumber of Chromosomes (4)\n";
		std::cout << "-x:\tnumber of markers per chromosome (1000)\n";
		std::cout << "-q:\ttotal number of Quantitative trait loci. (50)\n";
		std::cout << "-eq:\ttotal number of QTL responding the the Environment (half of -q if --plasticity flag included)\n";
		std::cout << "-f:\tmaximum Fecundity. (4)\n";
		std::cout << "-e:\tmaximum number of Encounters during mating. (50)\n";
		std::cout << "-a:\tnumber of alleles (2).\n";
		std::cout << "-p:\tnumber of populations. (1)\n";
		std::cout << "-i:\tnumber of initial generations. (10000)\n";
		std::cout << "-g:\tnumber of experimental generations (2000).\n";
		std::cout << "-mu:\tmaximum mutation rate (0.0002).\n";
		std::cout << "-v:\tviability selection strength against courters and parents (50).\n";
		std::cout << "-r:\tRecombination rate. (0.2) \n";
		std::cout << "-asd:\tAllelic Standard Deviation (0.5)\n";
		std::cout << "-prs:\tParental male reproductive success (8)\n";
		std::cout << "-nprs:\tNon-parental male reproductive success (4)\n";
		std::cout << "-crs:\tCourter male reproductive success (8)\n";
		std::cout << "-ncrs:\tNon-courter male reproductive success (4)\n";
		std::cout << "-sprop:\tSupergene proportion of a chromosome (0.1). NOTE: must be > the total number of qtls.\n";
		std::cout << "-sperm-r:\tSperm competition r. If r = 1, paternity is determined through a fair raffle. \n\tIf 0 < r < 1, the parental male has a higher share of paternity and the sneakers fertilize with penalty r (rs2/(s1+s2); default 0.5)\n";
		std::cout << "-surv-noparent:\tSurvival probability of eggs without a parent (0.1)\n";
		std::cout << "-surv-parent:\tSurvival probability with a parent (0.9)\n";
		std::cout << "-mm:\tMax number of Mates. This includes the chosen male (3)\n";
		std::cout << "-qpc:\tQTLs per chromosome (-q/-c, aka an even distribution of QTLs among chromosomes). To specify different numbers per chromosome separate them by commas.\n\tExample: for 4 chromosomes, input: 4,2,3,0\n";
		std::cout << "--gene-network:\tModel a plastic morph, where genotype is determined by interactions between genes.\n";
		std::cout << "--env-cue:\tModel a gene network (--gene-network) that incorporates social information as an environmental cue.\n";
		std::cout << "--freq-dependent-preference:\tInclude this flag if preferences should be based on the frequency of male morphs (defaults to independent).\n";
		std::cout << "--condition-dependent-preference:\tInclude this flag if female preferences are determined by female condition (defaults to independent).\n";
		std::cout << "--courter:\tInclude this flag if males should have the courter trait (a trait affecting mating probabilities). Defaults to preference for randomly chosen morph unless other flags included. \n";
		std::cout << "--parent:\tInclude this flag if males should have the parental trait (a trait affecting offspring survival). Defaults to preference for randomly chosen morph unless other flags included. \n";
		std::cout << "--freq-dependent-courter:\tIf the courter trait is experiencing frequency dependent selection.\n";
		std::cout << "--freq-dependent-parent:\tIf the parent trait is experiencing frequency dependent selection.\n";
		std::cout << "--condition-dependent-courter:\tIf the courter trait is influenced by male condition.\n";
		std::cout << "--condition-dependent-parent:\tIf the parent trait is influenced by male condition.\n";
		std::cout << "--independent-pref:\tSpecifies an independent female preference (defaults to Gaussian preference for randomly chosen morph unless other flags included). \n";
		std::cout << "--correlated-pref:\tSpecifies a female preference correlated with the male courter trait (defaults to Gaussian preference for randomly chosen morph unless other flags included).\n";
		std::cout << "--random-mating:\tSpecifies no female choice (default: true).\n";
		std::cout << "--all-sneak:\tSpecifies all males sneak, not just sneakers (default: false).\n";
		std::cout << "--supergene:\tSpecifies whether the QTLs are grouped together in a supergene that has reduced recombination.\n";
		std::cout << "--polygyny:\tAllows males to mate multiply (default: true).\n";
		std::cout << "--courter-conditional:\tIf the courter trait has no genetic basis and is determined randomly or through environmental effects.\n";
		std::cout << "--parent-conditional:\tIf the parent trait has no genetic basis and is determined randomly or through environmental effects.\n";
		std::cout << "--thresholds-evolve:\tIf the thresholds are allowed to evolve (i.e., they have a quantitative genetic basis).\n";
		std::cout << "--thresholds-in-supergene:\tThe thresholds have a genetic basis and their loci are in the supergene.\n";
		std::cout << "--verbose:\tOutputs info about every step during every initial generation -- good for debugging";
		std::cout << "--no-genetics:\tRemoves genetic architecture; traits encoded by heritable unlinked additive genetic variance.\n";
		std::cout << "--linked-additive:\t(default) Traits are determined by genome-wide additive genetic variance distributed among chromosomes.\n";
		std::cout << "--viability:\tIf included, viability selection acts on offspring. If not, viability selection is turned off.\n";
		std::cout << "--density-independent:\tTurns off density-dependent selection. This is not recommended as it will likely lead to population crashes.\n";
		std::cout << "--optimize:\tOutput time steps for initial generations, for optimizing the code. (default is false)\n";
		std::cout << "--same-base:\tStart each replicate population with the same base population (default is true)\n";
        std::cout << "--output-vcf:\tInclude vcf output for all genotypes of individuals (default is false)\n";
		std::cout << "--ae-vcf:\tInclude vcf output for the allelic effects for QTLs for all individuals in generation 0 (default is false)\n";
		std::cout << "--log-file:\tSave output to logfile instead of std::cout\n";
		std::cout << "-h or --help:\tPrint this help message.\n";
	}

	void dependent_params(string qpc ="")
	{
		int j, qtl_counter;
		//assign values based on input
		if (no_genetics)
		{
			gene_network = supergene = false;
			num_chrom = 1;
			num_markers = 0;
			recombination_rate = 0;
		}
		if (qpc.size() > 0)
		{
			istringstream ss(qpc);
			qtl_counter = 0;
			for (j = 0; j < num_chrom; j++)
			{
				string temp;
				getline(ss, temp, ',');
				qtl_per_chrom[j] = atoi(temp.c_str());
				qtl_counter = qtl_counter+ qtl_per_chrom[j];
			}
		}
		else
		{
			qtl_counter = 0;
			for (j = 0; j < num_chrom; j++)
			{
				qtl_per_chrom.push_back(num_qtl / num_chrom);//default is an even distribution
				qtl_counter = qtl_counter +qtl_per_chrom[j];
			}
		}
		num_qtl = qtl_counter;
		//set other parameters to adhere to model requirements/assumptions
		if (num_init_gen < 2)
			num_init_gen = 2;
		if (gene_network)
		{
			if (num_env_qtl == 0)
				num_env_qtl = num_qtl / 2;
			else
				num_env_qtl = num_env_qtl / num_chrom;
			if (!court_trait && !parent_trait && !ind_pref && !cor_prefs)//there needs to be a genetically based trait!
				court_trait = true;
		}
		if (ind_pref || cor_prefs || FD_pref || CD_pref)
		{//if there are female preferences but neither male trait is specified,
			if (!court_trait && !parent_trait && !parent_conditional) //then the courter trait exists without genetic basis
				courter_conditional = true;
			random_mating = false;
		}
		if (CD_pref)
		{
			if (!ind_pref && !cor_prefs)
				ind_pref = true;
			if (!court_trait && !parent_trait && !parent_conditional)
				courter_conditional = true;
			random_mating = false;
		}
		if (env_cue)
			gene_network = true;
		if (!gene_network)
			env_cue = false;
		if (thresholds_in_supergene)
			thresholds_evolve = supergene = true;
		if (thresholds_evolve)
		{//then there must be a trait
			if (!court_trait && !parent_trait && !parent_conditional) //then the courter trait exists without genetic basis
				courter_conditional = true;
		}
		//these are non-overlapping traits
		if (court_trait)
			courter_conditional = false;
		if (courter_conditional)
			court_trait = false;
		if (parent_trait)
			parent_conditional = false;
		if (parent_conditional)
			parent_trait = false;
		if (cor_prefs)
			ind_pref = false;
		if (ind_pref)
			cor_prefs = false;
		if (court_trait || courter_conditional)// if there's a courtship trait then there must be a female preference
		{
			if (!cor_prefs && !ind_pref)
				random_mating = false;
		}
		//if (parent_conditional || parent_trait)
		//{
		//	if (!courter_conditional && !court_trait)// if there's no courtship trait and no female preferences specified
		//	{										//then females choose a nest randomly. 
		//		if (!ind_pref && !cor_prefs)
		//			random_mating = true;
		//	}
		//}
		if (!no_genetics & !gene_network & !supergene & !courter_conditional & !parent_conditional)
			linked_additive = true;
		else
			linked_additive = false;
		if(per_fem_mating && density_dependent)
			per_fem_mating = false;
	}

	bool parse_parameters(int argc, char*argv[])
	{
		int j;
		string qpc;
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
			if (tempstring1 == "-h" || tempstring1 == "--help")
			{
				help_message();
				run_program = false;
			}
			else
			{
				run_program = true;
				for (j = 1; j < argc; j++)
				{
					tempstring1 = argv[j];
					if (tempstring1.substr(0, 2) != "--")
					{
						tempstring2 = argv[j + 1];
						j++;
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
						if (tempstring1 == "-prs")
							rs_p = atoi(tempstring2.c_str());
						if (tempstring1 == "-nprs")
							rs_np = atoi(tempstring2.c_str());
						if (tempstring1 == "-crs")
							rs_c = atoi(tempstring2.c_str());
						if (tempstring1 == "-ncrs")
							rs_nc = atoi(tempstring2.c_str());
						if (tempstring1 == "-sprop")
							supergene_prop = atof(tempstring2.c_str());
						if (tempstring1 == "-qpc")
							qpc = tempstring2;
						if (tempstring1 == "-surv-noparent")
							egg_surv_noparent = atof(tempstring2.c_str());
						if (tempstring1 == "-surv-parent")
							egg_surv_parent = atof(tempstring2.c_str());
						if (tempstring1 == "-sperm-r")
							sperm_comp_r = atof(tempstring2.c_str());
						if (tempstring1 == "-mm")
							max_num_mates = atoi(tempstring2.c_str());
						if (tempstring1 == "-v")
							via_sel_strength = atof(tempstring2.c_str());
					}
					else
					{
						if (tempstring1 == "--gene-network")
							gene_network = true;
						if (tempstring1 == "--env-cue")
							env_cue = true;
						if (tempstring1 == "--freq-dependent-preference")
							FD_pref = true;
						if (tempstring1 == "--condition-dependent-preference")
							CD_pref = true;
						if (tempstring1 == "--courter")
							court_trait = true;
						if (tempstring1 == "--parent")
							parent_trait = true;
						if (tempstring1 == "--freq-dependent-courter")
							FD_court = court_trait = true;
						if (tempstring1 == "--freq-dependent-parent")
							FD_parent = parent_trait = true;
						if (tempstring1 == "--condition-dependent-courter")
							CD_court = court_trait = polygyny = true;
						if (tempstring1 == "--condition-dependent-parent")
							CD_parent = parent_trait = polygyny = true;
						if (tempstring1 == "--independent-pref")
							ind_pref = true;
						if (tempstring1 == "--correlated-pref")
							cor_prefs = court_trait = true;
						if (tempstring1 == "--all-sneak")
							all_sneak = true;
						if (tempstring1 == "--random-mating")
						{
							cor_prefs = ind_pref = false;
							random_mating = true;
						}
						if (tempstring1 == "--supergene")
							supergene = true;
						if (tempstring1 == "--polygyny")
							polygyny = true;
						if (tempstring1 == "--courter-conditional")
							courter_conditional = true;
						if (tempstring1 == "--parent-conditional")
							parent_conditional = true;
						if (tempstring1 == "--thresholds-evolve")
							thresholds_evolve = true;
						if (tempstring1 == "--thresholds-in-supergene")
							thresholds_in_supergene = true;
						if (tempstring1 == "--verbose")
							verbose = true;
						if (tempstring1 == "--no-genetics")
							no_genetics = true;
						if (tempstring1 == "--linked-additive")
							linked_additive = true;
						if (tempstring1 == "--optimize")
							optimize = true;
						if (tempstring1 == "--same-base")
							same_base = true;
                        if(tempstring1 == "--output-vcf")
                            output_vcf = true;
						if(tempstring1 == "--ae-vcf")
                            ae_vcf = true;
						if (tempstring1 == "--log-file")
							log_file = true;
						if(tempstring1 == "--density-independent")
							density_dependent = false;
						if(tempstring1 == "--viability")
							viability_selection = true;
					}
				}
				
			}
		}

		dependent_params(qpc);

		return run_program;
	}

	void output_parameters()
	{
		ofstream param_out;
		string param_out_name = base_name + "_parameters.txt";
		param_out.open(param_out_name);
		param_out << "Carrying capacity:\t" << carrying_capacity;
		param_out << "\nNum Sampled:\t" << num_sampled;
		param_out << "\nNum Chrom:\t" << num_chrom;
		param_out << "\nNum markers:\t" << num_markers;
		param_out << "\nTotal Num QTL:\t" << num_qtl;
		for (int j = 0; j < num_chrom; j++)
			param_out << "\n" << qtl_per_chrom[j] << " QTL on Chrom " << j;
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
		param_out << "\nGaussian preference mean:\t" << gaussian_pref_mean;
		param_out << "\nViability selection strength:\t" << via_sel_strength;
		param_out << "\nSperm competition r:\t" << sperm_comp_r;
		param_out << "\nEgg survival without parent:\t" << egg_surv_noparent;
		param_out << "\nEgg survival with parent:\t" << egg_surv_parent;
		param_out << "\nMaximum number of mates:\t" << max_num_mates;
		if (gene_network)
			param_out << "\nGene Network";
		if (env_cue)
			param_out << "\nGene network with environmental cues";
		if (FD_pref)
			param_out << "\n--freq-dependent-preference";	
		if (CD_pref)
			param_out << "\n--condition-dependent-preference";
		if (court_trait) {
			param_out << "\n--courter";
			param_out << "\nReproductive success for courter:\t" << rs_c
				<< "\nReproductive success for non-courter:\t" << rs_nc;
		}
		if (parent_trait) {
			param_out << "\n--parent";
			param_out << "\nReproductive success for parent:\t" << rs_p
				<< "\nReproductive success for non-parent:\t" << rs_np;
		}
		if (FD_court)
			param_out << "\n--freq-dependent-courter";	 
		if (FD_parent)
			param_out << "\n--freq-dependent-parent";
		if (CD_court)
			param_out << "\n--condition-dependent-courter";
		if (CD_parent)
			param_out << "\n--condition-dependent-parent";
		if (CD_court || CD_parent)
			param_out << "\nCondition dependent adjustment amount:\t" << cond_adj;
		if (ind_pref)
			param_out << "\n--independent-pref";
		if (cor_prefs)
			param_out<< "\n--correlated-pref";
		if (polygyny)
			param_out << "\n--polygyny";
		if (supergene)
			param_out << "\n--supergene\n--" << supergene_prop;
		if (parent_conditional)
			param_out << "\n--parent_conditional";
		if (courter_conditional)
			param_out << "\n--courter_conditional";
		if (thresholds_evolve)
			param_out << "\n--thresholds_evolve";
		if (thresholds_in_supergene)
			param_out << "\n--thresholds_in_supergene";
		if (verbose)
			param_out << "\n--verbose";
		if (no_genetics)
			param_out << "\n--no-genetics";
		if (linked_additive)
			param_out << "\n--linked-additive";
		if (all_sneak)
			param_out << "\n--all-sneak";
        if(output_vcf)
            param_out << "\n--output-vcf";
		if(ae_vcf)
            param_out << "\n--ae-vcf";
		if(density_dependent)
			param_out << "\ndensity dependent";
		else
			param_out << "\ndensity independent";
		if(per_fem_mating)
			param_out << "\nper-female mating (deprecated)";
		if(viability_selection)
			param_out <<"\nviability selection";
        if(log_file)
            param_out << "\n--log-file ";
		param_out.close();
	}
};

string determine_date()
{
	time_t now = std::time(0);
	struct tm ltm;
	const size_t size = 1024;
	char buffer[size];

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__))
	localtime_s(&ltm, &now);
	
#else    // POSIX
	localtime_r(&now, &ltm);
	
#endif
	
	std::strftime(buffer, size, "%Y%m%d", &ltm);
	return buffer;
}

class individual
{
public:
	double courter_trait, parent_trait, female_pref,pref_trait, ind_cthresh, ind_pthresh;
	vector<chromosome> maternal, paternal;
	vector<tracker> courter_int, parent_int, pref_int, cthresh_int, pthresh_int, courter_Y, parent_Y, pref_Y,cthresh_Y,pthresh_Y;
	vector<int> courter_Z, parent_Z, pref_Z, cthresh_Z, pthresh_Z;
	vector<double> courter_x, parent_x, pref_x, cthresh_x, pthresh_x;
	bool female, alive, courter, parent;
	int mate_found, pot_rs, mate_id,lifetime_rs, mom,dad;

	individual()
	{
		courter_trait = parent_trait = female_pref = pref_trait = ind_cthresh = ind_pthresh = double();
		maternal = paternal = vector<chromosome>();
		courter_int = parent_int = pref_int =cthresh_int = pthresh_int, courter_Y = parent_Y = pref_Y = cthresh_Y=pthresh_Y =vector<tracker>();
		courter_Z = parent_Z = pref_Z = cthresh_Z = pthresh_Z = vector<int>();
		courter_x = parent_x = pref_x = cthresh_x = pthresh_x = vector<double>();
		mate_found = pot_rs = mate_id = lifetime_rs =mom = dad = int();
		female = alive = parent = courter = bool();
	}

	//reset
	void reset(bool living)
	{
		mate_id = dad = mom = int();
		pot_rs = mate_found = lifetime_rs = 0;
		alive = living;
	}
	//phenotype functions
	void assign_threshold_gt(parameters gp, int id, vector<double>&tempalleles, bool assign_courter)
	{
		int k, kk;
		for (k = 0; k < gp.num_chrom; k++)
		{
			for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
			{
				if (assign_courter)
				{
					maternal[k].courter_thresh[kk] = tempalleles[id%gp.num_alleles];
					paternal[k].courter_thresh[kk] = tempalleles[id%gp.num_alleles];
				}
				else
				{
					maternal[k].parent_thresh[kk] = tempalleles[id%gp.num_alleles];
					paternal[k].parent_thresh[kk] = tempalleles[id%gp.num_alleles];

				}
			}
		}
	}
	void gene_interactions(parameters gp, vector<tracker> & Y, vector<tracker> & in, vector<double>&X,double env_cue, bool stability = true)
	{
		int k, kk, t;
		double diff, x, xlast;
		alive = true;
		for (k = 0; k < gp.num_qtl; k++)
		{
			xlast = 1;
			x = 0;
			for (t = 0; t < 20; t++)
			{

				for (kk = 0; kk < gp.num_qtl; kk++)
				{
					if (k < gp.num_env_qtl)
						x = x + (Y[k].per_locus[kk] * in[k].per_locus[kk] * xlast) + env_cue;
					else
						x = x + (Y[k].per_locus[kk] * in[k].per_locus[kk] * xlast);
				}
				x = 1 / (1 + exp(-1 * x));
				diff = x - xlast;
				xlast = x;
			}
			if (stability)
			{
				if (abs(diff) > 0.01)
					alive = false;
			}
			if (k >= gp.num_env_qtl)
				X[k - gp.num_env_qtl] = x;
		}
	}
	void initialize_network(parameters gp, vector<tracker>&in,vector<tracker>&Y,vector<double>& x,vector<int>&Z)
	{
		int jj, jjj;
		for (jj = 0; jj < gp.num_qtl; jj++)
		{//the first ones are the environmentally regulated ones.
			in.push_back(tracker());//SxS matrix of gene interactions
			Y.push_back(tracker());	//SxS matrix of whether genes regulate each other					
			if (jj >= gp.num_env_qtl)
			{
				Z.push_back(1);		//(S-E) vector of whether gene j contributes to trait i. 
				x.push_back(double());//(S-E) vector of expression levels for genes.
			}
			for (jjj = 0; jjj < gp.num_qtl; jjj++)
			{
				Y[jj].per_locus.push_back(1);
				in[jj].per_locus.push_back(randnorm(0, 0.5));
			}
		}
	}
	//overloaded phenotype functions
	void determine_threshold(parameters gp, double cthresh, double pthresh)
	{
		int k,kk;
		if (gp.thresholds_evolve)
		{
			ind_cthresh = ind_pthresh = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
				{
					if (gp.courter_conditional || gp.court_trait)
						ind_cthresh = ind_cthresh + maternal[k].courter_thresh[kk] + paternal[k].courter_thresh[kk];
					if (gp.parent_conditional || gp.parent_trait)
						ind_pthresh = ind_pthresh + maternal[k].parent_thresh[kk] + paternal[k].parent_thresh[kk];
				}
			}
		}
		else
		{
			if (gp.courter_conditional || gp.court_trait)
				ind_cthresh = cthresh;
			if (gp.parent_conditional || gp.parent_trait)
				ind_pthresh = pthresh;
		}
	}
	void determine_threshold(parameters gp, double cthresh, double pthresh, vector<tracker>& cthresh_env, vector<tracker>& pthresh_env, double env_cue = 0, bool stability = true)
	{
		int k, kk;
		if (gp.thresholds_evolve)//sanity check
		{
			if(!gp.gene_network)//again, sanity check
			{
				determine_threshold(gp, cthresh, pthresh);
			}
			else
			{
				if (gp.parent_trait || gp.parent_conditional)
				{
					ind_pthresh = 0;
					gene_interactions(gp, pthresh_Y, pthresh_int, pthresh_x, env_cue,stability);
					int qtl_index = 0;
					for (k = 0; k < gp.num_chrom; k++)
					{
						for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
						{
							if (pthresh_env[k].per_locus[kk] < 0)//if it's not an environmental qtl
								ind_pthresh = ind_pthresh + (pthresh_Z[qtl_index] * (maternal[k].parent_thresh[kk] + paternal[k].parent_thresh[kk])*pthresh_x[qtl_index]);
							qtl_index++;
						}
					}
				}
				if (gp.courter_conditional || gp.court_trait)
				{
					ind_cthresh = 0;
					gene_interactions(gp, cthresh_Y, cthresh_int, cthresh_x, env_cue);
					int qtl_index = 0;
					for (k = 0; k < gp.num_chrom; k++)
					{
						for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
						{
							if (cthresh_env[k].per_locus[kk] < 0)//if it's not an environmental qtl
								ind_cthresh = ind_cthresh + (cthresh_Z[qtl_index] * (maternal[k].courter_thresh[kk] + paternal[k].courter_thresh[kk])*cthresh_x[qtl_index]);
							qtl_index++;
						}
					}
				}
			}
		}
		else
		{
			if (gp.courter_conditional || gp.court_trait)
				ind_cthresh = cthresh;
			if (gp.parent_conditional || gp.parent_trait)
				ind_pthresh = pthresh;
		}
	}
	void calc_courter_trait(parameters gp)//for non-env effects
	{
		int k, kk;
		courter_trait = 0;
		if (!gp.gene_network)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
					courter_trait = courter_trait + maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk];
			}
		}
		else
		{
			std::cout << "\nWARNING: Courter traits not calculated because environmentally responsive QTL are turned on but not specified in function calls.\n";
		}
	}
	void calc_courter_trait(parameters gp, vector<tracker>& courter_env, double env_cue = 0, bool stability = true)
	{
		int k, kk;
		courter_trait = 0;
		if (!gp.gene_network)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
					courter_trait = courter_trait + maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk];
			}
		}
		else
		{
			gene_interactions(gp, courter_Y, courter_int,courter_x, env_cue, stability);
			int trait_qtl_index= 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)//loop through courter_env_qtls
				{
					if (courter_env[k].per_locus[kk] >= gp.num_env_qtl)//if it's not an env qtl
					{
						courter_trait = courter_trait + (courter_Z[trait_qtl_index] * (maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk])*courter_x[trait_qtl_index]);
						trait_qtl_index++;
					}
				}
			}
		}
	}
	void calc_parent_trait(parameters gp)//for non-environmental effects
	{
		int k, kk;
		parent_trait = 0;
		if (!gp.gene_network)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
					parent_trait = parent_trait + maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk];
			}
		}
		else
		{
			std::cout << "\nWARNING: Parent traits not calculated because environmentally responsive QTL are turned on but not specified in function calls.\n";
		}
	}
	void calc_parent_trait(parameters gp, vector<tracker>& parent_env, double env_cue = 0, bool stability = true)
	{
		int k, kk;
		if(gp.parent_trait)
		{
			parent_trait = 0;
			if (!gp.gene_network)
			{
				for (k = 0; k < gp.num_chrom; k++)
				{
					for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
						parent_trait = parent_trait + maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk];
				}
			}
			else
			{
				gene_interactions(gp, parent_Y, parent_int, parent_x, env_cue, stability);
				int trait_qtl_index = 0;
				for (k = 0; k < gp.num_chrom; k++)
				{
					for (kk = 0; kk < gp.qtl_per_chrom[k]; kk++)
					{
						if (parent_env[k].per_locus[kk] >= gp.num_env_qtl)//if it's not an environmental qtl
						{
							parent_trait = parent_trait + (parent_Z[trait_qtl_index] * (maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk])*parent_x[trait_qtl_index]);
							trait_qtl_index++;
						}
					}
				}
			}
		}
	}
	void calc_preference_trait(parameters gp, double threshold)
	{
		int k, kk;
		female_pref = 0;
		if (!gp.gene_network)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < maternal[k].pref_ae.size(); kk++)
					female_pref = female_pref + maternal[k].pref_ae[kk] + paternal[k].pref_ae[kk];
			}
			pref_trait = female_pref;
		}
		else
		{
			std::cout << "\nWARNING: Preference traits not calculated because environmentally responsive QTL are turned on but not specified in function calls.\n";
		}
		if (female_pref < threshold)
			female_pref = 0;
		else
			female_pref = 1;
	}
	void calc_preference_trait(parameters gp, double threshold, vector<tracker>& pref_env, double env_cue = 0, bool stability = true)
	{
		int k, kk;
		female_pref = 0;
		if (!gp.gene_network)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < maternal[k].pref_ae.size(); kk++)
					female_pref = female_pref + maternal[k].pref_ae[kk] + paternal[k].pref_ae[kk];
			}
			pref_trait = female_pref;
		}
		else
		{
			gene_interactions(gp,pref_Y,pref_int,pref_x, env_cue, stability);
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < maternal[k].pref_ae.size(); kk++)
				{
					if (pref_env[k].per_locus[kk] >= gp.num_env_qtl)//if it's not an environmental qtl
					{
						female_pref = female_pref + (pref_Z[qtl_index] * (maternal[k].pref_ae[kk] + paternal[k].pref_ae[kk])*pref_x[qtl_index]);
						qtl_index++;
					}
				}
			}
			pref_trait = female_pref;
		}
		if (female_pref < threshold)
			female_pref = 0;
		else
			female_pref = 1;
	}
	void assign_court_morph(parameters gp)
	{
		if (courter_trait < ind_cthresh)
			courter = false;
		else
			courter = true;
		if (!female)
		{
			if (courter)
				pot_rs = gp.rs_c;
			else
				pot_rs = gp.rs_nc;
		}
	}
	void assign_parent_morph(parameters gp)
	{
		if (parent_trait < ind_pthresh)
			parent = false;
		else
			parent = true;
		if (!female)
		{
			if (parent)
				pot_rs = gp.rs_p;
			else
				pot_rs = gp.rs_np;
		}
	}
	void assign_conditional_traits(parameters gp)
	{
		if (gp.courter_conditional)
			courter_trait = randnorm(0, gp.allelic_std_dev);
		if (gp.parent_conditional)
			parent_trait = randnorm(0, gp.allelic_std_dev);
		
	}
	void assign_random_prefs()
	{
		if (genrand() > 0.5)
			female_pref = 0;
		else
			female_pref = 1;
	}
	void update_traits(parameters gp, double courter_thresh, double parent_thresh, double pref_thresh)//non-environmental
	{
		if (gp.thresholds_evolve)
			determine_threshold(gp,courter_thresh, parent_thresh);
		if (gp.court_trait)
		{
			calc_courter_trait(gp);
			assign_court_morph(gp);
		}
		if (gp.parent_trait)
		{
			calc_parent_trait(gp);
			assign_parent_morph(gp);
		}
		if (gp.ind_pref || gp.cor_prefs)
			calc_preference_trait(gp, pref_thresh);
		else
		{
			if (!gp.random_mating)
				assign_random_prefs();
		}
	}
	void update_traits(parameters gp, double courter_thresh, double parent_thresh, double pref_thresh, 
		double env_cue, vector<tracker>& court_env,vector<tracker>& parent_env, vector<tracker>& pref_env, vector<tracker>& pthresh_env, vector<tracker>& cthresh_env)//non-environmental
	{
		if (gp.thresholds_evolve)
			determine_threshold(gp,courter_thresh, parent_thresh,cthresh_env,pthresh_env,env_cue);
		if (gp.court_trait)
		{
			calc_courter_trait(gp, court_env, env_cue);
			assign_court_morph(gp);
		}
		if (gp.parent_trait)
		{
			calc_parent_trait(gp,parent_env,env_cue);
			assign_parent_morph(gp);
		}
		if (gp.ind_pref || gp.cor_prefs)
			calc_preference_trait(gp, pref_thresh,pref_env,env_cue);
		else
		{
			if (!gp.random_mating)
				assign_random_prefs();
		}
	}

	//life cycle functions
	int find_qtl_index(parameters gp, int chrom, int qtl_loc)
	{
		int c, qtl_index;
		qtl_index = 0;
		for (c = 0; c < chrom; c++)
		{
			qtl_index = qtl_index + gp.qtl_per_chrom[c];
		}
		qtl_index = qtl_index + qtl_loc;
		return qtl_index;
	}
	void mutation(parameters gp, vector<tracker>& court_qtl, vector<tracker>& parent_qtl, vector<tracker>& pref_qtl,
		vector<tracker>& courter_thresh_qtl, vector<tracker>& parent_thresh_qtl, 
		vector<tracker>& courter_env_qtl, vector<tracker>& parent_env_qtl, vector<tracker>& cthresh_env_qtl, vector<tracker>& pthresh_env_qtl,vector<tracker>&pref_env_qtl)
	{
		int mm, irand, irand2, irand3, qtl_index;
		double rnd1, rnd2;
		double IndMutationRate;
		double MutSD;
		bool mutated;
		double YorZ = ((gp.num_qtl - gp.num_env_qtl)*(gp.num_qtl - gp.num_env_qtl)) / double(gp.num_qtl*gp.num_qtl);
		double alpha, mu_mean, sigma_mu;
		alpha = 0.2;
		mu_mean = 0.02;
		sigma_mu = 0.5; 
		MutSD = sqrt(gp.mutational_var);
		IndMutationRate = gp.mutation_rate * 2 * gp.num_markers*gp.num_chrom;

		rnd1 = genrand();
		mutated = false;
		if (rnd1 < IndMutationRate)
		{
			irand = randnum(gp.num_chrom);
			irand2 = randnum(gp.num_markers);
			rnd2 = genrand();//to choose maternal or paternal
			if (rnd2 < 0.5)//affects maternal chromosome	
			{
				while (!mutated) 
				{
					irand3 = randnum(gp.num_alleles);
					if (!maternal[irand].loci[irand2] == irand3)
					{
						maternal[irand].loci[irand2] = irand3;
						mutated = true;
					}
				}
				for (mm = 0; mm < gp.qtl_per_chrom[irand]; mm++)
				{
					if (gp.court_trait)
					{
						if (court_qtl[irand].per_locus.size() > 0)
						{
							if (court_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, courter_Y, courter_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, courter_Z, courter_env_qtl[irand].per_locus[mm], maternal[irand].courter_ae[mm], alpha, sigma_mu);
								}
								else
									maternal[irand].courter_ae[mm] = maternal[irand].courter_ae[mm] + randnorm(0, MutSD);
							}
						}
					}
					if (gp.parent_trait)
					{
						if (parent_qtl[irand].per_locus.size() > 0)
						{
							if (parent_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, parent_Y, parent_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, parent_Z, parent_env_qtl[irand].per_locus[mm], maternal[irand].parent_ae[mm], alpha, sigma_mu);
								}
								else
									maternal[irand].parent_ae[mm] = maternal[irand].parent_ae[mm] + randnorm(0, MutSD);
							}
						}
					}
					if (gp.thresholds_evolve)
					{
						if (gp.courter_conditional || gp.court_trait)
						{
							if (courter_thresh_qtl[irand].per_locus.size() > 0)
							{
								if (courter_thresh_qtl[irand].per_locus[mm] == irand2)
								{
									if (gp.gene_network)
									{
										qtl_index = find_qtl_index(gp, irand, mm);
										if (genrand() < YorZ)
											mut_env_Y(gp, qtl_index, cthresh_Y, cthresh_int, alpha, sigma_mu);
										else
											mut_env_Z(gp, qtl_index, cthresh_Z, cthresh_env_qtl[irand].per_locus[mm], maternal[irand].courter_thresh[mm], alpha, sigma_mu);
									}
									else
										maternal[irand].courter_thresh[mm] = maternal[irand].courter_thresh[mm] + randnorm(0, MutSD);
								}
							}
						}
						if (gp.parent_conditional || gp.parent_trait)
						{
							if (parent_thresh_qtl[irand].per_locus.size() > 0)
							{
								if (parent_thresh_qtl[irand].per_locus[mm] == irand2)
								{
									if (gp.gene_network)
									{
										qtl_index = find_qtl_index(gp, irand, mm);
										if (genrand() < YorZ)
											mut_env_Y(gp, qtl_index, pthresh_Y, pthresh_int, alpha, sigma_mu);
										else
											mut_env_Z(gp, qtl_index, pthresh_Z, pthresh_env_qtl[irand].per_locus[mm], maternal[irand].parent_thresh[mm], alpha, sigma_mu);
									}
									else
										maternal[irand].parent_thresh[mm] = maternal[irand].parent_thresh[mm] + randnorm(0, MutSD);
								}
							}
						}
					}
				}
				
				for (mm = 0; mm < maternal[irand].pref_ae.size(); mm++)
				{
					if (gp.ind_pref || gp.cor_prefs)
					{
						if (pref_qtl[irand].per_locus.size() > 0)
						{
							if (pref_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, pref_Y, pref_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, pref_Z, pref_env_qtl[irand].per_locus[mm], maternal[irand].pref_ae[mm], alpha, sigma_mu);
								}
								else
									maternal[irand].pref_ae[mm] = maternal[irand].pref_ae[mm] + randnorm(0, MutSD);
							}
						}
					}
				}
			}
			else//affects paternal chromosome
			{
				while (!mutated) 
				{
					irand3 = randnum(gp.num_alleles);
					if (!paternal[irand].loci[irand2] == irand3)
					{
						paternal[irand].loci[irand2] = irand3;
						mutated = true;
					}
				}
				for (mm = 0; mm < gp.qtl_per_chrom[irand]; mm++)
				{
					if (gp.court_trait)
					{
						if (court_qtl[irand].per_locus.size() > 0)
						{
							if (court_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, courter_Y, courter_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, courter_Z, courter_env_qtl[irand].per_locus[mm], paternal[irand].courter_ae[mm], alpha, sigma_mu);
								}
								else
									paternal[irand].courter_ae[mm] = paternal[irand].courter_ae[mm] + randnorm(0, MutSD);
							}
						}						
					}
					if (gp.parent_trait)
					{
						if (parent_qtl[irand].per_locus.size() > 0)
						{
							if (parent_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, parent_Y, parent_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, parent_Z, parent_env_qtl[irand].per_locus[mm], paternal[irand].parent_ae[mm], alpha, sigma_mu);
								}
								else
									paternal[irand].parent_ae[mm] = paternal[irand].parent_ae[mm] + randnorm(0, MutSD);
							}
						}
					}
					if (gp.thresholds_evolve)
					{
						if (gp.courter_conditional || gp.court_trait)
						{
							if (courter_thresh_qtl[irand].per_locus.size() > 0)
							{
								if (courter_thresh_qtl[irand].per_locus[mm] == irand2)
								{
									if (gp.gene_network)
									{
										qtl_index = find_qtl_index(gp, irand, mm);
										if (genrand() < YorZ)
											mut_env_Y(gp, qtl_index, cthresh_Y, cthresh_int, alpha, sigma_mu);
										else
											mut_env_Z(gp, qtl_index, cthresh_Z, cthresh_env_qtl[irand].per_locus[mm], paternal[irand].courter_thresh[mm], alpha, sigma_mu);
									}
									else
										paternal[irand].courter_thresh[mm] = paternal[irand].courter_thresh[mm] + randnorm(0, MutSD);
								}
							}
						}
						if (gp.parent_conditional || gp.parent_trait)
						{
							if (parent_thresh_qtl[irand].per_locus.size() > 0)
							{
								if (parent_thresh_qtl[irand].per_locus[mm] == irand2)
								{
									if (gp.gene_network)
									{
										qtl_index = find_qtl_index(gp, irand, mm);
										if (genrand() < YorZ)
											mut_env_Y(gp, qtl_index, pthresh_Y, pthresh_int, alpha, sigma_mu);
										else
											mut_env_Z(gp, qtl_index, pthresh_Z, pthresh_env_qtl[irand].per_locus[mm], paternal[irand].parent_thresh[mm], alpha, sigma_mu);
									}
									else
										paternal[irand].parent_thresh[mm] = paternal[irand].parent_thresh[mm] + randnorm(0, MutSD);
								}
							}
						}
					}
				}
				
				for (mm = 0; mm < maternal[irand].pref_ae.size(); mm++)
				{
					if (gp.ind_pref || gp.cor_prefs)
					{
						if (pref_qtl[irand].per_locus.size() > 0)
						{
							if (pref_qtl[irand].per_locus[mm] == irand2)
							{
								if (gp.gene_network)
								{
									qtl_index = find_qtl_index(gp, irand, mm);
									if (genrand() < YorZ)
										mut_env_Y(gp, qtl_index, pref_Y, pref_int, alpha, sigma_mu);
									else
										mut_env_Z(gp, qtl_index, pref_Z, pref_env_qtl[irand].per_locus[mm], paternal[irand].pref_ae[mm], alpha, sigma_mu);
								}
								else
									paternal[irand].pref_ae[mm] = paternal[irand].pref_ae[mm] + randnorm(0, MutSD);
							}
						}
					}
				}
			}
		}//end of if	
	}//mutation
	//mutate gene network
	void mut_env_Y(parameters gp, int locus, vector<tracker>&Y, vector<tracker>&in, double alpha, double sigma_mu)
	{
		int rand_loc;
		rand_loc = randnum(gp.num_qtl);//randomly choose an interaction to mutate
		if (genrand() <= alpha)//then it add/removes interactions by mutating Y.
		{
			if (Y[locus].per_locus[rand_loc] == 0)//add interaction
				Y[locus].per_locus[rand_loc] = 1;
			else//remove interaction
				Y[locus].per_locus[rand_loc] = 0;
		}
		else //then it affects weight of interaction between loci
		{
			in[locus].per_locus[rand_loc] = in[locus].per_locus[rand_loc] + randnorm(0, sigma_mu);
		}
	}
	void mut_env_Z(parameters gp,int qtl_index, vector<int> &Z, double env_tracker, double& ae, double alpha, double sigma_mu)
	{
		if (genrand() <= alpha && env_tracker > -1 && qtl_index > gp.num_env_qtl)//then it add/removes interactions by mutating Z.
		{
			if (Z[qtl_index - gp.num_env_qtl] == 0)//add interaction
				Z[qtl_index - gp.num_env_qtl] = 1;
			else//remove interaction
				Z[qtl_index - gp.num_env_qtl] = 0;
		}
		else //then it affects weight of the interaction on the trait
		{
			ae = ae + randnorm(0, sigma_mu);
		}
	}
	
	//DEPRECATED DON'T USE
	//original way of mutating gene network
	int mut_env_Y(parameters gp, vector<tracker>&Y, vector<tracker>&in, double alpha, double sigma_mu)
	{
		int rand_loc1, rand_loc2;
		rand_loc1 = randnum(gp.num_qtl);
		rand_loc2 = randnum(gp.num_qtl);
		if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
		{
			if (Y[rand_loc1].per_locus[rand_loc2] == 0)//add interaction
				Y[rand_loc1].per_locus[rand_loc2] = 1;
			else//remove interaction
				Y[rand_loc1].per_locus[rand_loc2] = 0;
		}
		else //then it affects weight of interaction
		{
			in[rand_loc1].per_locus[rand_loc2] = in[rand_loc1].per_locus[rand_loc2] + randnorm(0, sigma_mu);
		}
	}
	int mut_env_Z(parameters gp, double locus, vector<int> &Z,double maternal_ae, double paternal_ae, double alpha, double sigma_mu)
	{
		if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
		{
			if (Z[locus] == 0)//add interaction
				Z[locus] = 1;
			else//remove interaction
				Z[locus] = 0;
		}
		else //then it affects weight of the interaction on 
		{
			if (genrand() < 0.5)
				maternal_ae = maternal_ae + randnorm(0, sigma_mu);
			else
				paternal_ae = paternal_ae + randnorm(0, sigma_mu);
		}
	}
	void mutation_env(parameters gp)
	{
		//do something.
		int j, jj, rand_loc1, mut_count;
		double alpha, mu_mean, sigma_mu, num_mutations, YorZ;
		alpha = 0.2;
		mu_mean = 0.02;
		sigma_mu = 0.5;
		num_mutations = poissonrand(mu_mean);
		mut_count = 0;
		YorZ = ((gp.num_qtl- gp.num_env_qtl)*(gp.num_qtl- gp.num_env_qtl)) / (gp.num_qtl*gp.num_qtl);
		while (mut_count < num_mutations)//mut_count never gets incremented
		{
			if (gp.court_trait)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					mut_env_Y(gp, courter_Y, courter_int, alpha, sigma_mu);
					mut_count++;
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl - gp.num_env_qtl);
					int index = 0;
					int chrom_id, loc_id;
					for (int j = 0; j < gp.num_chrom; j++)
					{
						for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						{
							if (index == rand_loc1)
							{
								chrom_id = j;
								loc_id = jj;
							}
							index++;
						}
					}
					mut_env_Z(gp, rand_loc1, courter_Z, maternal[chrom_id].courter_ae[loc_id], paternal[chrom_id].courter_ae[loc_id], alpha, sigma_mu);
					mut_count++;
				}
			}
			if (gp.parent_trait)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					mut_env_Y(gp, parent_Y, parent_int, alpha, sigma_mu);
					mut_count++;
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl - gp.num_env_qtl);
					int index = 0;
					int chrom_id, loc_id;
					for (int j = 0; j < gp.num_chrom; j++)
					{
						for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						{
							if (index == rand_loc1)
							{
								chrom_id = j;
								loc_id = jj;
							}
							index++;
						}
					}
					mut_env_Z(gp, rand_loc1, parent_Z, maternal[chrom_id].parent_ae[loc_id], paternal[chrom_id].parent_ae[loc_id], alpha, sigma_mu);
					mut_count++;
				}
			}
			if (gp.ind_pref)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					mut_env_Y(gp, pref_Y, pref_int, alpha, sigma_mu);
					mut_count++;
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl - gp.num_env_qtl);
					int index = 0;
					int chrom_id, loc_id;
					for (int j = 0; j < gp.num_chrom; j++)
					{
						for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
						{
							if (index == rand_loc1)
							{
								chrom_id = j;
								loc_id = jj;
							}
							index++;
						}
					}
					mut_env_Z(gp, rand_loc1, pref_Z, maternal[chrom_id].pref_ae[loc_id], paternal[chrom_id].pref_ae[loc_id], alpha, sigma_mu);
					mut_count++;
				}
			}//ind prefs
			if (gp.thresholds_evolve)
			{
				if (gp.courter_conditional || gp.court_trait)
				{
					if (genrand() >= YorZ)//then we'll be working with Y
					{
						mut_env_Y(gp, cthresh_Y, cthresh_int, alpha, sigma_mu);
						mut_count++;
					}
					else // then we'll work with Z
					{
						rand_loc1 = randnum(gp.num_qtl - gp.num_env_qtl);
						int index = 0;
						int chrom_id, loc_id;
						for (int j = 0; j < gp.num_chrom; j++)
						{
							for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
							{
								if (index == rand_loc1)
								{
									chrom_id = j;
									loc_id = jj;
								}
								index++;
							}
						}
						mut_env_Z(gp, rand_loc1, cthresh_Z, maternal[chrom_id].courter_thresh[loc_id], paternal[chrom_id].courter_thresh[loc_id], alpha, sigma_mu);
						mut_count++;
					}
				}
				if (gp.parent_conditional || gp.parent_trait)
				{
					if (genrand() >= YorZ)//then we'll be working with Y
					{
						mut_env_Y(gp, pthresh_Y, pthresh_int, alpha, sigma_mu);
						mut_count++;
					}
					else // then we'll work with Z
					{
						rand_loc1 = randnum(gp.num_qtl - gp.num_env_qtl);
						int index = 0;
						int chrom_id, loc_id;
						for (int j = 0; j < gp.num_chrom; j++)
						{
							for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
							{
								if (index == rand_loc1)
								{
									chrom_id = j;
									loc_id = jj;
								}
								index++;
							}
						}
						mut_env_Z(gp, rand_loc1, pthresh_Z, maternal[chrom_id].parent_thresh[loc_id], paternal[chrom_id].parent_thresh[loc_id], alpha, sigma_mu);
						mut_count++;
					}
				}
			}
			if (gp.cor_prefs && gp.court_trait)//it's the same as courter
			{
				for (j = 0; j < gp.num_qtl; j++)
				{
					for (jj = 0; jj < gp.num_qtl; jj++)
					{
						pref_Y[j].per_locus[jj] = courter_Y[j].per_locus[jj];
						pref_int[j].per_locus[jj] = courter_int[j].per_locus[jj];
					}
				}
				for (j = 0; j < gp.num_chrom; j++)
				{
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
					{
						maternal[j].pref_ae[jj] = maternal[j].courter_ae[jj];
						paternal[j].pref_ae[jj] = paternal[j].courter_ae[jj];
					}
				}
			}
			if (gp.cor_prefs && !gp.court_trait)
			{
				for (j = 0; j < gp.num_qtl; j++)
				{
					for (jj = 0; jj <  gp.num_qtl; jj++)
					{
						pref_Y[j].per_locus[jj] = parent_Y[j].per_locus[jj];
						pref_int[j].per_locus[jj] = parent_int[j].per_locus[jj];
					}
				}
				for (j = 0; j < gp.num_chrom; j++)
				{
					for (jj = 0; jj < maternal[j].pref_ae.size(); jj++)
					{
						maternal[j].pref_ae[jj] = maternal[j].parent_ae[jj];
						paternal[j].pref_ae[jj] = paternal[j].parent_ae[jj];
					}
				}
			}
		}
	}
};
