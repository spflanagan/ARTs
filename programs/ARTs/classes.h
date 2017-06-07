#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <ctime>
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
	vector<double> courter_ae, parent_ae, pref_ae;


	chromosome()
	{
		loci = vector<int>();
		courter_ae = parent_ae = pref_ae = vector<double>();
	}


};

class parameters
{
public:
	int carrying_capacity, num_sampled, num_chrom, num_markers, num_qtl, num_env_qtl, max_fecund, max_encounters, num_alleles;
	int num_pops, num_init_gen, num_exp_gen, num_ld_comparisons, rs_c, rs_nc, rs_p, rs_np;
	double mutation_rate, mutational_var, recombination_rate, allelic_std_dev, gaussian_pref_mean, cond_adj, via_sel_strength, supergene_prop;
	string base_name;
	bool court_trait, parent_trait, env_effects, cor_prefs, ind_pref, FD_pref, CD_pref, FD_court, FD_parent,CD_court, CD_parent, polygyny, cor_mal_traits;
	bool supergene, var_recomb;
	vector <int> qtl_per_chrom;

	parameters()
	{
		carrying_capacity = num_sampled = num_chrom = num_markers = num_qtl = max_fecund = max_encounters = num_alleles = int();
		num_pops = num_init_gen = num_exp_gen = int();
		mutation_rate =recombination_rate = allelic_std_dev = double();
		env_effects = court_trait = parent_trait = cor_prefs = ind_pref = FD_pref = CD_pref = FD_court = FD_parent = CD_court = CD_parent = polygyny = cor_mal_traits = supergene =  bool();
		var_recomb = bool();
		qtl_per_chrom = vector<int>();
	}

	void set_defaults()
	{
		num_init_gen = 10000;//10000 normally
		num_exp_gen = 2000;//2000 normally
		num_pops = 1; //default 1
		carrying_capacity = 1000;//1000
		num_sampled = 50;//default 50
		num_chrom = 4;//default 4
		num_markers = 1000;//1000
		num_qtl = 50;//default: 50
		for (int j = 0; j < num_chrom; j++)
		{
			qtl_per_chrom.push_back(num_qtl / num_chrom);//default is an even distribution
		}
		num_env_qtl = 0;//default 0
		max_fecund = 4;//default 4
		max_encounters = 50;//default 50
		num_alleles = 2; //biallelic to start
		mutation_rate = 0.0002;//default 0.0002
		mutational_var = 0;//default 0
		recombination_rate = 0.2;//0.2
		allelic_std_dev = 0.5;//default 0.5
		supergene_prop = 0.1; //default 0.1 (10% of num_markers = 100)
		supergene = var_recomb = false; //default: false
		court_trait = false; //default: false
		parent_trait = false;//default false
		env_effects = false;//default false
		cor_prefs= ind_pref= FD_pref= CD_pref= FD_court= FD_parent= CD_court= CD_parent = false;//no selection
		cor_mal_traits = false;//default false
		polygyny = false;//default false
		base_name = "../../results/arts";//default: "../../results/arts"
		num_ld_comparisons = 100;//default 100
		rs_c = 8;//default 8
		rs_nc = 4;//default 4
		rs_p = 8;//default 8
		rs_np = 4;//default 4
		gaussian_pref_mean = 0;//default 0
		via_sel_strength = 50;//unsure what value to put here (currently 50)
		cond_adj = 0.1;//amount to add/subtract to condition dependent traits (default 0.1)
	}

	void help_message()
	{
		cout << "\n\t\tHELP MENU\n";
		cout << "\nSimulation model of the evolution of alternative reproductive tactics with complex genetic architectures.\n";
		cout << "Below are the parameters to input to the model. (Defaults in parentheses)\n";
		cout << "-b:\tBase for output file names (arts)\n";
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
		cout << "-prs:\tParental male reproductive success (8)\n";
		cout << "-nprs:\tNon-parental male reproductive success (4)\n";
		cout << "-crs:\tCourter male reproductive success (8)\n";
		cout << "-ncrs:\tNon-courter male reproductive success (4)\n";
		cout << "-sprop:\tSupergene proportion of a chromosome (0.1). NOTE: must be > the total number of qtls.\n";
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
		cout << "--random-mating:\tSpecifies no female choice (default: true).\n";
		cout << "--supergene:\tSpecifies whether the QTLs are grouped together in a supergene.\n";
		cout << "--var-recomb:\tVariable recombination rate.\n";
		cout << "-h or --help:\tPrint this help message.\n";
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
			if (tempstring1 == "-h" || tempstring1 == "--help")
			{
				help_message();
				run_program = false;
			}
			else
			{
				run_program = true;
				for (j = 1; j < argc - 1; j++)
				{
					tempstring1 = argv[j];
					if (tempstring1.substr(0, 2) != "--")
					{
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
						if (tempstring1 == "-prs")
							rs_p = atof(tempstring2.c_str());
						if (tempstring1 == "-nprs")
							rs_np = atof(tempstring2.c_str());
						if (tempstring1 == "-crs")
							rs_c = atof(tempstring2.c_str());
						if (tempstring1 == "-ncrs")
							rs_nc = atof(tempstring2.c_str());
						if (tempstring1 == "-sprop")
							supergene_prop = atof(tempstring2.c_str());
					}
					else
					{
						if (tempstring1 == "--plasticity")
							env_effects = true;
						if (tempstring1 == "--freq-dependent-preference")
							FD_pref = ind_pref = true;
						if (tempstring1 == "--condition-dependent-preference")
							CD_pref = ind_pref = true;
						if (tempstring1 == "--courter")
							court_trait = ind_pref = true;
						if (tempstring1 == "--parent")
							parent_trait = true;
						if (tempstring1 == "--freq-dependent-courter")
							FD_court = court_trait = true;
						if (tempstring1 == "--freq-dependent-parent")
							FD_parent = parent_trait = true;
						if (tempstring1 == "--condition-dependent-courter")
							CD_court = court_trait = true;
						if (tempstring1 == "--condition-dependent-parent")
							CD_parent = parent_trait = true;
						if (tempstring1 == "--independent-pref")
							ind_pref = true;
						if (tempstring1 == "--correlated-pref")
							cor_prefs = court_trait = true;
						if (tempstring1 == "--random-mating")
							cor_prefs = ind_pref = false;
						if (tempstring1 == "--supergene")
							supergene = true;
					}
				}
				if (env_effects)
				{
					if (num_env_qtl == 0)
						num_env_qtl = num_qtl / 2;
					else
						num_env_qtl = num_env_qtl / num_chrom;
				}
				if (ind_pref || cor_prefs)
					court_trait = true;
			}
		}
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
		if (env_effects)
			param_out << "\nPlasticity";
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
			
		param_out.close();
	}
};

string determine_date()
{
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__))
//	string date;
	time_t now = time(0);
	tm ltm;
//	localtime_s(&ltm, &now);
//	int yr, mn, dy;
//	string month, day;
//	yr = 1900 + ltm.tm_year;
//	mn = 1 + ltm.tm_mon;
//	month = to_string(mn);
//	if (month.size() == 1)
//		month = "0" + month;
//	dy = ltm.tm_mday;
//	day = to_string(dy);
//	if (day.size() == 1)
//		day = "0" + day;
//	date = to_string(yr) + month + day;
//	return date;
//#else    // LINUX
	const size_t size = 1024;
	char buffer[size];
	std::strftime(buffer, size, "%Y%D%M", &ltm);
	
	return buffer;
//#endif
}

class individual
{
public:
	double courter_trait, parent_trait, female_pref;
	vector<chromosome> maternal, paternal;
	vector<tracker> courter_int, parent_int, pref_int, courter_Y, parent_Y, pref_Y;
	vector<int> courter_Z, parent_Z, pref_Z, courter_x, parent_x, pref_x;
	bool female, alive, courter, parent;
	int mate_found, pot_rs;

	individual()
	{
		courter_trait = parent_trait = female_pref = double();
		maternal = paternal = vector<chromosome>();
		courter_int = parent_int = pref_int = courter_Y = parent_Y = pref_Y = vector<tracker>();
		courter_Z = parent_Z = pref_Z = courter_x = parent_x = pref_x = vector<int>();
		mate_found = pot_rs = int();
		female = alive = parent = courter = bool();
	}

	//phenotype functions
	void courter_gene_interactions(parameters gp, double env_cue)
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
						x = x + (courter_Y[k].per_locus[kk] * courter_int[k].per_locus[kk] * xlast) + env_cue;
					else
						x = x + (courter_Y[k].per_locus[kk] * courter_int[k].per_locus[kk] * xlast);
				}
				x = 1 / (1 + exp(-1 * x));
				diff = x - xlast;

			}
			if (diff <= 0.01)
				alive = false;
			if (k > gp.num_env_qtl)
				courter_x[k] = x;
		}
	}
	void parent_gene_interactions(parameters gp, double env_cue)
	{
		int k, kk, t;
		double diff, x, xlast;

		for (k = 0; k < gp.num_qtl; k++)
		{
			xlast = 1;
			x = 0;
			for (t = 0; t < 20; t++)
			{
				for (kk = 0; kk < gp.num_qtl; kk++)
				{
					if (k < gp.num_env_qtl)
						x = x + (parent_Y[k].per_locus[kk] * parent_int[k].per_locus[kk] * xlast) + env_cue;
					else
						x = x + (parent_Y[k].per_locus[kk] * parent_int[k].per_locus[kk] * xlast);
				}
				x = 1 / (1 + exp(-1 * x));
				diff = x - xlast;
			}
			if (diff <= 0.01)
				alive = false;
			if (k > gp.num_env_qtl)
				parent_x[k] = x;
		}
	}
	void pref_gene_interactions(parameters gp, double env_cue)
	{
		int k, kk, t;
		double diff, x, xlast;

		for (k = 0; k < gp.num_qtl; k++)
		{
			xlast = 1;
			x = 0;
			for (t = 0; t < 20; t++)
			{
				for (kk = 0; kk <gp.num_qtl; kk++)
				{
					if (k < gp.num_env_qtl)
						x = x + (pref_Y[k].per_locus[kk] * pref_int[k].per_locus[kk] * xlast) + env_cue;
					else
						x = x + (pref_Y[k].per_locus[kk] * pref_int[k].per_locus[kk] * xlast);
				}
				x = 1 / (1 + exp(-1 * x));
				diff = x - xlast;
			}
			if (diff <= 0.01)
				alive = false;
			if (k > gp.num_env_qtl)
				parent_x[k] = x;
		}
	}
	void calc_courter_trait(parameters gp, double env_cue = 0)
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
			courter_gene_interactions(gp, env_cue);
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = gp.num_env_qtl; kk < gp.num_qtl; kk++)
				{
					courter_trait = courter_trait + (courter_Z[qtl_index] * (maternal[k].courter_ae[kk] + paternal[k].courter_ae[kk])*courter_x[qtl_index]);
					qtl_index++;
				}
			}
		}
	}
	void calc_parent_trait(parameters gp, double env_cue = 0)
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
			parent_gene_interactions(gp, env_cue);
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = gp.num_env_qtl; kk < gp.num_qtl; kk++)
				{
					parent_trait = parent_trait + (parent_Z[qtl_index] * (maternal[k].parent_ae[kk] + paternal[k].parent_ae[kk])*parent_x[qtl_index]);
					qtl_index++;
				}
			}
		}
	}
	void calc_preference_trait(parameters gp, double threshold, double env_cue = 0)
	{
		int k, kk, kkk;
		female_pref = 0;
		if (!gp.env_effects)
		{
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = 0; kk < gp.num_alleles; kk++)
					female_pref = female_pref + maternal[k].pref_ae[kk] + paternal[k].pref_ae[kk];
			}
		}
		else
		{
			pref_gene_interactions(gp, env_cue);
			int qtl_index = 0;
			for (k = 0; k < gp.num_chrom; k++)
			{
				for (kk = gp.num_env_qtl; kk < gp.num_qtl; kk++)
				{
					female_pref = female_pref + (pref_Z[qtl_index] * (maternal[k].pref_ae[kk] + paternal[k].pref_ae[kk])*pref_x[qtl_index]);
					qtl_index++;
				}
			}
		}
		if (female_pref < threshold)
			female_pref = 0;
		else
			female_pref = 1;
	}
	void assign_court_morph(parameters gp, double threshold)
	{
		if (courter_trait < threshold)
			courter = false;
		else
			courter = true;
		if (courter)
			pot_rs = gp.rs_c;
		else
			pot_rs = gp.rs_nc;
	}
	void assign_parent_morph(parameters gp, double threshold)
	{
		if (parent_trait < threshold)
			parent = false;
		else
			parent = true;
		if (parent)
			pot_rs = gp.rs_p;
		else
			pot_rs = gp.rs_np;
	}
	void update_traits(parameters gp, double court_thresh, double parent_thresh, double env_cue = 0)
	{
		if (gp.court_trait)
		{
			calc_courter_trait(gp, env_cue);
			assign_court_morph(gp, court_thresh);
		}
		if (gp.parent_trait)
		{
			calc_parent_trait(gp, env_cue);
			assign_parent_morph(gp, parent_thresh);
		}
		if (gp.ind_pref || gp.cor_prefs)
			calc_preference_trait(gp, env_cue);
	}

	//life cycle functions
	void mutation(parameters gp, vector<tracker>& court_qtl, vector<tracker>& parent_qtl, vector<tracker>& pref_qtl)
	{
		int m, mm, irand, irand2, gg, ggg, irand3, locus;
		double rnd1, rnd2;
		double IndMutationRate;
		double MutSD;
		bool mutated;

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
				while (!mutated) {
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
						if (court_qtl[irand].per_locus[mm] == irand2)
							maternal[irand].courter_ae[mm] =
							maternal[irand].courter_ae[mm] + randnorm(0, MutSD);
					}
					if (gp.parent_trait)
					{
						if (parent_qtl[irand].per_locus[mm] == irand2)
							maternal[irand].parent_ae[mm] =
							maternal[irand].parent_ae[mm] + randnorm(0, MutSD);
					}
					if (gp.ind_pref || gp.cor_prefs)
					{
						if (pref_qtl[irand].per_locus[mm] == irand2)
							maternal[irand].pref_ae[mm] =
							maternal[irand].pref_ae[mm] + randnorm(0, MutSD);
					}
				}
			}
			else//affects paternal chromosome
			{
				while (!mutated) {
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
						if (court_qtl[irand].per_locus[mm] == irand2)
							paternal[irand].courter_ae[mm] =
							paternal[irand].courter_ae[mm] + randnorm(0, MutSD);
					}
					if (gp.parent_trait)
					{
						if (parent_qtl[irand].per_locus[mm] == irand2)
							paternal[irand].parent_ae[mm] =
							paternal[irand].parent_ae[mm] + randnorm(0, MutSD);
					}
					if (gp.ind_pref || gp.cor_prefs)
					{
						if (pref_qtl[irand].per_locus[mm] == irand2)
							paternal[irand].pref_ae[mm] =
							paternal[irand].pref_ae[mm] + randnorm(0, MutSD);
					}
				}
			}
		}//end of if	
	}//mutation
	void mutation_env(parameters gp)
	{
		//do something.
		int j, jj, rand_loc1, rand_loc2, mut_count;
		double alpha, mu_mean, sigma_mu, num_mutations, YorZ;
		alpha = 0.2;
		mu_mean = 0.02;
		sigma_mu = 0.5;
		num_mutations = poissonrand(mu_mean);
		mut_count = 0;
		YorZ = ((gp.num_qtl- gp.num_env_qtl)*(gp.num_qtl- gp.num_env_qtl)) / (gp.num_qtl*gp.num_qtl);
		while (mut_count < num_mutations)
		{
			if (gp.court_trait)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					rand_loc1 = randnum(gp.num_qtl);
					rand_loc2 = randnum(gp.num_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (courter_Y[rand_loc1].per_locus[rand_loc2] == 0)//add interaction
							courter_Y[rand_loc1].per_locus[rand_loc2] = 1;
						else//remove interaction
							courter_Y[rand_loc1].per_locus[rand_loc2] = 0;
					}
					else //then it affects weight of interaction
					{
						courter_int[rand_loc1].per_locus[rand_loc2] = courter_int[rand_loc1].per_locus[rand_loc2] + randnorm(0, sigma_mu);
					}
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl- gp.num_env_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (courter_Z[rand_loc1] == 0)//add interaction
							courter_Z[rand_loc1] = 1;
						else//remove interaction
							courter_Z[rand_loc1] = 0;
					}
					else //then it affects weight of interaction
					{
						int index = 0;
						for (int j = 0; j < gp.num_chrom; j++)
						{
							for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
							{
								if (index == rand_loc1)
								{
									if (genrand() < 0.5)
										maternal[j].courter_ae[jj] = maternal[j].courter_ae[jj] + randnorm(0, sigma_mu);
									else
										paternal[j].courter_ae[jj] = paternal[j].courter_ae[jj] + randnorm(0, sigma_mu);
								}
								index++;
							}
						}
					}
				}
			}
			if (gp.parent_trait)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					rand_loc1 = randnum(gp.num_qtl);
					rand_loc2 = randnum(gp.num_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (parent_Y[rand_loc1].per_locus[rand_loc2] == 0)//add interaction
							parent_Y[rand_loc1].per_locus[rand_loc2] = 1;
						else//remove interaction
							parent_Y[rand_loc1].per_locus[rand_loc2] = 0;
					}
					else //then it affects weight of interaction
					{
						parent_int[rand_loc1].per_locus[rand_loc2] = parent_int[rand_loc1].per_locus[rand_loc2] + randnorm(0, sigma_mu);
					}
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl-gp.num_env_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (parent_Z[rand_loc1] == 0)//add interaction
							parent_Z[rand_loc1] = 1;
						else//remove interaction
							parent_Z[rand_loc1] = 0;
					}
					else //then it affects weight of interaction
					{
						int index = 0;
						for (int j = 0; j < gp.num_chrom; j++)
						{
							for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
							{
								if (index == rand_loc1)
								{
									if (genrand() < 0.5)
										maternal[j].parent_ae[jj] = maternal[j].parent_ae[jj] + randnorm(0, sigma_mu);
									else
										paternal[j].parent_ae[jj] = paternal[j].parent_ae[jj] + randnorm(0, sigma_mu);
								}
								index++;
							}
						}
					}
				}
			}
			if (gp.ind_pref)
			{
				if (genrand() >= YorZ)//then we'll be working with Y
				{
					rand_loc1 = randnum(gp.num_qtl);
					rand_loc2 = randnum(gp.num_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (pref_Y[rand_loc1].per_locus[rand_loc2] == 0)//add interaction
							pref_Y[rand_loc1].per_locus[rand_loc2] = 1;
						else//remove interaction
							pref_Y[rand_loc1].per_locus[rand_loc2] = 0;
					}
					else //then it affects weight of interaction
					{
						pref_int[rand_loc1].per_locus[rand_loc2] = pref_int[rand_loc1].per_locus[rand_loc2] + randnorm(0, sigma_mu);
					}
				}
				else // then we'll work with Z
				{
					rand_loc1 = randnum(gp.num_qtl-gp.num_env_qtl);
					if (genrand() <= alpha)//then it add/removes interactions by mutating Y or Z.
					{
						if (pref_Z[rand_loc1] == 0)//add interaction
							pref_Z[rand_loc1] = 1;
						else//remove interaction
							pref_Z[rand_loc1] = 0;
					}
					else //then it affects weight of interaction
					{
						int index = 0;
						for (int j = 0; j < gp.num_chrom; j++)
						{
							for (int jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
							{
								if (index == rand_loc1)
								{
									if (genrand() < 0.5)
										maternal[j].pref_ae[jj] = maternal[j].pref_ae[jj] + randnorm(0, sigma_mu);
									else
										paternal[j].pref_ae[jj] = paternal[j].pref_ae[jj] + randnorm(0, sigma_mu);
								}
								index++;
							}
						}
					}
				}
			}//ind prefs
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
					for (jj = 0; jj < gp.qtl_per_chrom[j]; jj++)
					{
						maternal[j].pref_ae[jj] = maternal[j].parent_ae[jj];
						paternal[j].pref_ae[jj] = paternal[j].parent_ae[jj];
					}
				}
			}
		}
	}
};