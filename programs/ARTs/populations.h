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
	vector<chromosome> courter_qtls,parent_qtls,pref_qtls,courter_env_qtls,parent_env_qtls,pref_env_qtls, maf, hs;

	population()
	{
		population_size = num_mal = num_fem = num_progeny = sex_trait = max_num_migrants = migrant_index = int();
		theta = mean_fem_traits = mean_mal_traits = vector<double>();
		sex_theta = sex_omega = double();
		adults = progeny = vector<individual>();
		courter_env_qtls = courter_qtls = parent_env_qtls = parent_qtls = pref_env_qtls = pref_qtls = maf = hs = vector<chromosome>();
		sexual_selection = false;
		extinct = false;
		sgenrand(time(0));
	}
};