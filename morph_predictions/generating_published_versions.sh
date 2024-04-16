#!/bin/bash

# First the frequency list must be created
Rscript -e create_freqs_list.R

# to run the baseline with 10000 gens and no survival
nohup ./run_morph_gens_parallel.sh "morph_results_Ns_10000.csv" 0 10 > no_survival.log 2>&1 &

# to run with some survival - 10% survive (on the cluster, so j=25)
nohup ./run_morph_gens_parallel.sh "morph_results_Ns_10000_ten-survive.csv" 0.1 25 > some_survival.log 2>&1 &

# to run with more survival - 50% survive (on the cluster, so j=25)
nohup ./run_morph_gens_parallel.sh "morph_results_Ns_10000_fifty-survive.csv" 0.5 25 > more_survival.log 2>&1 &
