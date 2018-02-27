#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### This script will run the programs in the background and produce a log file in the logs/ directory ###

###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=10
NO_GENETICS=true
CONDITIONAL=false
COND_NFDS=false
GENETIC_ARCH=false
EVOLVING=false

## move to the correct directories - now you can run it from anywhere ##
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`

### --- RUN THE PARAMETER COMBINATIONS --- ###
echo "Running ${NUMREPS} reps with default parameters with courter, parent, and both traits."
echo "It will run the following scenarios:"
if [ "$NO_GENETICS" = true ]; then printf "\t%s\n" "NO_GENETICS"; fi
if [ "$CONDITIONAL" = true ]; then printf "\t%s\n" "CONDITIONAL"; fi
if [ "$COND_NFDS" = true ]; then printf "\t%s\n" "COND_NFDS"; fi
if [ "$GENETIC_ARCH" = true ]; then printf "\t%s\n" "GENETIC_ARCH"; fi
if [ "$EVOLVING" = true ]; then printf "\t%s\n" "EVOLVING"; fi
echo "The program will run in the background."
echo "Check the status with htop or by looking at logs/002_${DATE}.log"

for i in `seq 1 $NUMREPS`; do

	echo "Starting Rep ${i} of $NUMREPS"
	#No genetic architectures, just additive genetic variance
	 if [ "$NO_GENETICS" = true ]; then
		./ARTs --courter --no-genetics -b ../../results/courter-nogenetics_${i} --verbose
	    ./ARTs --parent --no-genetics -b ../../results/parent-nogenetics_${i} --verbose
		./ARTs --courter --no-genetics --parent -b ../../results/parent-courter-nogenetics_${i} --verbose
		./ARTs --courter --no-genetics --freq-dependent-preference -b ../../results/courter-nogenetics-nfds_${i} --verbose
		./ARTs --parent --no-genetics --freq-dependent-preference -b ../../results/parent-nogenetics-nfds_${i} --verbose
		./ARTs --courter --no-genetics --parent --freq-dependent-preference -b ../../results/parent-courter-nogenetics-nfds_${i} --verbose
	fi 

	#Random traits
	if [ "$CONDITIONAL" = true ]; then
		./ARTs --courter-conditional -b ../../results/courter-conditional_${i}
		./ARTs --parent-conditional -b ../../results/parent-conditional_${i}
		./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional_${i}
	fi

	#Frequency dependent selection
	if [ "$COND_NFDS" = true ]; then
		./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds_${i}
		./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds_${i}
		./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds_${i}
	fi

	#with a genetic architecture
	if [ "$GENETIC_ARCH" = true ]; then
		./ARTs --courter -b ../../results/courter_${i}
		./ARTs --parent -b ../../results/parent_${i}
		./ARTs --courter --parent -b ../../results/parent-courter_${i}
		 
		./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds_${i}
		./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds_${i}
		./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds_${i}
	fi

	#Evolving thresholds
	if [ "$EVOLVING" = true ]; then
		./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds_${i}
		./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds_${i}
		./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds_${i}
		 
		./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds_${i}
		./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds_${i}
		./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds_${i}
	fi
done >> ../../logs/002_${DATE}.log 2>&1 &
