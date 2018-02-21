#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
### and to ensure that all of the processes and outputs are working properly.

#DETERMINE WHAT SHOULD RUN
NO_GENETICS=true
CONDITIONAL=false
COND_NFDS=false
GENETIC_ARCH=false
EVOLVING=false
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/ARTs"

#move to the correct directories - now you can run it from whereever
cd $DIR
cd $PROGDIR


#No genetic architectures, just additive genetic variance
 if [ "$NO_GENETICS" = true ]; then
	./ARTs --courter --no-genetics -b ../../results/courter-nogenetics
    ./ARTs --parent --no-genetics -b ../../results/parent-nogenetics
	./ARTs --courter --no-genetics --parent -b ../../results/parent-courter-nogenetics
	./ARTs --courter --no-genetics --freq-dependent-preference -b ../../results/courter-nogenetics-nfds
	./ARTs --parent --no-genetics --freq-dependent-preference -b ../../results/parent-nogenetics-nfds
	./ARTs --courter --no-genetics --parent --freq-dependent-preference -b ../../results/parent-courter-nogenetics-nfds
fi &

#Random traits
if [ "$CONDITIONAL" = true ]; then
	./ARTs --courter-conditional -b ../../results/courter-conditional
	./ARTs --parent-conditional -b ../../results/parent-conditional
	./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional
fi &

#Frequency dependent selection
if [ "$COND_NFDS" = true ]; then
	./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds
	./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds
	./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds
fi &

#with a genetic architecture
if [ "$GENETIC_ARCH" = true ]; then
	./ARTs --courter -b ../../results/courter
	./ARTs --parent -b ../../results/parent
	./ARTs --courter --parent -b ../../results/parent-courter
	 
	./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds
	./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds
	./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds
fi &

#Evolving thresholds
if [ "$EVOLVING" = true ]; then
	./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds
	./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds
	./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds
	 
	./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds
	./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds
	./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds
fi &