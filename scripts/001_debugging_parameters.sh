#!/bin/bash

## This file runs each of the possible genetic architectures 
## and each of the possible flags
## to identify and remove bugs.

### This script will run the programs in the background and produce a log file in the logs/ directory ###

ARCHS=(	'courter-conditional'
		'parent-conditional'
		'no-genetics'
		'linked-additive'
		'supergene'
		'gene-network'
		'env-cue')
SCENARIOS=( 'freq-dependent-preference'
			'condition-dependent-preference'
			'correlated-pref'
			'independent-pref'
			'thresholds-evolve'
			'thresholds-in-supergene'
			'polygyny'
			'freq-dependent-courter'
			'condition-dependent-courter'
			'freq-dependent-parent'
			'condition-dependent-parent')

###---Move to the correct directory---###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`

###---RUN THE COMBINATIONS---###
for arch in $ARCHS
do
	#run one without any other parameters
	./ARTs --${arch} --courter --parent -b ../../results/001_${arch}

	#run all of the other parameters
	for var in $SCENARIOS
	do	
		./ARTs --${arch} --courter --parent --${var} -b ../../results/001_${var}_${arch}
	done
	
done >> ../../logs/001_${DATE}.log 2>&1 &
