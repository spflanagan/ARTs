#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### This script will run the programs in the background and produce a log file in the logs/ directory ###

###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=10
NO_GENETICS=true
CONDITIONAL=false
COND_NFDS=false
GENETIC_ARCH=true
EVOLVING=false
SUPERGENE=true
INDEP_PREF=false
FDS_PREF=false

## move to the correct directories - now you can run it from anywhere ##
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`

### --- OUTPUT THE CHOICES --- ###
{
echo "Running ${NUMREPS} reps with default parameters with courter, parent, and both traits."
echo "It will run the following scenarios:"
if [ "$NO_GENETICS" = true ]; then printf "\t%s\n" "NO_GENETICS"; fi
if [ "$CONDITIONAL" = true ]; then printf "\t%s\n" "CONDITIONAL"; fi
if [ "$COND_NFDS" = true ]; then printf "\t%s\n" "COND_NFDS"; fi
if [ "$GENETIC_ARCH" = true ]; then printf "\t%s\n" "GENETIC_ARCH"; fi
if [ "$SUPERGENE" = true ]; then printf "\t%s\n" "SUPERGENE"; fi
if [ "$EVOLVING" = true ]; then printf "\t%s\n" "EVOLVING"; fi
if [ "$INDEP_PREF" = true ]; then printf "\t%s\n" "INDEP_PREF"; fi
if [ "$FDS_PREF" = true ]; then printf "\t%s\n" "FDS_PREF"; fi
echo "The program will run in the background."
echo "Check the status with htop or by looking at logs/002_x_${DATE}.log"
} | tee ../../logs/002_${DATE}.log 2>&1

### --- RUN THE PARAMETER COMBINATIONS --- ###
#No genetic architectures, just additive genetic variance
if [ "$NO_GENETICS" = true ]; then
	script1=$(mktemp /tmp/.script.XXXXX)
	cat >$script1 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter --no-genetics -b ../../results/courter_unlinked_${i} --verbose --same-base -p 4
		./ARTs --parent --no-genetics -b ../../results/parent_unlinked_${i} --verbose --same-base -p 4
		./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_${i} --verbose --same-base -p 4
		if [ "$INDEP_PREF" = true ]; then
			./ARTs --courter --no-genetics --independent-pref -b ../../results/courter-pref-nogenetics_${i} --verbose --same-base -p 4
			./ARTs --parent --no-genetics --independent-pref -b ../../results/parent-pref-nogenetics_${i} --verbose --same-base -p 4
			./ARTs --courter --no-genetics --independent-pref --parent -b ../../results/parent-courter-pref-nogenetics_${i} --verbose --same-base -p 4
		fi
		if [ "$FDS_PREF" = true ]; then
			./ARTs --courter --no-genetics --freq-dependent-preference -b ../../results/courter-nogenetics-nfds_${i} --verbose
			./ARTs --parent --no-genetics --freq-dependent-preference -b ../../results/parent-nogenetics-nfds_${i} --verbose
			./ARTs --courter --no-genetics --parent --freq-dependent-preference -b ../../results/parent-courter-nogenetics-nfds_${i} --verbose
		fi
	done >> ../../logs/002_NOGEN_${DATE}.log 2>&1
END
	
else
	script1=`mktemp /tmp/.script.XXXXX`;
	cat >$script1 <<END
	echo "NO_GENETICS not run" >> ../../logs/002_NOGEN_${DATE}.log 2>&1
END
	

fi 

#Random traits
if [ "$CONDITIONAL" = true ]; then
	script2=`mktemp /tmp/.script.XXXXX`;
	cat >$script2 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter-conditional -b ../../results/courter-conditional_${i} --same-base -p 4
		./ARTs --parent-conditional -b ../../results/parent-conditional_${i} --same-base -p 4
		./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional_${i} --same-base -p 4
	done >> ../../logs/002_CONDITIONAL_${DATE}.log 2>&1
END
else
	script2=`mktemp /tmp/.script.XXXXX`;
	cat >$script2 <<END
	echo "CONDITIONAL not run" >> ../../logs/002_CONDITIONAL_${DATE}.log 2>&1
END
fi 

#Frequency dependent selection
if [ "$COND_NFDS" = true ]; then
	script3=`mktemp /tmp/.script.XXXXX`;
	cat >$script3 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds_${i} --same-base -p 4
		./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds_${i} --same-base -p 4
		./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds_${i} --same-base -p 4
	done >> ../../logs/002_CONDNFDS_${DATE}.log 2>&1
END
else
	script3=`mktemp /tmp/.script.XXXXX`;
	cat >$script3 <<END
	echo "COND_NFDS not run" >> ../../logs/002_CONDNFDS_${DATE}.log 2>&1
END
fi 

#with a genetic architecture
if [ "$GENETIC_ARCH" = true ]; then
	script4=`mktemp /tmp/.script.XXXXX`;
	cat >$script4 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter -b ../../results/courter_linked_${i} --verbose --same-base -p 4
		./ARTs --parent -b ../../results/parent_linked_${i} --verbose --same-base -p 4
		./ARTs --courter --parent -b ../../results/parent-courter_linked_${i} --verbose --same-base -p 4
		if [ "$FDS_PREF" = true ]; then
			./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds_${i} --verbose
			./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds_${i} --verbose
			./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds_${i} --verbose
		fi
		if [ "$INDEP_PREF" = true ]; then
			./ARTs --courter --independent-pref -b ../../results/courter-pref_${i} --verbose
			./ARTs --parent --independent-pref -b ../../results/parent-pref_${i} --verbose
			./ARTs --courter --independent-pref --parent -b ../../results/parent-courter-pref_${i} --verbose
		fi
	done >> ../../logs/002_GENETiCS_${DATE}.log 2>&1
END
else
	script4=`mktemp /tmp/.script.XXXXX`;
	cat >$script4 <<END
	echo "GENETIC_ARCH not run" >> ../../logs/002_GENETiCS_${DATE}.log 2>&1
END
fi 

#with a genetic architecture
if [ "$SUPERGENE" = true ]; then
	script5=`mktemp /tmp/.script.XXXXX`;
	cat >$script5 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter --supergene -b ../../results/courter_supergene_${i} --verbose --same-base -p 4
		./ARTs --parent --supergene -b ../../results/parent_supergene_${i} --verbose --same-base -p 4
		./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_${i} --verbose --same-base -p 4
		if [ "$FDS_PREF" = true ]; then
			./ARTs --courter --supergene --freq-dependent-preference -b ../../results/courter_supergene_nfds_${i} --verbose
			./ARTs --parent --supergene --freq-dependent-preference -b ../../results/parent_supergene_nfds_${i} --verbose
			./ARTs --courter --parent --supergene --freq-dependent-preference -b ../../results/parent-courter_supergene_nfds_${i} --verbose
		fi
		if [ "$INDEP_PREF" = true ]; then
			./ARTs --courter --independent-pref --supergene -b ../../results/courter-pref_supergene_${i} --verbose
			./ARTs --parent --independent-pref --supergene -b ../../results/parent-pref_supergene_${i} --verbose
			./ARTs --courter --independent-pref --supergene --parent -b ../../results/parent-courter-pref_supergene_${i} --verbose
		fi
	done >> ../../logs/002_SUPERGENE_${DATE}.log 2>&1
END
else
	script5=`mktemp /tmp/.script.XXXXX`;
	cat >$script5 <<END
	echo "SUPERGENE not run" >> ../../logs/002_SUPERGENE_${DATE}.log 2>&1
END
fi 

#Evolving thresholds
if [ "$EVOLVING" = true ]; then
	script6=`mktemp /tmp/.script.XXXXX`;
	cat >$script6 <<END
	for i in `seq -s '' 1 $NUMREPS`; do
		./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds_${i} --same-base -p 4
		./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds_${i} --same-base -p 4
		./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds_${i} --same-base -p 4
		 
		./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds_${i} --same-base -p 4
		./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds_${i} --same-base -p 4
		./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds_${i} --same-base -p 4
	done >> ../../logs/002_EVOLVING_${DATE}.log 2>&1 
END
else
	script6=`mktemp /tmp/.script.XXXXX`;
	cat >$script6 <<END
	echo "EVOLVING not run" >> ../../logs/002_EVOLVING_${DATE}.log 2>&1
END
fi 


chmod u+rx $script1 $script2 $script3 $script4 $script5 $script6
echo $script1 $script2 $script3 $script4 $script5 $script6
$script1 &!
PID1=$!
$script2 &!
PID2=$!
$script3 &!
PID3=$!
$script4 &!
PID4=$!
$script5 &!
PID5=$!
$script6 &!
PID6=$!
wait $PID1 $PID2 $PID3 $PID4 $PID5 $PID6


#concatenate the log files & remove the intermediate log files
cat ../../logs/002_*_${DATE}.log >> ../../logs/002_${DATE}.log && rm ../../logs/002_*_${DATE}.log
#and the temporary scripts
/bin/rm $script1 $script2 $script3 $script4 $script5 $script6

# generate the report
#R -e "rmarkdown::render('../../scripts/002_expectationTests.Rmd')"



