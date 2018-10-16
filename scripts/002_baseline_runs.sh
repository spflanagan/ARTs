#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###

#  nohup ./scripts/002_baseline_runs.sh > ./logs/002_DATE.log 2>&1 &


###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=5
#NUM_THREADS=8
NO_GENETICS=false
NO_GENETICS=true
CONDITIONAL=false
COND_NFDS=false
GENETIC_ARCH=true
EVOLVING=false
SUPERGENE=true
INDEP_PREF=false
FDS_PREF=false
NUM_COMMANDS=0

### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"  pwd )"
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
} #| tee ../../logs/002_${DATE}.log 2>1


### --- FIGURE OUT HOW MANY TO RUN --- ###
if [ "$NO_GENETICS" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+1)); fi
if [ "$CONDITIONAL" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+1)); fi
if [ "$COND_NFDS" = true ]; then printf NUM_COMMANDS=$((NUM_COMMANDS+1)); fi
if [ "$GENETIC_ARCH" = true ]; then printf NUM_COMMANDS=$((NUM_COMMANDS+1)); fi
if [ "$SUPERGENE" = true ]; then printf NUM_COMMANDS=$((NUM_COMMANDS+1)); fi
if [ "$EVOLVING" = true ]; then printf NUM_COMMANDS=$((NUM_COMMANDS+1)); fi

#NUM_ITER=$((NUM_THREADS/NUMREPS/NUM_COMMANDS))


### --- RUN THE PARAMETER COMBINATIONS --- ###

for iters in `seq 6`; then
#No genetic architectures, just additive genetic variance
if [ "$iters" = 1 ]; then
    if [ "$NO_GENETICS" = true ]; then

        for i in `seq ${NUMREPS}`; do
            ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_${i} --verbose --same-base -p 4
            ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_${i} --verbose --same-base -p 4
            ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_${i} --verbose --same-base -p 4
            if [ "$INDEP_PREF" = true ]; then
                ./ARTs --courter --no-genetics --independent-pref -b ../../results/courter-pref-nogenetics_${i} --verbose --same-base -p 4
                ./ARTs --parent --no-genetics --independent-pref -b ../../results/parent-pref-nogenetics_${i} --verbose --same-base -p 4
                ./ARTs --courter --no-genetics --independent-pref --parent -b ../../results/parent-courter-pref-nogenetics_${i} --verbose --same-base -p 4
            fi
            if [ "$FDS_PREF" = true ]; then
                ./ARTs --courter --no-genetics --freq-dependent-preference -b ../../results/courter-nogenetics-nfds_${i} --verbose -p 4
                ./ARTs --parent --no-genetics --freq-dependent-preference -b ../../results/parent-nogenetics-nfds_${i} --verbose -p 4
                ./ARTs --courter --no-genetics --parent --freq-dependent-preference -b ../../results/parent-courter-nogenetics-nfds_${i} --verbose -p 4
            fi
        done

    fi
fi

    #Random traits
if [ "$iters" = 2 ]; then
    if [ "$CONDITIONAL" = true ]; then

        for i in `seq 1 $NUMREPS`; do
            ./ARTs --courter-conditional -b ../../results/courter-conditional_${i} --same-base -p 4
            ./ARTs --parent-conditional -b ../../results/parent-conditional_${i} --same-base -p 4
            ./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional_${i} --same-base -p 4
        done

    fi
fi
    #Frequency dependent selection
if [ "$iters" = 3 ]; then
    if [ "$COND_NFDS" = true ]; then

        for i in `seq 1 $NUMREPS`; do
            ./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds_${i} --same-base -p 4
            ./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds_${i} --same-base -p 4
            ./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds_${i} --same-base -p 4
        done >> ../../logs/002_CONDNFDS_${DATE}.log 2>1

    fi
fi
    #with a genetic architecture
if [ "$iters" = 4 ]; then
    if [ "$GENETIC_ARCH" = true ]; then

        for i in `seq 1 $NUMREPS`; do
            ./ARTs --courter -b ../../results/courter_linked_${i} --verbose --same-base -p 4
            ./ARTs --parent -b ../../results/parent_linked_${i} --verbose --same-base -p 4
            ./ARTs --courter --parent -b ../../results/parent-courter_linked_${i} --verbose --same-base -p 4
            if [ "$FDS_PREF" = true ]; then
                ./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds_${i} --verbose -p 4
                ./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds_${i} --verbose -p 4
                ./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds_${i} --verbose -p 4
            fi
            if [ "$INDEP_PREF" = true ]; then
                ./ARTs --courter --independent-pref -b ../../results/courter-pref_${i} --verbose -p 4
                ./ARTs --parent --independent-pref -b ../../results/parent-pref_${i} --verbose -p 4
                ./ARTs --courter --independent-pref --parent -b ../../results/parent-courter-pref_${i} --verbose -p 4
            fi
        done

    fi
fi
    #with a genetic architecture
if [ "$iters" = 5 ]; then
    if [ "$SUPERGENE" = true ]; then

        for i in `seq 1 $NUMREPS`; do
            ./ARTs --courter --supergene -b ../../results/courter_supergene_${i} --verbose --same-base -p 4
            ./ARTs --parent --supergene -b ../../results/parent_supergene_${i} --verbose --same-base -p 4
            ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_${i} --verbose --same-base -p 4
            if [ "$FDS_PREF" = true ]; then
                ./ARTs --courter --supergene --freq-dependent-preference -b ../../results/courter_supergene_nfds_${i} --verbose -p 4
                ./ARTs --parent --supergene --freq-dependent-preference -b ../../results/parent_supergene_nfds_${i} --verbose -p 4
                ./ARTs --courter --parent --supergene --freq-dependent-preference -b ../../results/parent-courter_supergene_nfds_${i} --verbose -p 4
            fi
            if [ "$INDEP_PREF" = true ]; then
                ./ARTs --courter --independent-pref --supergene -b ../../results/courter-pref_supergene_${i} --verbose -p 4
                ./ARTs --parent --independent-pref --supergene -b ../../results/parent-pref_supergene_${i} --verbose -p 4
                ./ARTs --courter --independent-pref --supergene --parent -b ../../results/parent-courter-pref_supergene_${i} --verbose -p 4
            fi
        done
    fi
fi
    #Evolving thresholds
if [ "$iters" = 6 ]; then
    if [ "$EVOLVING" = true ]; then

        for i in `seq 1 $NUMREPS`; do
            ./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds_${i} --same-base -p 4
            ./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds_${i} --same-base -p 4
            ./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds_${i} --same-base -p 4

            ./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds_${i} --same-base -p 4
            ./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds_${i} --same-base -p 4
            ./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds_${i} --same-base -p 4
        done
    fi
fi
done

# generate the report
#Rscript -e "rmarkdown::render('../../docs/002_expectationTests.Rmd')"



