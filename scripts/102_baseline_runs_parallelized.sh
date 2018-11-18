#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###

# in scripts/
#  nohup ./102_baseline_runs.sh > ./logs/002_DATE.log 2>&1 &


###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=5
NUM_THREADS=25
NO_GENETICS=true
CONDITIONAL=false
COND_NFDS=false
GENETIC_ARCH=true
EVOLVING=false
SUPERGENE=true
INDEP_PREF=false
FDS_PREF=false


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
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
echo "Check the status with htop or by looking at logs/002_${DATE}.log"
} | tee ../../logs/002_${DATE}.log 2>1


### --- FIGURE OUT HOW MANY TO RUN --- ###
#NUM_COMMANDS=0
# if [ "$NO_GENETICS" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$CONDITIONAL" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$COND_NFDS" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$GENETIC_ARCH" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$SUPERGENE" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$EVOLVING" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+6)); fi
# if [ "$INDEP_PREF" = true ] && [ "$NO_GENETICS" = true  ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$FDS_PREF" = true ] && [ "$NO_GENETICS" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$INDEP_PREF" = true ] && [ "$GENETIC_ARCH" = true  ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$FDS_PREF" = true ] && [ "$GENETIC_ARCH" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$INDEP_PREF" = true ] && [ "$SUPERGENE" = true  ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi
# if [ "$FDS_PREF" = true ] && [ "$SUPERGENE" = true ]; then NUM_COMMANDS=$((NUM_COMMANDS+3)); fi

N=$((NUM_THREADS))


### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do
    #No genetic architectures, just additive genetic variance
    if [ "$NO_GENETICS" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_novs_${i} --verbose --same-base -p 4 &
        if [ "$INDEP_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics --independent-pref -b ../../results/courter-pref-nogenetics_novs_${i} --verbose --same-base -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --no-genetics --independent-pref -b ../../results/parent-pref-nogenetics_novs_${i} --verbose --same-base -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics --independent-pref --parent -b ../../results/parent-courter-pref-nogenetics_novs_${i} --verbose --same-base -p 4 &
        fi
        if [ "$FDS_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics --freq-dependent-preference -b ../../results/courter-nogenetics-nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --no-genetics --freq-dependent-preference -b ../../results/parent-nogenetics-nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --no-genetics --parent --freq-dependent-preference -b ../../results/parent-courter-nogenetics-nfds_novs_${i} --verbose -p 4 &
        fi
    fi


    #Random traits
    if [ "$CONDITIONAL" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional -b ../../results/courter-conditional_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent-conditional -b ../../results/parent-conditional_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional_novs_${i} --same-base -p 4 &
    fi

    #Frequency dependent selection
    if [ "$COND_NFDS" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds_novs_${i} --same-base -p 4 &
    fi

    #with a genetic architecture
    if [ "$GENETIC_ARCH" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter -b ../../results/courter_linked_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent -b ../../results/parent_linked_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --parent -b ../../results/parent-courter_linked_novs_${i} --verbose --same-base -p 4 &
        if [ "$FDS_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds_novs_${i} --verbose -p 4 &
        fi
        if [ "$INDEP_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --independent-pref -b ../../results/courter-pref_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --independent-pref -b ../../results/parent-pref_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --independent-pref --parent -b ../../results/parent-courter-pref_novs_${i} --verbose -p 4 &
        fi
    fi

    #with a supergene genetic architecture
    if [ "$SUPERGENE" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --supergene -b ../../results/courter_supergene_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --supergene -b ../../results/parent_supergene_novs_${i} --verbose --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_novs_${i} --verbose --same-base -p 4 &
        if [ "$FDS_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --supergene --freq-dependent-preference -b ../../results/courter_supergene_nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --supergene --freq-dependent-preference -b ../../results/parent_supergene_nfds_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --parent --supergene --freq-dependent-preference -b ../../results/parent-courter_supergene_nfds_novs_${i} --verbose -p 4 &
        fi
        if [ "$INDEP_PREF" = true ]; then
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --independent-pref --supergene -b ../../results/courter-pref_supergene_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --independent-pref --supergene -b ../../results/parent-pref_supergene_novs_${i} --verbose -p 4 &
            ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --independent-pref --supergene --parent -b ../../results/parent-courter-pref_supergene_novs_${i} --verbose -p 4 &
        fi  
    fi

    #Evolving thresholds
    if [ "$EVOLVING" = true ]; then
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds_novs_${i} --same-base -p 4 &

        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds_novs_${i} --same-base -p 4 &
        ((ii=ii%N)); ((ii++==0)) && wait; ./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds_novs_${i} --same-base -p 4 &
    fi
done | tee ../../logs/002_${DATE}.log 2>1

wait
# generate the report
Rscript -e "rmarkdown::render('../../docs/002_expectationTests.Rmd')"



