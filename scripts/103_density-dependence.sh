#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###

# in ARTs/
#  nohup ./scripts/103_density-dependence.sh > ./logs/103_DDMMYY.log 2>&1 &


###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=5

NO_GENETICS=true
GENETIC_ARCH=true
SUPERGENE=true

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
if [ "$GENETIC_ARCH" = true ]; then printf "\t%s\n" "GENETIC_ARCH"; fi
if [ "$SUPERGENE" = true ]; then printf "\t%s\n" "SUPERGENE"; fi
echo "The program will run in the background."
echo "Check the status with htop or by looking at logs/${DATE}.log"
} 


### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do
    #No genetic architectures, just additive genetic variance
    if [ "$NO_GENETICS" = true ]; then
		#with viability selection
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_${i} --verbose --viability --same-base -p 4 &
		#without viability selection
		./ARTs --courter --no-genetics -b ../../results/courter_unlinked_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_novs_${i} --verbose --same-base -p 4 &
    fi


    #with a genetic architecture
    if [ "$GENETIC_ARCH" = true ]; then
		#with viability selection
        ./ARTs --courter -b ../../results/courter_linked_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --parent -b ../../results/parent_linked_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_${i} --verbose --viability --same-base -p 4 &
		#without viability selection
        ./ARTs --courter -b ../../results/courter_linked_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --parent -b ../../results/parent_linked_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_novs_${i} --verbose --same-base -p 4 &
	fi

    #with a supergene genetic architecture
    if [ "$SUPERGENE" = true ]; then
		#with viability selection
        ./ARTs --courter --supergene -b ../../results/courter_supergene_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --parent --supergene -b ../../results/parent_supergene_${i} --verbose --viability --same-base -p 4 &
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_${i} --verbose --viability --same-base -p 4 &
		#without viability selection
        ./ARTs --courter --supergene -b ../../results/courter_supergene_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --parent --supergene -b ../../results/parent_supergene_novs_${i} --verbose --same-base -p 4 &
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_novs_${i} --verbose --same-base -p 4 &
    fi
	
	#wait for all processes in the loop to finish before starting next loop
	wait
done

wait
# generate the report
Rscript -e "rmarkdown::render('../../docs/103_density-dependence.Rmd')" &



