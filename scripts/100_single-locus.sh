#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###

# in ARTs/
#  nohup ./scripts/104_sensitivity-analysis.sh > ./logs/104_DDMMYY.log 2>&1 &


###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=5

#Parameters of interest include 
#--polygyny
POLYGYNY=false


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`



### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do

    # with polygyny
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 4 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4 --polygyny &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 4 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 8 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 8 -crs 8 -ncrs 4 --viability --same-base -p 4 --polygyny &
	# without polygyny
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 4 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4 &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 4 -crs 4 -ncrs 4 --viability --same-base -p 4 &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 8 -crs 4 -ncrs 4 --viability --same-base -p 4 &
	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -prs 4 -nprs 8 -crs 8 -ncrs 4 --viability --same-base -p 4 &
    

done