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

    if [ "$POLYGYNY" = true ]; then
		./ARTs --courter --no-genetics --parent -b ../../results/single_locus/parent-courter_unlinked_1locus_polygyny_${i} --verbose -q 1 -x 1 -c 1 --viability --same-base -p 4 --polygyny &
		
    else
    	./ARTs --courter --no-genetics --parent -b ../../results/single_locus/parent-courter_unlinked_1locus_${i} --verbose -q 1 -x 1 -c 1 --viability --same-base -p 4 &
    fi

done