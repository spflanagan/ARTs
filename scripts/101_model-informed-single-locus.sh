#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###
# simply run this script - it reliees on gnu parallel to spawn these jobs



###----DETERMINE WHAT SHOULD RUN----###
NUMREPS=5

#Parameters of interest include 
## for high diversity
#-crs
CRS=5
#-ncrs
NRS=7
#-sperm-r
C=0.5
## for low diversity
#-crs
LCRS=6
#-ncrs
LNRS=6
#-sperm-r
LC=0.5


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR


### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do

    # expected high diversity
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_monogamy_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_monogamy_nm_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_nm_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_RM_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_monogamy_RM_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_v0_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_highDiversity_monogamy_v0_allSneak_${i} --all-sneak -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4" >> "parallel_cmds.sh"
	
	# expected low diversity
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_monogamy_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_monogamy_nm_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_nm_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_RM_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_monogamy_RM_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_v0_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_monogamy_v0_allSneak_${i} --all-sneak -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4" >> "parallel_cmds.sh"

done

parallel -k < "parallel_cmds.sh"

