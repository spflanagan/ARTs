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
CRS=4
#-ncrs
NRS=8
#-sperm-r
C=0.5
## for low diversity
#-crs
LCRS=8
#-ncrs
LNRS=4
#-sperm-r
LC=0.5


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

rm parallel_cmds.sh
### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do

    # expected high diversity
		echo "./ARTs --courter --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_qtls_${i} -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -surv-noparent 0 -surv-parent 1 --viability --same-base -mm 4 -p 4 --polygyny" >> "parallel_cmds.sh"
		echo "./ARTs --courter --parent -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_nm_qtls_${i} -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -surv-noparent 0 -surv-parent 1 --viability --same-base -mm 4 -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
		
		echo "./ARTs --courter --parent --supergene -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_supergene_${i} -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -surv-noparent 0 -surv-parent 1 --viability --same-base -mm 4 -p 4 --polygyny" >> "parallel_cmds.sh"
		echo "./ARTs --courter --parent --supergene -b ../../results/single_locus/pcu_1locus_highDiversity_polygyny_nm_supergene_${i} -crs ${CRS} -ncrs ${NRS} -sperm-r ${C} -surv-noparent 0 -surv-parent 1 --viability --same-base -mm 4 -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	
	
	# expected low diversity
	
	echo "./ARTs --courter --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_qtls_${i} -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -surv-noparent 0 -surv-parent 1 --viability -mm 4 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --parent -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_nm_qtls_${i} -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -surv-noparent 0 -surv-parent 1 --viability -mm 4 --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	
	echo "./ARTs --courter --parent --supergene -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_supergene_${i} -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -surv-noparent 0 -surv-parent 1 --viability -mm 4 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --parent --supergene -b ../../results/single_locus/pcu_1locus_lowDiversity_polygyny_nm_supergene_${i} -crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} -surv-noparent 0 -surv-parent 1 --viability -mm 4 --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	

done

parallel -k < "parallel_cmds.sh"

