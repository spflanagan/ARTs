#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.

### --- TO RUN THE SCRIPT --- ###
# simply run this script - it reliees on gnu parallel to spawn these jobs



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


### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do

    # with polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4 --polygyny " >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	
	# without polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --same-base -p 4" >> "parallel_cmds.sh"

    # with polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nm_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nm_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nm_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nm_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nm_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --same-base -p 4 --polygyny --allow-no-mating" >> "parallel_cmds.sh"
	
	# without polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_nm_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_nm_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_nm_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_nm_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_nm_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --same-base -p 4 --allow-no-mating" >> "parallel_cmds.sh"

	# random mating, polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_RM_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --polygyny --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_RM_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --random-mating --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_RM_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base --random-mating -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_RM_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base --random-mating -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_RM_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --random-mating --same-base -p 4 --polygyny" >> "parallel_cmds.sh"

	# random mating without polygyny
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_RM_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_RM_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 --viability --same-base -p 4--random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_RM_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_RM_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_RM_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 --viability --same-base -p 4 --random-mating" >> "parallel_cmds.sh"

	# with polygyny, no viability
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 -v 0 --same-base -p 4 --polygyny " >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_polygyny_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 -v 0 --same-base -p 4 --polygyny" >> "parallel_cmds.sh"
	
	# without polygyny, no viability selection
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -v 0 --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_equalFitness_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -sperm-r 1 -v 0 --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_onlySpermComp_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -v 0 --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 4 -ncrs 4 -v 0 --same-base -p 4" >> "parallel_cmds.sh"
	echo "./ARTs --courter --no-genetics --parent -b ../../results/single_locus/pcu_1locus_nestBinary_sneakRSadvantage_courterRSadvantage_${i} --verbose -q 1 -x 1 -c 1 -surv-noparent 0 -surv-parent 1 -crs 8 -ncrs 4 -v 0 --same-base -p 4" >> "parallel_cmds.sh"


done

parallel -k < "parallel_cmds.sh"

