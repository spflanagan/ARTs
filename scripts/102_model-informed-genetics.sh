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
NRS=8
#-sperm-r
C=0.7
## for low diversity
#-crs
LCRS=8
#-ncrs
LNRS=5
#-sperm-r
LC=0.5

# Things to vary
SUPERGENE_PROP_VARS='0.05 0.25 0.5'
NUM_CHROM='2 4 8'
NUM_QTL='8 16 32 64'

### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

rm parallel_cmds.sh
### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do

	for q in $NUM_QTL; do
		
		for c in $NUM_CHROM; do
			# QTLs
			
			# expected high diversity 
			echo "./ARTs --courter --parent \
				-b ../../results/qtls/highDiversity_qtls_q${q}_c${c}_${i} \
				-crs ${CRS} -ncrs ${NRS} -sperm-r ${C} \
				-surv-noparent 0 -surv-parent 1 --viability \
				--same-base -mm 4 -p 4 --polygyny --output-vcf \
				-q ${q} -c ${c}" >> "parallel_cmds.sh"
			echo "./ARTs --courter --parent \
				-b ../../results/qtls/highDiversity_qtls_nm_q${q}_c${c}_${i} \
				-crs ${CRS} -ncrs ${NRS} -sperm-r ${C} \
				-surv-noparent 0 -surv-parent 1 --viability \
				--same-base -mm 4 -p 4 --polygyny --allow-no-mating --output-vcf \
				-q ${q} -c ${c}" >> "parallel_cmds.sh"	
			
			# expected low diversity
			echo "./ARTs --courter --parent \
				-b ../../results/qtls/lowDiversity_polygyny_qtls_q${q}_c${c}_${i} \
				-crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} \
				-surv-noparent 0 -surv-parent 1 --viability \
				-mm 4 --same-base -p 4 --polygyny --output-vcf \
				-q ${q} -c ${c}" >> "parallel_cmds.sh"
			echo "./ARTs --courter --parent \
				-b ../../results/qtls/lowDiversity_qtls_nm_q${q}_c${c}_${i} \
				-crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} \
				-surv-noparent 0 -surv-parent 1 --viability \
				-mm 4 --same-base -p 4 --polygyny --allow-no-mating --output-vcf \
				-q ${q} -c ${c}" >> "parallel_cmds.sh"
				
			# supergenes
			for p in $SUPERGENE_PROP_VARS; do
				
				# expected high diversity 
				echo "./ARTs --courter --parent --supergene \
					-b ../../results/supergene/highDiversity_supergene_prop${p}_q${q}_c${c}_${i} \
					-crs ${CRS} -ncrs ${NRS} -sperm-r ${C} \
					-surv-noparent 0 -surv-parent 1 --viability \
					--same-base -mm 4 -p 4 --polygyny --output-vcf \
					-q ${q} -c ${c} -sprop ${p}" >> "parallel_cmds.sh"
				echo "./ARTs --courter --parent --supergene \
					-b ../../results/supergene/highDiversity_supergene_nm_prop${p}_q${q}_c${c}_${i} \
					-crs ${CRS} -ncrs ${NRS} -sperm-r ${C} \
					-surv-noparent 0 -surv-parent 1 --viability \
					--same-base -mm 4 -p 4 --polygyny --allow-no-mating --output-vcf \
					-q ${q} -c ${c} -sprop ${p}" >> "parallel_cmds.sh"	
				
				# expected low diversity
				echo "./ARTs --courter --parent --supergene \
					-b ../../results/supergene/lowDiversity_supergene_polygyny_prop${p}_q${q}_c${c}_${i} \
					-crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} \
					-surv-noparent 0 -surv-parent 1 --viability \
					-mm 4 --same-base -p 4 --polygyny --output-vcf \
					-q ${q} -c ${c} -sprop ${p}" >> "parallel_cmds.sh"
				echo "./ARTs --courter --parent --supergene \
					-b ../../results/supergene/lowDiversity_supergene_nm_prop${p}_q${q}_c${c}_${i} \
					-crs ${LCRS} -ncrs ${LNRS} -sperm-r ${LC} \
					-surv-noparent 0 -surv-parent 1 --viability \
					-mm 4 --same-base -p 4 --polygyny --allow-no-mating --output-vcf \
					-q ${q} -c ${c} -sprop ${p}" >> "parallel_cmds.sh"
			done
		done
	done
done

parallel -k -j 10 < "parallel_cmds.sh"

