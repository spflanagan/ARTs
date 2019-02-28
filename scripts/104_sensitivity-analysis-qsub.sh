#!/bin/bash

## This file runs some initial tests of the ARTs model to verify that the results are in accordance with expectations
## and to ensure that all of the processes and outputs are working properly.
## Additionally, these runs will highlight which parameters play a major role in the model outcomes

### --- TO RUN THE SCRIPT --- ###

#qsub -t 1-NUMREPS -cwd -S /bin/bash 104_sensitivity-analysis-qsub.sh
#where the thread number determines which rep it is


### --- SGE SETTINGS --- ###

#$ -M sarah.flanagan@canterbury.ac.nz
#$ -m abe
#$ -cwd
#$ -o /home/spf50/jobs/
#$ -e /home/spf50/jobs/

###----DETERMINE WHAT SHOULD RUN----###
#NUMREPS=5

#Parameters of interest include 
#--polygyny
POLYGYNY=false
#-surv-parent 0.9
PARENT_SURV=false 
PARENT_SURV_VARS='0.8 0.7 0.6 0.5'
#-surv-noparent 0.1
NONPARENT_SURV=false
NONPARENT_SURV_VARS='0.2 0.3 0.4 0.5'
#-crs 8, -ncrs 4
RS=false
CRS=('4' '8' '8' '2')
NCRS=('8' '8' '2' '8')
#-f 4
FECUNDITY=true
FECUNDITY_VARS='2 8'
#-v 50
VIABILITY=true
VIABILITY_VARS='25 75 100'
#-e 50
ENCOUNTERS=true
ENCOUNTERS_VARS='25 75 100'
#-mm 3
MAXMATES=true
MAXMATES_VARS='6 12'
#-sprop 0.1
SUPERGENE_PROP=true
SUPERGENE_PROP_VARS='0.05 0.25 0.5'
#-sperm-r 0.5
SPERMR=true
SPERMR_VARS='0.25 0.75'
#-mu 0.0002
MUTATION=true
MU='0.0001 0.0004 0.001'


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR




### --- RUN THE PARAMETER COMBINATIONS --- ###

i=$SGE_TASK_ID
    
if [ "$POLYGYNY" = true ]; then
	#with no genetic architecture
    ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
	#with genetic architecture
    ./ARTs --courter -b ../../results/courter_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --parent -b ../../results/parent_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --courter --parent -b ../../results/parent-courter_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    #with supergene
    ./ARTs --courter --supergene -b ../../results/courter_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --parent --supergene -b ../../results/parent_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
    ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny
fi

if [ "$PARENT_SURV" = true ]; then
 	 
	for j in $PARENT_SURV_VARS
	do
		#with no genetic architecture
		./ARTs --parent --no-genetics -b ../../results/parent_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		#with genetic architecture
		./ARTs --parent -b ../../results/parent_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --parent -b ../../results/parent-courter_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		#with supergene
		./ARTs --parent --supergene -b ../../results/parent_supergene_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4
		
	done
fi

if [ "$NONPARENT_SURV" = true ]; then
	
	for j in $NONPARENT_SURV_VARS
	do
		#with no genetic architecture
		./ARTs --parent --no-genetics -b ../../results/parent_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		#with genetic architecture
		./ARTs --parent -b ../../results/parent_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --parent -b ../../results/parent-courter_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		#with supergene
		./ARTs --parent --supergene -b ../../results/parent_supergene_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4
		
	done
fi

if [ "$RS" = true ]; then
	for ((idx=0; idx<${#CRS[@]}; ++idx)); do
		crs=${CRS[idx]}
		ncrs=${NCRS[idx]}
		#with no genetic architecture
		./ARTs --courter --no-genetics -b ../../results/courter_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4
		./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4
		#with genetic architecture
		./ARTs --courter -b ../../results/courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4
		./ARTs --courter --parent -b ../../results/parent-courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} -surv-noparent ${j} --verbose --viability --same-base -p 4
		#with supergene
		./ARTs --courter --supergene -b ../../results/courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4
		./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4
		
	done
fi


if [ "$FECUNDITY" = true ]; then
	for j in $FECUNDITY_VARS
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --parent -b ../../results/parent_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j}
    done
fi

if [ "$VIABILITY" = true ]; then
	for j in $VIABILITY_VARS
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --parent -b ../../results/parent_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j}
    done
fi


if [ "$ENCOUNTERS" = true ]; then
	for j in $ENCOUNTER_VARS
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --parent -b ../../results/parent_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j}
    done
fi

if [ "$MAXMATES" = true ]; then
	for j in $MAXMATES_VARS
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --parent -b ../../results/parent_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j}
    done
fi

if [ "$SUPERGENE_PROP" = true ]; then
	for j in $SUPERGENE_PROP_VARS
	do
		#with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j}
    done
fi

if [ "$SPERMR" = true ]; then
	for j in $SPERMR_VARS
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --parent -b ../../results/parent_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j}
    done
fi

if [ "$MUTATION" = true ]; then
	for j in $MU
	do
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --parent -b ../../results/parent_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --parent --supergene -b ../../results/parent_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j}
    done
fi




