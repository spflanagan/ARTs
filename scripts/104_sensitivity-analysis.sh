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
POLYGYNY=true
#-surv-parent 0.9
PARENT_SURV=true 
PARENT_SURV_VARS='0.8 0.7 0.6 0.5'
#-surv-noparent 0.1
NONPARENT_SURV=true
NONPARENT_SURV_VARS='0.2 0.3 0.4 0.5'
#-crs 8, -ncrs 4
RS=true
CRS=('4' '8' '8' '2')
NCRS=('8' '8' '2' '8')
#-f 4
#-v 50
#-e 50
#-mm 3
#-sprop 0.1
#-sperm-r 0.5
#-mu 0.0002


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
		#with no genetic architecture
        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
		#with genetic architecture
        ./ARTs --courter -b ../../results/courter_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --parent -b ../../results/parent_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --courter --parent -b ../../results/parent-courter_linked_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        #with supergene
        ./ARTs --courter --supergene -b ../../results/courter_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --parent --supergene -b ../../results/parent_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_polygyny_${i} --verbose --viability --same-base -p 4 --polygyny &
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
	if [ "$PARENT_SURV" = true ]; then
	 	 
    	for j in $PARENT_SURV_VARS
    	do
			#with no genetic architecture
			./ARTs --courter --no-genetics -b ../../results/courter_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent --no-genetics -b ../../results/parent_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --courter -b ../../results/courter_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent -b ../../results/parent_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
			./ARTs --courter --supergene -b ../../results/courter_supergene_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent --supergene -b ../../results/parent_supergene_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			#wait for all processes in the loop statement to finish before starting next
			wait
		done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
	if [ "$NONPARENT_SURV" = true ]; then
		
		for j in $NONPARENT_SURV_VARS
    	do
			#with no genetic architecture
			./ARTs --courter --no-genetics -b ../../results/courter_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent --no-genetics -b ../../results/parent_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --courter -b ../../results/courter_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent -b ../../results/parent_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
			./ARTs --courter --supergene -b ../../results/courter_supergene_npsurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --parent --supergene -b ../../results/parent_supergene_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#wait for all processes in the loop statement to finish before starting next
			wait
		done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
	if [ "$RS" = true ]; then
		for ((idx=0; idx<${#CRS[@]}; ++idx)); do
			crs=${CRS[idx]}
			ncrs=${NCRS[idx]}
			#with no genetic architecture
			./ARTs --courter --no-genetics -b ../../results/courter_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --parent --no-genetics -b ../../results/parent_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --courter -b ../../results/courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --parent -b ../../results/parent_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
			./ARTs --courter --supergene -b ../../results/courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --parent --supergene -b ../../results/parent_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			#wait for all processes in the loop statement to finish before starting next
			wait
		done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
done




