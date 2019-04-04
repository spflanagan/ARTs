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
#-surv-parent 0.9 -surv-noparent 0.1
SURVIVAL=false 
PARENT_SURV_VARS=('0.8' '0.7' '0.6' '0.5')
NONPARENT_SURV_VARS=('0.2' '0.3' '0.4' '0.5')
#-crs 8, -ncrs 4
RS=false #goes set by set
CRS=('2' '4' '8')
NCRS=('2' '4' '8')
# Combinations of both research allocation types!
RESEARCH_ALLOCATION=true #only runs parent-courter combos
#-f 4
FECUNDITY=false
FECUNDITY_VARS='2 8'
#-v 50
VIABILITY=false
VIABILITY_VARS='25 75 100'
#-e 50
ENCOUNTERS=false
ENCOUNTERS_VARS='25 75 100'
#-mm 3
MAXMATES=false
MAXMATES_VARS='6 12'
#-sprop 0.1
SUPERGENE_PROP=false
SUPERGENE_PROP_VARS='0.05 0.25 0.5'
#-sperm-r 0.5
SPERMR=false
SPERMR_VARS='0.25 0.75'
#-mu 0.0002
MUTATION=false
MU='0.0001 0.0004 0.001'


### --- MOVE TO THE CORRECT DIRECTORIES --- ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR"
PROGDIR="../programs/ARTs"
cd $DIR
cd $PROGDIR

DATE=`date +%Y%m%d`



### --- RUN THE PARAMETER COMBINATIONS --- ###

for i in `seq ${NUMREPS}`; do


	if [ "$RESEARCH_ALLOCATION" = true ]; then
	    for ((par=0; par<${#PARENT_SURV_VARS[@]}; ++par)); do
	        parent=${PARENT_SURV_VARS[par]}
	        for ((np=0; np<${#NONPARENT_SURV_VARS[@]}; ++np)); do
	            nonparent=${NONPARENT_SURV_VARS[np]}
	            for ((c=0; c<${#CRS[@]}; ++c)); do
	                crs=${CRS[c]}
	                for ((nc=0; nc<${#NCRS[@]}; ++nc)); do
	                    ncrs=${NCRS[nc]}
	                    echo "c${crs}_nc${ncrs}_p${parent}_np${nonparent}_${i}"
	                    #with no genetic architecture
	                    ./ARTs --courter --no-genetics --parent -b /data/people/spf50/parent-courter_unlinked_RA_c${crs}_nc${ncrs}_p${parent}_np${nonparent}_${i} -crs ${crs} -ncrs ${ncrs} -surv-parent ${parent} -surv-noparent ${nonparent} --verbose --viability --same-base -p 4
	                    #with genetic architecture
	                    ./ARTs --courter --parent -b /data/people/spf50/parent-courter_linked_RA_c${crs}_nc${ncrs}_p${parent}_np${nonparent}_${i} -crs ${crs} -ncrs ${ncrs} -surv-parent ${parent} -surv-noparent ${nonparent} --verbose --viability --same-base -p 4
	                    #with supergene
	                    ./ARTs --courter --parent --supergene -b /data/people/spf50/parent-courter_supergene_RA_c${crs}_nc${ncrs}_p${parent}_np${nonparent}_${i} -crs ${crs} -ncrs ${ncrs} -surv-parent ${parent} -surv-noparent ${nonparent} --verbose --viability --same-base -p 4
	                    gzip /data/people/spf50/*_RA*.txt
	                done
	            done
	        done
	    done
	fi

    
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
			./ARTs --parent --no-genetics -b ../../results/parent_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --parent -b ../../results/parent_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_psurv${j}_${i} -surv-parent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
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
			./ARTs --parent --no-genetics -b ../../results/parent_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --parent -b ../../results/parent_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_npsurv${j}_${i} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
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
			./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			#with genetic architecture
			./ARTs --courter -b ../../results/courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent -b ../../results/parent-courter_linked_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} -surv-noparent ${j} --verbose --viability --same-base -p 4 &
			#with supergene
			./ARTs --courter --supergene -b ../../results/courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_crs${crs}_ncrs${ncrs}_${i} -crs ${crs} -ncrs ${ncrs} --verbose --viability --same-base -p 4 &
			#wait for all processes in the loop statement to finish before starting next
			wait
		done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait


	if [ "$FECUNDITY" = true ]; then
		for j in $FECUNDITY_VARS
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_fecundity${j}_${i} --verbose --viability --same-base -p 4 -f ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
	if [ "$VIABILITY" = true ]; then
		for j in $VIABILITY_VARS
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_viability${j}_${i} --verbose --viability --same-base -p 4 -v ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait


	if [ "$ENCOUNTERS" = true ]; then
		for j in $ENCOUNTER_VARS
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_encounters${j}_${i} --verbose --viability --same-base -p 4 -e ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait

	if [ "$MAXMATES" = true ]; then
		for j in $MAXMATES_VARS
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_maxmates${j}_${i} --verbose --viability --same-base -p 4 -mm ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait

	if [ "$SUPERGENE_PROP" = true ]; then
		for j in $SUPERGENE_PROP_VARS
    	do
			#with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_prop${j}_${i} --verbose --viability --same-base -p 4 -sprop ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait

	if [ "$SPERMR" = true ]; then
		for j in $SPERMR_VARS
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_spermr${j}_${i} --verbose --viability --same-base -p 4 -sperm-r ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait

	if [ "$MUTATION" = true ]; then
		for j in $MU
    	do
			#with no genetic architecture
	        ./ARTs --courter --no-genetics -b ../../results/courter_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --parent --no-genetics -b ../../results/parent_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --courter --no-genetics --parent -b ../../results/parent-courter_unlinked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
			#with genetic architecture
	        ./ARTs --courter -b ../../results/courter_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --parent -b ../../results/parent_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --courter --parent -b ../../results/parent-courter_linked_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        #with supergene
	        ./ARTs --courter --supergene -b ../../results/courter_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --parent --supergene -b ../../results/parent_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	        ./ARTs --courter --parent --supergene -b ../../results/parent-courter_supergene_mu${j}_${i} --verbose --viability --same-base -p 4 -mu ${j} &
	    done
    fi
	#wait for all processes in the if statement to finish before starting next if statement
	wait
	
done




