#!/bin/bash

# #No genetic architectures
# ./ARTs --courter-conditional -b ../../results/courter-conditional
# ./ARTs --parent-conditional -b ../../results/parent-conditional
# ./ARTs --courter-conditional --parent-conditional -b ../../results/parent-courter-conditional

# #Frequency dependent selection
# ./ARTs --courter-conditional --freq-dependent-preference -b ../../results/courter-conditional_nfds
# ./ARTs --parent-conditional --freq-dependent-preference -b ../../results/parent-conditional_nfds
# ./ARTs --courter-conditional --parent-conditional --freq-dependent-preference -b ../../results/parent-courter-conditional_nfds

# #with a genetic architecture
# ./ARTs --courter -b ../../results/courter
# ./ARTs --parent -b ../../results/parent
#./ARTs --courter --parent -b ../../results/parent-courter

./ARTs --courter --freq-dependent-preference -b ../../results/courter_nfds
./ARTs --parent --freq-dependent-preference -b ../../results/parent_nfds
./ARTs --courter --parent --freq-dependent-preference -b ../../results/parent-courter_nfds

#Evolving thresholds
./ARTs --courter-conditional --thresholds-evolve -b ../../results/courter-conditional_thresholds
./ARTs --parent-conditional --thresholds-evolve -b ../../results/parent-conditional_thresholds
./ARTs --courter-conditional --parent-conditional --thresholds-evolve -b ../../results/parent-courter-conditional_thresholds

./ARTs --courter --thresholds-evolve -b ../../results/courter_thresholds
./ARTs --parent --thresholds-evolve -b ../../results/parent_thresholds
./ARTs --courter --parent --thresholds-evolve -b ../../results/parent-courter_thresholds

