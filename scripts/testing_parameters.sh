#!/bin/bash

./ARTs --courter -b ../../results/testing_courter
./ARTs --courter --parent -b ../../results/testing_both
./ARTs --courter --parent --freq-dependent-preference -b ../../results/testing_both_fdpref
./ARTs --courter --parent --independent-pref -b ../../results/testing_both_prefs
./ARTs --courter --parent --supergene -b ../../results/testing_both_supergene
./ARTs --courter-conditional --thresholds-evolve -b ../../results/testing_nongen_courter_thresh
./ARTs --parent-conditional --thresholds-evolve -b ../../results/testing_nongen_parent_thresh
./ARTs --courter-conditonal --parent-conditional --independent-pref -b ../../results/testing_both_nongen_prefs
./ARTs --courter --parent --thresholds-in-supergene -b ../../results/testing_thresh_supergene
./ARTs --courter --gene-network -b ../../results/testing_courter_network
./ARTs --courter --gene-network --env-cue -b ../../results/testing_courter_enetwork
./ARTs --parent --gene-network -b ../../results/testing_parent_network
./ARTs --parent --gene-network --env-cue -b ../../results/testing_parent_enetwork
./ARTs --courter --parent --independent-pref --gene-network -b ../../results/testing_all_network
./ARTs --courter --parent --independent-pref --gene-network --env-cue -b ../../results/testing_all_enetwork
