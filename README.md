[![DOI](https://zenodo.org/badge/83730380.svg)](https://zenodo.org/badge/latestdoi/83730380)

# Alternative Reproductive Tactics Models

This repository contains code to model the evolution of alternative reproductive tactics and the code to analyse the outputs of the models.
The repository supports a manuscript submitted to Proceedings of the Royal Society B, which is investigating the effects of explicit genetic architecture on evolutionary dynamics of alternative reproductive tactics. 
The full text of that manuscript is found in `docs/ARTs_ms.Rmd` (which uses `docs/preamble.sty`, `proceedings-of-the-royal-society-b.csl`, and `docs/references.bib`, `R/formatting.R`, and figures in `figs/`). 
A supplementary document, `docs/ARTs_supplement.Rmd`, is also available in this repository.

The `figs/` directory contains an image of the model overview and the figures for the main text.

Relevant model outputs have been archived on zenodo (doi: 10.5281/zenodo.8248198), and the README for those data is in `results/README_results.md`.

Two programs are contained in this repository: a baseline analytical model (written in R) and a simulation-based model (written in C++). I will describe those independently and then describe the contents associated with analysis and writing.

## Baseline analytical model

The baseline model is contained in `morph_predictions/`, which also contains its own README to describe that model. 

Most of the files in that directory are R scripts that run different elements of the baseline model and/or create the shiny interface. One file, `expectations_list.RDS`, contains outputs from the model to rapidly load the shiny interface.

The R scripts have the following functions:

- `app.R`: Creates interactive plots and runs the shiny app
- `check_freqs.R`: This script contains a function to check that the input frequencies provided to the function make sense (i.e., sum to 1)
- `expectation_calculations.R`: This script runs the numerical iterations to generate predictions over some number of generations. It will also plot some outputs if desired.
- `morph_gens_ns.R`: This script contains functions that perform different steps in the life cycle and wrapper scripts to combine those functions in a single generation and then also across multiple generations.
- `morph_predictions.R`: This function generates the predictions for a single generation based on some initial frequencies. 

## Simulation-based model

The simulation-based model can be found in `programs/ARTs/`, and has the following files in addition to a README for that particular model:

- `lifecycle.cpp`: This is where most of the work happens
- `classes.h`: A header file containing class objects (such as 'adults')
- `populations.h`: A header file containing most of the details of life cycle components, as well as initialising functions.
- `random_numbers.h`: A header file containing random number generation functions.

Compilation of the program requires gcc version 11.4.0, and can be done using the included `makefile`.

Running the simulation model is most efficient on a command-line system, and the scripts used to run it are found in `scripts/`. They are the following:

- `100_single-locus.sh`
- `101_model-informed-single-locus.sh`
- `101_model-informed-single-locus-tradeoffs.sh`
- `102_model-informed-genetics.sh`

## Analysis of the model outputs

The model outputs were all analysed in R. The overall logic and processes for various analyses are found in Rmarkdown documents in `docs/':

- `100_baseline-models.Rmd`: Includes details of model output exploration of the baseline analytical model and the simulation model without explicit genetic architectures.
- `101_explicit-genetics-outcomes.Rmd`: Describes exploration of model results from models with either QTLs or supergenes as explicit genetic architectures.
- `102_all-model-outcomes.Rmd`: This document synthesises results that were explored in more depth in `100_baseline-models.Rmd` and `101_explicit-genetics-outcomes`, creating ternary diagrams to summarise results. Figure 1 in the main text is created in this document.
- `103_genomics-ARTs.Rmd`: This document explores population genomic statistics and genome-wide association studies of outputs from the simulation models with explicit genetic architectures. 

These rely on custom functions and other bits of reusable code found in `R/`:

- `formatting.R`: Contains common formatting information for all plots
- `freq_functions.R`: Contains functions to extract frequencies of morphs from raw model output
- `genomic_summaries.R`: Contains functions to summarise genomic data generated from analysis of the vcf files created by the simulation models
- `gwas_preparation.R`: Contains functions to pre-process data to use with the GWAS analysis as well as functions to process GWAS outputs to relate them back to relevant information from the simulation models.
- `102_run_genomic_summaries.R`: Identifies whether QTLs are in outlier peaks for population genetics statistics. Uses functions in `genomic_summaries.R`. 
- `103_run_gwas.R`: Runs genome-wide association studies to identify loci associated with the courter and parent traits. Uses functions in `gwas_preparation.R`.
- `104_summarize-LD.R`: This file compiles information on the linkage disequilibrium analysis and calculates summary statistics for each locus in each simulation.
- `105_summarize-TajimaD.R`: Assesses Tajima D output from vcftools using functions in `genomic_summaries.R`.


The analysis was conducted in the order of the numbers in the names of the files.


## Contributors

The work was conceptualised and results interpreted by by Sarah Flanagan and Suzanne Alonzo. All code and analysis was conducted by Sarah Flanagan (spflanagan.phd@gmail.com).

## Funding

This work was conducted in part while SPF was a postdoctoral fellow at the National Institute for Mathematical and Biological Synthesis, sponsored by the National Science Foundation through NSF Award #DBI-1300426, with additional support from The University of Tennessee, Knoxville. 
SPF was partially funded by the Marsden fund grant number UOC1904 while working on this research. 
SHA acknowledges support from funding from the National Science Foundation (Grants IOS-1522491 and IOS-1655297) and the University of California Santa Cruz.

## Licensing

The code is licensed with a GNU GPLv3 license (details in `LICENSE`). Please cite our paper (once it comes out) if you use and/or modify the code!

