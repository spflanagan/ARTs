# Alternative Reproductive Tactics Models

This repository contains code to model the evolution of alternative reproductive tactics and the code to analyse the outputs of the models.
The repository supports a manuscript submitted to Proceedings of the Royal Society B, which is investigating the effects of explicit genetic architecture on evolutionary dynamics of alternative reproductive tactics. 
The full text of that manuscript is found in `docs/ARTs_ms.Rmd` (which uses `docs/preamble.sty`, `proceedings-of-the-royal-society-b.csl`, and `docs/references.bib`, `R/formatting.R`, and figures in `figs/`). 
A supplementary document, `docs/ARTs_supplement.Rmd`, is also available in this repository.

The `figs/` directory contains an image of the model overview.

Two programs are contained in this repository: a baseline analytical model (written in R) and a simulation-based model (written in C++). I will describe those independently and then describe the contents associated with analysis and writing.

## Baseline analytical model

The baseline model is contained in `morph_predictions/`, which also contains its own README to describe that model. 

Most of the files in that directory are R scripts that run different elements of the baseline model and/or create the shiny interface. One file, `expectations_list.RDS`, contains outputs from the model to rapidly load the shiny interface.

The R scripts have the following functions:

- `app.R`
- `check_freqs.R`
- `expectation_calculations.R`
- `morph_gens.R`
- `morph_gens_ns.R`
- `morph_predictions.R`

## Simulation-based model

The simulation-based model can be found in `programs/ARTs/`, and has the following files in addition to a README for that particular model:

- `lifecycle.cpp`
- `classes.h`
- `populations.h`
- `random_numbers.h`
- `makefile`

Running the simulation model is most efficient on a command-line system, and the scripts used to run it are found in `scripts/`. They are the following:

- `100_single-locus.sh`
- `101_model-informed-single-locus.sh`
- `101_model-informed-single-locus-tradeoffs.sh`
- `102_model-informed-genetics.sh`

## Analysis of the model outputs

The model outputs were all analysed in R. The overall logic and processes for various analyses are found in Rmarkdown documents in `docs/':

- `100_ARTs-notebook.Rmd`
- `100_model-comparisons.Rmd`
- `100_single-locus.Rmd`
- `101_ARTs-genetic-makeups.Rmd`
- `102_genomics-ARTs.Rmd`
- `testing_ternary_plots.Rmd`

These rely on custom functions and other bits of reusable code found in `R/`:

- `001_expectationTests.R`
- `002_freq_functions.R`
- `003_genomic_summaries.R`
- `102_run_genomic_summaries.R`
- `103_run_gwas.R`
- `104_summarize-LD.R`
- `105_summarize-TajimaD.R`
- `formatting.R`

The analysis was conducted in the order of the numbers

## Versions



## Contributors

The work was conceptualised and results interpreted by by Sarah Flanagan and Suzanne Alonzo. All code and analysis was conducted by Sarah Flanagan (spflanagan.phd@gmail.com).

## Funding

This work was conducted in part while SPF was a postdoctoral fellow at the National Institute for Mathematical and Biological Synthesis, sponsored by the National Science Foundation through NSF Award #DBI-1300426, with additional support from The University of Tennessee, Knoxville. 
SPF was partially funded by the Marsden fund grant number UOC1904 while working on this research. 
SHA acknowledges support from funding from the National Science Foundation (Grants IOS-1522491 and IOS-1655297) and the University of California Santa Cruz.

## Licensing

The code is licensed with a GNU GPLv3 license (details in `LICENSE`). Please cite our paper if you use and/or modify the code!

