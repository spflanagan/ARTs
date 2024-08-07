# Archived outputs from models of alternative reproductive tactics

This zenodo repository (doi: 10.5281/zenodo.8248198) contains outputs from models of the evolution of alternative reproductive tactics, with accompanying code available on github (https://github.com/spflanagan/ARTs).
The repository supports a manuscript submitted to Proceedings of the Royal Society B, which is investigating the effects of explicit genetic architecture on evolutionary dynamics of alternative reproductive tactics. 

The analysis used two separate programs: a baseline analytical model (written in R) and a simulation-based model (written in C++). Outputs from both of these models are archived here.

## Baseline model

- `morph_results_Ns.RDS`: An R data file containing a data.frame with the results of the baseline analytical model. It contains 11 columns:
  - initial_CP = frequency of the courter/parent morph in the initial generation of the model
  - initial_CN = frequency of the courter/non-parent morph in the initial generation of the model
  - initial_NP = frequency of the non-courter/parent morph in the initial generation of the model
  - initial_NN = frequency of the non-courter/non-parent morph in the initial generation of the model
  - CP = frequency of the courter/parent morph in the final generation of the model
  - CN = frequency of the courter/non-parent morph in the final generation of the model
  - NP = frequency of the non-courter/parent morph in the final generation of the model
  - NN = frequency of the non-courter/non-parent morph in the final generation of the model
  - r = the relative reproductive investment parameter
  - c = the sperm competition coefficient
  - num_sneak = the number of males allowed to sneak fertilisations within a single clutch
- `morph_results_10000_equalStart.RDS`: An R data file containing a data.frame with the results of the baseline analytical model after it had been run for 10,000 generations. It contains the same 11 columns as `morph_results_Ns.RDS`.

## Simulation model

The simulation model was run with a variety of parameter combinations, which have been summarised in various files and the archived files are provided.

### Raw outputs

These results were generated by running the scripts `scripts/101_model-informed-single-locus.sh`, `101_model-informed-single-locus-tradeoffs.sh`, and `102_model-informed-genetics.sh`. Each parameter combination was run multiple times and generated the same sets of files:

- `*_parameters.txt`: Outputs the parameter settings for that run in a text file with each parameter on its own line
- `*_log.txt`: A text file with the log outputs from the model - will note whether any errors occurred in that run.
- `*_traits.txt`: A tab-delimited file containing the trait data for every individual in generation 0 and generation 12000. The columns are:
	- Gen: generation
	- Pop: population ID within simulations that started with identical starting conditions
	- Individual: A numerical index for the individual whose information is output.
	- Sex: whether the individual is MALE or FEMALE
	- Courter: A boolean value for whether the individual has the courter trait (1) or the non-courter trait (0). Note that females can carry the courter trait (i.e., have a value of 1 in this column) but do not express the courter trait.
	- CourtTrait: The actual trait value for the individual, which is the sum of allelic effects at courter QTLs. 
	- Parent: A boolean value for whether the individual has the parent trait (1) or the non-parent trait (0). Note that females can carry the parent trait (i.e., have a value of 1 in this column) but do not express the parent trait.
	- ParentTrait: The actual trait value for the individual, which is the sum of allelic effects at parent QTLs.
	- Preference: A boolean value for whether the individual prefers the courting male (1) or non-courting males (0). For the simulation runs here, all individuals have a preference for courters.
	- PrefTrait: The trait value for the preference trait if the trait has a genetic basis. In all iterations of the model shared here, the preference trait was not genetically inherited.
	- MateFound: A count of how many mates the individual was able to obtain.
	- PotRS: The potential reproductive success of the individual, based on their fecundity (which is set by the parameter settings in the model).
	- LifetimeRS: The realised reproductive success of the individual based on the number of matings and the number of offpsring produced.
	- Alive: A boolean value tracking whether the individual died or survived to mate. 
- `*_summary.txt`: A tab-delimited file summarising the final frequencies of various morphs and other demographic parameters for each generation of the model, and for each population (when populations were initiated with identical starting parameters). It contains the following columns:
	- Generation: Generation number (an integer)
	- Pop: Numerical population ID. All populations were initiatlised with identical conditions within a single file.
	- PopSize: The population size (i.e., number of adults)
	- NumMal: Number of adult males in the population
	- NumFem: Number of adult females in the population
	- NumProgeny: The number of progeny produced
	- ParentThresh: The population-level threshold for the parent trait to switch from parent to non-parent (this is the mean allelic effects in Gen 0).
	- ParentFreq: Frequency of the parent trait in the population
	- ParentAEmean: Mean allelic effects of the Parent QTLs
	- ParentAEsd: Standard deviation in allelic effects of the Parent QTLs
	- ParentW: Relative fitness of parent males
	- NonParentW: Relative fitness of non-parent males
	- CourterThresh: The population-level threshold for the courter trait to switch from courter to non-courter (this is the mean allelic effects in Gen 0).
	- CourterFreq: Frequency of the courter trait in the population
	- CourterAEmean: Mean allelic effects of courter QTLs
	- CourterAEsd: Standard deviation of allelic effects of courter QTLs
	- CourterW: Relative fitness of courting males
	- NonCourterW: Relative fitness of non-courting males
	- FreqNcNp: Frequency of non-courting/non-parent (NN) morph
	- FreqCNp: Frequency of courting/non-parent (CN) morph
	- FreqNcP: Frequency of non-courting/parent (NP) morph
	- Freq CP: Frequency of courting/parent (CP) morph
	- PrefThresh: The population-level threshold for the preference trait to switch from preferring the courting to non-courting males (not relevant to these simulations)
	- PrefFreq: Frequency of the preference for the courting male in th epopulation (not relevant to these simulations)
	- NumRandMate: Number of females that randomly mated (i.e., did not find a partner with the preferred trait)

The runs with explicit genetic architectures also have the following files:

- `*_qtlinfo.txt`: A tab-delimited file summarising the location of each type of QTL. Each column is a different QTL and each row is a different population. If they are initialised to be identical, the QTL information will be the same for each population. The format of the QTL location information is a the chromosome number as an integer (starting at 0), followed by a decimal, and the following numbers are the location among the marker loci. So, 0.850 refers to a QTL on chromosome 0 at marker location 850 (out of 1000). 
- `*_allelic-effects.txt`: A tab-delimited file containing the allelic effects for each QTL. These are the additive contributions each QTL makes towards the trait, and these mutate if a mutation occurs at the location of the QTL (which is recorded in the corresponding `*_qtlinfo.txt` file).
- `*_markers.txt`: A tab-delimited file summarising the allele frequency at each marker locus.
- `*vcf`: A variant call format file for the population in the final generation of the simulations. See standard formats for this type of file online (e.g., https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- `*Tajima.D`: The vcftools output format containing Tajima's D statistics, which was generated using the vcf file. See the vcftools manual for more details (https://vcftools.github.io/man_latest.html)
- `*LD.geno.ld`: The vcftools output format containing pairwise linkage disequilibrium statistics, which was generated using the vcf file. See the vcftools manual for more details (https://vcftools.github.io/man_latest.html)
- `*_gt.csv`: A comma-separated file summarising the genotype information for each individual at all marker loci. It is similar to vcf file format, with the following columns:
	- Marker: marker ID
	- Chrom: Chromosome ID or number (starting at 0)
	- Position: Locaiton on the chromosome (starting at 0)
	- REF: Reference allele
	- ALT: Alternative allele
	- The remaining columns are each individual's genotype 
- `*_pheno.csv`: A comma-separated file summarising the phenotypes for a population with the following columns:
	- ID: individual ID
	- CourtTrait: Whether the individual is a courter (2) or a non-courter (1)
	- ParentTrait: Whether the individual is a parent (2) or a non-parent (1)
	- Sex: Whether the individual is a male (MAL) or a female (FEM)
	- Morph: The morph of the individual (CP = courter/parent, C = courter, P = non-courter/parent, N = non-courter/non-parent; females can have preferences as well but this is not relevant to these datasets)


These above outputs, run with various parameter settings, are found in the following `tar.gz` files:

- `stochasticity.tar.gz`: Contains the output for the parameter combinations without explicit genetic architectures but including demographic stochasticity, recombination, and mutation. Their file names contain the following labels:
  - lowDiversity, highDiversity, or highDiversityStrict: indicates the parameter settings used as described in the main text. The 'Strict' case is described in the supplement.
  - polygyny or monogamy: Specifies if males were allowed to mate with multiple females.
  - RM: if present, it means that females mated randomly if no suitable males were found
  - nm: if present, it means that no random mating was permitted if a suitable male was not found
  - v0: if present, no viability selection was imposed on males
  - _\d: the final number in the name represents the replicate number that was run with the particular parameter combination
- `qtls.tar.gz`: Contains the output for parameter combinations run at high diversity, low diversity, various mating systems, and various genetic parameters. Their file names contain the following labels:
  - highDiversity or lowDiversity: indicates the parameter settings used as described in the main text.
  - polygyny or monogamy: Specifies if males were allowed to mate with multiple females. (if this is absent, polygyny is the default)
  - nm: if present, it means that no random mating was permitted if a suitable male was not found
  - q\d+: The number of QTLs underlying each trait
  - c\d+: The number of chromosomes in the simulation
  - _\d+: the final number in the name represents the replicate number that was run with the particular parameter combination
  - pop_\d+: For some genomic analyses, the population number from within a replicate is also included. Within a replicate, all populations had identical starting conditions.  
- `supergene.tar.gz`: Contains the output for parameter combinations run at high diversity, low diversity, various mating systems, and various supergene genetic parameters. Their file names contain the following labels:
  - highDiversity or lowDiversity: indicates the parameter settings used as described in the main text.
  - polygyny or monogamy: Specifies if males were allowed to mate with multiple females. (if this is absent, polygyny is the default)
  - nm: if present, it means that no random mating was permitted if a suitable male was not found
  - prop\d.\d+: specifies the proportion of a chromosome the supergene stretched across. The options are 0.05,0.25,0.5
  - q\d+: The number of QTLs underlying each trait
  - c\d+: The number of chromosomes in the simulation
  - _\d+: the final number in the name represents the replicate number that was run with the particular parameter combination
  - pop_\d+: For some genomic analyses, the population number from within a replicate is also included. Within a replicate, all populations had identical starting conditions.


### Summaries

The outputs in the zipped files listed above were analysed and summarised in a variety of ways, and those summaries were saved into files for ease of analysis and figure making. 

- `lowDiv.RDS`, `highDiv.RDS`, and `highDivStrict.RDS` are three R data files that summarise the final frequencies in the simulations witout explicit genetic archictecures (from the `stochastic.tar.gz` file). These files are produced in `docs/102_all-model-outcomes.Rmd`. All three are lists with names reflecting the replicate ID and population, and each list element contains three frequencies, in this order: 
  - frequency of the courter/parent morph in the final generation of the model
  - frequency of the non-courter/parent morph in the final generation of the model 
  - frequency of the non-courter/non-parent morph in the final generation of the model  
- `morph_freqs_summary.csv`: A comma-separated file containing the final frequencies of all four morphs for all analysed simulations with explicit genetic architectures (qtls vs supergenes). It contains a header row describing the five columns:
  - CP = frequency of the courter/parent morph
  - C = frequency of the courter/non-parent morph
  - P = frequency of the non-courter/parent morph
  - NON = frequency of the non-courter/non-parent morph
  - file = the name and relative path (on my computer) where the results are stored
- `selection_gradients.RDS`: An R data file containing a list of data.frames. Each list item contains the name of the trait file it summarises, which captures details on the simulation settings and replicate number. The summary statistics were generated in `docs/101_explicit-genetics-outcomes.Rmd`. Each data.frame contains one row per population for each of generation 0 and generation 12000, with the following columns:
  - court_cat_b0: intercept of the regression of courter trait values on lifetime reproductive success 
  - court_cat_b0_se: standard error of the intercept of the regression of courter trait values on lifetime reproductive success   
  - court_cat_b0_p: p-value associated with the intercept of the regression of courter trait values on lifetime reproductive success 
  - court_cat_b1: slope of the regression of courter trait values on lifetime reproductive success  
  - court_cat_b1_se: standard error of the slope of the regression of courter trait values on lifetime reproductive success    
  - court_cat_b1_p:  p-value associated with the slope of the regression of courter trait values on lifetime reproductive success 
  - parent_cat_b0: intercept of the regression of parent trait values on lifetime reproductive success  
  - parent_cat_b0_se: standard error of the intercept of the regression of parent trait values on lifetime reproductive success   
  - parent_cat_b0_p: p-value associated with the intercept of the regression of parent trait values on lifetime reproductive success  
  - parent_cat_b1: slope of the regression of parent trait values on lifetime reproductive success  
  - parent_cat_b1_se: standard error of the slope of the regression of parent trait values on lifetime reproductive success     
  - parent_cat_b1_p:  p-value associated with the slope of the regression of parent trait values on lifetime reproductive success 
  - cp_freq: frequency of the courter/parent morph in the final generation of the model
  - np_freq: frequency of the non-courter/parent morph in the final generation of the model 
  - cn_freq: frequency of the courter/non-parent morph in the final generation of the model
  - nn_freq: frequency of the non-courter/non-parent morph in the final generation of the model   
  - gen: which generation the results come from
  - pop: which identically-initialised population doe the results come from
  - supergene: A boolean value representing whether the replicate had a supergene (1) or did not (0)
  - highDiv: Whether the replicate used high-diversity parameter settings (1) or low diversity parameter settings (0)
- `fitness_dat.RDS`: An R data file containing a list of data.frames. Each list item is named as the trait file it summarises, which captures the details of the replicate settings. Each data.frame contains one row per population for generation 0 and generation 12000, with the following columns:
  - cp_muRS: Mean lifetime reproductive success for males of the courter/parent morph
  - np_muRS: Mean lifetime reproductive success for males of the non-courter/parent morph 
  - cn_muRS: Mean lifetime reproductive success for males of the courter/non-parent morph  
  - nn_muRS: Mean lifetime reproductive success for males of the non-courter/non-parent morph   
  - cp_semRS: standard error of the mean lifetime reproductive success for males of the courter/parent morph   
  - np_semRS: standard error of the mean lifetime reproductive success for males of the non-courter/parent morph  
  - cn_semRS: standard error of the mean lifetime reproductive success for males of the courter/non-parent morph    
  - nn_semRS: standard error of the mean lifetime reproductive success for males of the non-courter/non-parent morph      
  - cp_freq: frequency of the courter/parent morph in the final generation of the model
  - np_freq: frequency of the non-courter/parent morph in the final generation of the model 
  - cn_freq: frequency of the courter/non-parent morph in the final generation of the model
  - nn_freq: frequency of the non-courter/non-parent morph in the final generation of the model   
  - gen: which generation the results come from
  - pop: which identically-initialised population doe the results come from
  - supergene: A boolean value representing whether the replicate had a supergene (1) or did not (0)
  - highDiv: Whether the replicate used high-diversity parameter settings (1) or low diversity parameter settings (0)
- `GWAS_summary.csv`: A comma-separated file containing the genome-wide association study summary outputs as generated by `R/103_run_gwas.R`. It contains 9 columns and no header row:
  - Marker ID
  - Chromosome number
  - position on chromosome
  - which GWAS model it was an outlier in
  - The estimated $R^2$
  - the p-value for the GWAS test
  - A boolean value identifying whether the outlier is near a QTL (in a QTL peak, see `docs/103_genomics-ARTs.Rmd`)
  - Which trait was used in the association (Courter or Parent)
  - filename for the vcf file used in the analysis
- `ld_summary_all.csv`: A comma-separated file containing linkage disequilibrium summary statistics. It does not contain a header row or rownames, and has the following columns:
  - Chromosome Id
  - position on the chromosome
  - average $R^2$ value for that locus
  - average $R^2$ value for that locus but only when compared to QTLs
  - average $R^2$ value for that locus but only when compared to non-QTLs
  - filename of LD results (which inludes replicate and population IDs)
- `ld_summary_resampled.csv`: A comma-separated file containing rownames (first column) and 10 additional columns, as well as a header row. The file contains average linkage disequilibrium coefficients and averages from resampling. Thos resampled values were generated using `ld_summary_all.csv` and the analysis is detailed in `docs/103_genomics-ARTs.Rmd`. The rownames are the name of the file analysed (i.e., replicate and population IDs), and the subsequent rows are:
  - nqtl: number of QTL parameterised in the simulation
  - nqtl_actual: the actual number of QTLs in the simulations. This number sometimes differed if the same marker locus was randomly chosen more than once in the initialisation of the simulations (especially in supergene scenarios)
  - nchrom: number of chromosomes
  - prop: the proportion of a chromosome with a supergene (NA if supergenes were not modelled)
  - avgQ: Average pairwise linkage disequilibrium (LD) for QTLs
  - semQ: standard error of the mean pairwise LD for QTLs
  - avgNQ: Average pairwise LD for non-QTLs
  - semNQ: standard error of the mean pairwise LD for non-QTLs
  - avgResampled: average pairwise LD for $n$ randomly chosen loci, where $n$ is the number of QTLs in the simulation, averaged over 999 iterations of random re-sampling
  - semResampled: standard error of the mean pairwise LD for $n$ randomly chosen loci, where $n$ is the number of QTLs in the simulation, averaged over 999 iterations of random re-sampling
- `td_freqs.csv`: A comma-separated file containing summary statistics about Tajima's D values for each population and each replicate. Tajima's D values were estimated by vcftools and the processing of those outputs is recorded in `docs/103_ARTs-genomics`. Note that all items are quoted. The file contains a header row and 10 columns:
  - unnamed (i.e., rownames): Filename and relative path of the replicate and population
  - avgTD: average Tajima's D genome-wide
  - semTD: Standard error of the mean Tajima's D genome-wide
  - avgTD_qtls: average Tajima's D of only the true QTLs
  - semTD_qtls: Standard error of the mean Tajima's D at true QTLs
  - avgTD_neutral: average Tajima's D of only neutral markers
  - semTD_neutral: Standard error of the mean Tajima's D at neutral markers
  - t.t: test statistic from a Welch's corrected t-test comparing average Tajima's D at QTLs vs neutral markers
  - df.df: degrees of freedom from a Welch's corrected t-test comparing average Tajima's D at QTLs vs neutral markers
  - p: p-value from a Welch's corrected t-test comparing average Tajima's D at QTLs vs neutral markers
  - architecture: The genetic architecture in the simulation (QTLs vs supergene)
  - diversity: Whether the simulation was run with high diversity parameters or low diversity parameters
  - CP = frequency of the courter/parent morph in the final generation of the model
  - C = frequency of the courter/non-parent morph in the final generation of the model
  - P = frequency of the non-courter/parent morph in the final generation of the model
  - NON = frequency of the non-courter/non-parent morph in the final generation of the model
  - file = the name and relative path (on my computer) where the results are stored 
- `popgen_summary.csv`: A comma-separated file containing summary statistics from the population genomics analyses for each simulation. The file was generated in `docs/103_genomics-ARTs.Rmd`. It has the following columns:
  - sim: The name of the simulation 
  - nChrom: Number of chromosomes in the simulation
  - nQTLS: Number of QTLs in the simulation
  - propMale: The proportion of individuals in the population that were male
  - propP: Proportion of males in the population that were non-courter/parents
  - propCP: Proportion of males in the population that were courter/parents
  - propC: Proportion of males in the population that were courter/non-parents
  - propN: Proportion of males in the population that were non-courter/non-parents
  - propPref0: Proportion of individuals the preferred the non-courters (not relevant to this study)
  - gHtMean: Mean $H_t$ for all genomic markers
  - gHtMedian: Median $H_t$ for all genomic markers
  - gHtVar: Variance in $H_t$ for all genomic markers
  - qHtMean: Mean $H_t$ for QTLs only 
  - qHtMedian: Median $H_t$ for QTLs only  
  - qHtVar: Variance in $H_t$ for QTLs only 
  - gGpMeanMM: Mean $G'$ for the male-male comparison for all genomic markers
  - gGpMedianMM: Median $G'$ for the male-male comparison for all genomic markers
  - gGpVarMM: Variance in $G'$ for the male-male comparison for all genomic markers
  - qGpMeanMM: Mean $G'$ for the male-male comparison for QTLs only
  - qGpMedianMM: Median $G'$ for the male-male comparison for QTLs only
  - qGpVarMM: Variance in $G'$ for the male-male comparison for QTLs only
  - nPeaksMM: The number of $G'$ peaks identified for the male-male comparison
  - nSigPeaksMM: The number of significant $G'$ peaks in the male-male comparison
  - propQTLsInPeaksMM: Proportion of QTLs that were in $G'$ peaks for the male-male comparison
  - gGpMeanMF: Mean $G'$ for the male-female comparison for all genomic markers
  - gGpMedianMF: Median $G'$ for the male-female comparison for all genomic markers
  - gGpVarMF: Variance in $G'$ for the male-female comparison for all genomic markers
  - qGpMeanMF: Mean $G'$ for the male-female comparison for QTLs only
  - qGpMedianMF: Median $G'$ for the male-female comparison for QTLs only
  - qGpVarMF: Variance in $G'$ for the male-female comparison for QTLs only
  - nPeaksMF: The number of $G'$ peaks identified for the male-female comparison
  - nSigPeaksMF: The number of significant $G'$ peaks in the male-female comparison
  - propQTLsInPeaksMM: Proportion of QTLs that were in $G'$ peaks for the male-male comparison  - Architecture
  - Architecture: The genetic architecture in the simulation (QTLs vs supergene)
  - HetP: The p-value comparing mean heterozygosity between QTLs and genome-wide markers
  - Diversity: Whether the simulation was run with high diversity parameters or low diversity parameters

## Reproducing main figures

To reproduce **Figure 1** in the main text, the following files are required:

- `morph_results_Ns.RDS`
- `lowDiv.RDS`
- `highDiv.RDS`
- `morph_freqs_summary.csv`

These files are converted into appropriate formats for ternary plots and compared to each other in `docs/102_all-model-outcomes.Rmd`.

To reproduce **Figure 2** in the main text, the following files are required:

- `popgen_summary.csv`
- `td_freqs.csv`

The figure is produced in `docs/103_genomics-ARTs`.

