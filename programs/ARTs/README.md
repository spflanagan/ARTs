# Alternative Reproductive tactics simulation-based model

This model simulates the evolution of alternative reproductive tactics with explicit genetic architectures. 
It generates a population of adults, who have traits that are determined by additive genetic loci, which are anchored to chromosomes.
These adults mate, produce offspring, and their offspring experience viability selection as juveniles.
The model has non-overlapping generations.

The model has a large number of parameters to set, which can be obtained by reading the help menu (run `./ARTs -h`), which is provided below:


```
		HELP MENU

Simulation model of the evolution of alternative reproductive tactics with complex genetic architectures.
Below are the parameters to input to the model. (Defaults in parentheses)
-b:	Base for output file names (arts)
-K:	carrying capacity (1000)
-s:	number of individuals Sampled. (50)
-c:	number of Chromosomes (4)
-x:	number of markers per chromosome (1000)
-q:	total number of Quantitative trait loci. (50)
-eq:	total number of QTL responding the the Environment (half of -q if --plasticity flag included)
-f:	maximum Fecundity. (4)
-e:	maximum number of Encounters during mating. (50)
-a:	number of alleles (2).
-p:	number of populations. (1)
-i:	number of initial generations. (10000)
-g:	number of experimental generations (2000).
-mu:	maximum mutation rate (0.0002).
-v:	viability selection strength against courters and parents (50).
-r:	Recombination rate. (0.2) 
-asd:	Allelic Standard Deviation (0.5)
-prs:	Parental male reproductive success (only used if no courter trait)
-nprs:	Non-parental male reproductive success (only used if no courter trait)
-crs:	Courter male reproductive success (4)
-ncrs:	Non-courter male reproductive success (8)
-sprop:	Supergene proportion of a chromosome (0.1). NOTE: must be > the total number of qtls.
-sperm-r:	Sperm competition r. If r = 1, paternity is determined through a fair raffle. 
	If 0 < r < 1, the parental male has a higher share of paternity and the sneakers fertilize with penalty r (rs2/(s1+s2); default 0.5)
-surv-noparent:	Survival probability of eggs without a parent (0.1)
-surv-parent:	Survival probability with a parent (0.9)
-mm:	Max number of Mates. This includes the chosen male (3)
-qpc:	QTLs per chromosome (-q/-c, aka an even distribution of QTLs among chromosomes). To specify different numbers per chromosome separate them by commas.
	Example: for 4 chromosomes, input: 4,2,3,0
--gene-network:	Model a plastic morph, where genotype is determined by interactions between genes.
--env-cue:	Model a gene network (--gene-network) that incorporates social information as an environmental cue.
--freq-dependent-preference:	Include this flag if preferences should be based on the frequency of male morphs (defaults to independent).
--condition-dependent-preference:	Include this flag if female preferences are determined by female condition (defaults to independent).
--courter:	Include this flag if males should have the courter trait (a trait affecting mating probabilities). Defaults to preference for randomly chosen morph unless other flags included. 
--parent:	Include this flag if males should have the parental trait (a trait affecting offspring survival). Defaults to preference for randomly chosen morph unless other flags included. 
--freq-dependent-courter:	If the courter trait is experiencing frequency dependent selection.
--freq-dependent-parent:	If the parent trait is experiencing frequency dependent selection.
--condition-dependent-courter:	If the courter trait is influenced by male condition.
--condition-dependent-parent:	If the parent trait is influenced by male condition.
--independent-pref:	Specifies an independent female preference (defaults to Gaussian preference for randomly chosen morph unless other flags included). 
--correlated-pref:	Specifies a female preference correlated with the male courter trait (defaults to Gaussian preference for randomly chosen morph unless other flags included).
--random-mating:	Specifies no female choice (default: true).
--all-sneak:	Specifies all males sneak, not just sneakers (default: false).
--supergene:	Specifies whether the QTLs are grouped together in a supergene that has reduced recombination.
--polygyny:	Allows males to mate multiply (default: true).
--courter-conditional:	If the courter trait has no genetic basis and is determined randomly or through environmental effects.
--parent-conditional:	If the parent trait has no genetic basis and is determined randomly or through environmental effects.
--thresholds-evolve:	If the thresholds are allowed to evolve (i.e., they have a quantitative genetic basis).
--thresholds-in-supergene:	The thresholds have a genetic basis and their loci are in the supergene.
--verbose:	Outputs info about every step during every initial generation -- good for debugging--no-genetics:	Removes genetic architecture; traits encoded by heritable unlinked additive genetic variance.
--linked-additive:	(default) Traits are determined by genome-wide additive genetic variance distributed among chromosomes.
--viability:	If included, viability selection acts on offspring. If not, viability selection is turned off.
--density-independent:	Turns off density-dependent selection. This is not recommended as it will likely lead to population crashes.
--allow-no-mating:	If females cannot find an acceptable male, rather than mate randomly they will not nest at all. (default is false)
--optimize:	Output time steps for initial generations, for optimizing the code. (default is false)
--same-base:	Start each replicate population with the same base population (default is true)
--output-vcf:	Include vcf output for all genotypes of individuals (default is false)
--ae-vcf:	Include vcf output for the allelic effects for QTLs for all individuals in generation 0 (default is false)
--log-file:	Save output to logfile instead of std::cout
--debug:	Output additional information that could be useful in debugging to either log or std::cout
--save-markers:	Create a file with the frequencies of each marker in each generation (default is false; be warned, marker files become large).
-h or --help:	Print this help message.
```
