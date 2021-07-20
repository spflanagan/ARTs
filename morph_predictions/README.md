README

This is a shiny app to show the frequency dependent expectations of our model of alternative reproductive tactics.


It contains the following parameters:


- freqs: A vector containing the frequencies of the other morphs (CP,CS,NP,NS). The frequencies must be labelled or all four provided in the order above. If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
- Nm: number of males in the population. Default is 500.
- Nf: number of females in the population. Default is 500.
- r: Relative reproductive effort/contribution within each clutch for parents. Default is 2/3.
- c: Sperm competition coefficient. Default is 0.5
- ws: Sexual selection strength, aka female preference for male type. Default is 1 (unidirectional preference for courters).
- wn: Selection strength on nesting trait in males, aka nest survival. Default is 1 (parental male nests survive and non-parental nests all die).
- wv: Viability selection against courtship and nesting traits. Default is exp(-0.5/(2*50)).
