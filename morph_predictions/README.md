# README

This is a shiny app to show the frequency dependent expectations of our model of alternative reproductive tactics.


## The model

The model pertains to a population with four male reproductive morphs (courter-parents, courter-sneakers, noncourter-parents, and noncourter-sneakers). We are predicting which paramter combinations will allow multiple morphs to be maintained in the population. 

Here is the model:

<img src="frequency_dependence.png" alt="model" width="500"/>

It contains the following parameters:

**Adjustable:**

- initCP: Initial frequency of the courter-parent morph in the population. You can control this in the slider bar and see its effect when initNP==0 in the left-hand plot.
- initNP: Initial frequency of the noncourter-parent morph in the population. You can control this in the slider bar and see its effect when initCP==0 in the right-hand plot.
- initCN: Initial frequency of the courter-sneaker morph in the population. These values are plotted in the two graphs on the y-axis.
- initNN: Initial frequency of the noncourter-sneaker morph in the population. These values are plotted in the two graphs on the x-axis.
- c: Sperm competition coefficient. You can adjust this in the slider bar.
- r: Relative reproductive effort/contribution within each clutch for courting morphs, bounded between 0 and 1. This is the proportion of a female's eggs (aka female fecundity) that a given morph can fertilize, given his maximum amount of sperm. Default simulation conditions are 8/12 (0.6667) for non-courting morphs (noncourter-parent and noncourter-nonparent) and 4/12 (0.3333) for courting morphs (courter-parent and courter-nonparents). You can adjust this in the slider bar. 


**Not currently adjustable:**

- Nm: number of males in the population. Set to 500.
- Nf: number of females in the population. Set to 500.
- ws: Sexual selection strength, aka female preference for male type. Set to 1 (unidirectional preference for courters).
- wn: Selection strength on nesting trait in males, aka nest survival. Set to 1 (parental male nests survive and non-parental nests all die).
- wv: Viability selection against courtship and nesting traits. Set to exp(-0.5/(2*50)).


The output is an estimate of the diversity (Shannon's diversity, estimated using `diversity` in the package `vegan`) of the morphs in the population after 1000 generations. These estimates are plotteda s the values on the contour plots.




