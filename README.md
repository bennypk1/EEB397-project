# EEB397 Project

## Project Overview

This is the repository containing all the code and datasets used to complete a year-long independent research course: EEB397: Research Project in Ecology and Evolutionary Biology at the University of Toronto. The project explores crowding sensitivity (i.e. spatially-constrained density dependent mextinction mechanisms) in a resource species, and it's potential effects on higher-order consumers. To complete the project, pair-approximated ODEs were specified, largely based off of the work of Liao et al. (2017), Liao et al. (2013), and Pillai et al. (2010). Each ODE repredents a specific community motif (e.g. Omnivory ; see `ODEFinal.jl`). Once specified, models were graphically explored (using the functions in `plot_ODEFinal.jl`) and behaviour / results deemed interesting were systematically tested against diverse parameter combinations to verify their robustness. The projectg was primarily carried out in Julia, however some R code was used for parameter set generation and some post-hoc data analysis. All datasets are stored as CSV files.

## File Navigation

Scripts in the `LEGACY` directory are no longer in active use, but contain code used in previous versions of the
project. Scripts in `CONSTANTS` contain constant objects (such as lists of hex colors) useful for model visualization, for providing standardized model parameters, etc. Scripts in `SENSITIVITY ANALYSIS` are for running a sensitivity analysis accross multiple parameter combinations, as well as generating the parameter sets in the first place. Code scripts as well as datasets have self-explanatory names, and comments throughout describing their internal functionality.

## Licence
