# Best variable subset for modeling

This project implements a genetic algorithm to discover a good (if not best) subset of variables to use.

## Dependencies

* foreach

## Suggested

* doMC -- strongly recommended, for running the algorithm on multiple cores with greater speed
* devtools -- optional, for building the package
* testthat -- optional, for running tests to verify functionality

## Usage:

```
# OPTIONAL -- for parallelism
library(doMC)
registerDoMC(cores=4) # Where 4 corresponds to the number of cores in your machine to utilize

library(bestSubsetGA)
my_subset <- bestSubset()
dependent_variables <- my_data[, my_subset] # Take only the given columns
```

The `bestSubset()` function returns a vector of the designated clusters for each row in dataset.  See `help(bestSubset)` for more information or `tests/` for more examples.


## Installation:

```
# If you do not already have the "devtools" package installed
install.packages("devtools")

devtools::install_github("rdtaylor/bestSubsetGA")
```

Alternately, use the `build.sh` script in the root of the project.
