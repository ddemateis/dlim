# dlim

### Overview


This repository contains the source files for the ```dlim``` R package that implements distributed lag interaction models as described in the paper:

* Demateis D, Keller KP, Rojas-Rueda D, Kioumourtzoglou M-A, Wilson A (2024). “Penalized Distributed Lag Interaction Model: Air Pollution, Birth Weight and Neighborhood Vulnerability.” _Environmetrics_. ([DOI: 10.1002/env.2843](https://doi.org/10.1002/env.2843),
[arXiv preprint](https://arxiv.org/abs/2401.02939))

Simulation scripts from the manuscript can be found at https://github.com/ddemateis/dlim_simulations.
 
The ```dlim``` package estimates a distributed lag model with modification by a single continuous factor. If you are instead interested in a DLM with modification by a single categorical or binary variable see the [```bdlim```](https://anderwilson.github.io/bdlim/) package. If you are interested in distributed lag models with heterogeneity by with multiple modifiers see the heterogeneous distributed lag model in the [```dlmtree```](https://danielmork.github.io/dlmtree/) package.


### Installation


Install from CRAN:
```
install.packages("dlim")
```

Install from GitHub:
```
remotes::install_github("ddemateis/dlim")
```

A vignette is available on the package website or in the package using the following code.
```
vignette('dlimOverview')
```


