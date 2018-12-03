![PILAF-logo](logo.svg?sanitize=true)

## Overview

**PILAF** (Phylodynamics InfLuenzA Forecast) uses viral genetic data to forecast CDC influenza count time series. Viral sequence data is integrated into a dynamic model under the framework of phylodynamic inference. Integrated Nested Laplace Approximation (INLA) was used for computationally efficient inference.


## Install

```r
devtools::install_github("Mamie/PILAF")
```

## Vignette

- [simulation](vignettes/simulation.Rmd): An example showing how to simulate count time series and viral genealogy and perform forecast.

## License

PILAF is licensed under the MIT License. See LICENSE for further details.
