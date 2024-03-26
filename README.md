#  A spatially-explicit model for density that accounts for availability: a case study with Mojave desert tortoises

### Erin R. Zylstra, Linda J. Allison, Roy C. Averill-Murray, Vincent Landau, Nathaniel S. Pope, and Robert J. Steidl

### Ecosphere [10.1002/ecs2.4448](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.4448)

### Code/Data DOI: [![DOI](https://doi.org/10.5281/zenodo.7402616.svg)](https://doi.org/10.5281/zenodo.7402616)
_______________________________________________________________________________________________________________________________________

## Abstract:
Estimating population density and identifying those areas where density is changing through time are central to prioritizing conservation and management strategies. Obtaining reliable estimates of density and trends can be challenging, however, especially for long-lived species that are rare, have broad geographic distributions, and are difficult to detect reliably during field surveys. We developed a hierarchical model for distance-sampling data that characterizes spatial variation in density at two scales and simultaneously estimates regional trends while accounting for variation in detection probability and availability across surveys. We applied the model to data collected over a 20-year period (2001–2020) in an area that encompassed most of the geographic range of the Mojave desert tortoise (Gopherus agassizii). Density of adult tortoises varied with multiple biotic and abiotic features, including topography, aspect, geology, and seasonal precipitation and temperature regimes. Across the entire period and study area, the density of adult tortoises decreased by an average of 1.8% per year (95% CI = −3.5% to −0.2%). Trends varied geographically, however, with the steepest declines in the western part of the range (−4.1%, −6.9% to −1.3%). Accounting for habitat loss across our study area, the abundance of this threatened species declined by an estimated 129,000 adults (36%) between 2001 and 2020. Our modeling approach extends traditional distance-sampling frameworks by accounting for ecological and observational processes that could mask spatiotemporal variation in density and, at the same time, provides spatially explicit estimates to guide conservation and management strategies for tortoises and other rare species.
_______________________________________________________________________________________________________________________________________

## Code
1. [MDT_Import_Format_Data.R](MDT_Import_Format_Data.R): R script to import all the data needed for analysis (distance-sampling data, telemetry data, covariate data), format the data, and output everything into MDT_Data.Rdata.
2. [MDT_DensityModel.R](MDT_DensityModel.R): R script to load MDT_Data.Rdata, identify covariates to include in the model, and run the model in Stan (calling MDT_DensityModel.stan). This script is also used to summarize results and produce a figure with log-linear trends in each recovery unit.
3. [MDT_DensityModel.stan](MDT_DensityModel.stan): Stan model file.
4. [MDT_DensityModel_Predictions.R](MDT_DensityModel_Predictions.R): R script to predict tortoise density throughout the study area based on covariate values and parameter estimates from the model.  This script is also used to create a figure depicting the study area (transect and plot locations, recovery units, and tortoise conservation areas).

## Data
Survey data for Mojave desert tortoises are not provided here, as the species is listed as threatened under the U.S. Endangered Species Act. 

## Output
1. [PredictedDensity_2020.tif](PredictedDensity_2020.tif): Raster with predicted density in 2020 of adult Mojave desert tortoises in 1 sq km grid cells throughout much of their geographic range. Excludes the Upper Virgin River Recovery Unit, cells with >40% impervious surfaces, and cells with values of environmental covariates that were >10% outside the range of values at survey locations. 
