# Estimating density of Mojave desert tortoises based on distance-sampling surveys, 2001-2020

## Code
1. [MDT_Import_Format_Data.R](MDT_Import_Format_Data.R): R script to import all the data needed for analysis (distance-sampling data, telemetry data, covariate data), format the data, and output everything into MDT_Data.Rdata.
2. [MDT_DensityModel.R](MDT_DensityModel.R): R script to load MDT_Data.Rdata, identify covariates to include in the model, and run the model in Stan (calling MDT_DensityModel.stan). This script also used to summarize results, and produce figures (study area map, trends)
3. [MDT_DensityModel.stan](MDT_DensityModel.stan): Stan model file.
4. [MDT_DensityModel_Predictions.R](MDT_DensityModel_Predictions.R): R script to predict tortoise density throughout the study area based on covariate values and parameter estimates from the model.
