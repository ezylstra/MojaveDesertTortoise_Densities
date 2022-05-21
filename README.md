# Modeling rangewide densities of Mojave desert tortoises
Modeling density of Mojave desert tortoises based on distance-sampling surveys, 2001-2020

________________________________________________________________________________________________________

## Code
1. [MDT_Import_Format_Data.R](MDT_Import_Format_Data.R): R code to import all the data needed for analysis (distance-sampling data, telemetry data, covariate data), format the data, and output everything into MDT_Data.Rdata.
2. [MDT_DensityModel.R](MDT_DensityModel.R): R code to load MDT_Data.Rdata, identify covariates to include in the model, run the model in Stan, and view results.
3. [MDT_DensityModel.stan](MDT_DensityModel.stan): Stan model file.
4. [MDT_DensityModel_Predictions.R](MDT_DensityModel_Predictions.R): R code to predict tortoise density throughout the study area based on covariate values and parameter estimates from the model.
