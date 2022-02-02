
# General information 
This folder includes all the necessary code to 
-estimate WAHE indicator based on cross-sectional data, including health status defined across three single dimensions of health: 
1.chronic diseases (2 levels), 2.limitations in activity of daily living (GALI,3 levels), 3.self-rated health (5 levels), 
and simultaneously across these 3 dimensions (multiple, 2x3x5=30 levels)
-replicate figures and tables in the article "Well-being Adjusted Health Expectancy-a New Summary Measure of Population Health"

#Data
The original application of the WAHE indicator was based on the The EU Statistics on Income and Living Conditions (EU-SILC), 2018 cross-sectional data files. This data is only available for registered and approved users, more information on the process of obtaining access to the data: https://ec.europa.eu/eurostat/web/microdata

#Files included in the order to be run
"dataprep.r"- data preparation

"weights.r" - estimate well-being weights from coefficients of ordered probit regression models: dependent variable: well-being, indepndent variable: state of health, all models are estimated separately by country and sex and in all models we control for age and age-squared

"prevalence.r" <- estimate prevalence of health states by country, sex and age-group

"SMPH.r" - estimate HE and WAHE at age 15 based on the weights and prevalence files, also LE(15), read-in data for DALE (15), write out a file with all the SMPH by country and sex


"plotweights.r", "plotprevalence.r", "tables_figures.r" - replicate all the figures included in the article and supplementary material
