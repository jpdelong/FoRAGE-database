# FoRAGE-database
Code used to process functional response data contained in FoRAGE

This repository contains scripts and functions used to process functional response data contained in the FoRAGE database, which can be found on the Knowledge Network for Biodiversity site:

https://knb.ecoinformatics.org/view/doi:10.5063/F17H1GTQ

Files ending in .m are run in Matlab.
The main pipeline is ‘FR_fitting_pipeline_{date}.m’, which is run as a script. It calls the function ‘model_dataset.m’ to generate simulated full datasets from data reported in a paper as means with errors. The pipeline is set up to run in parallel on a local machine, so it also calls ‘par_save.m’ to save each dataset and parameter output as it cycles through all of the functional response curve data in the overall database.
