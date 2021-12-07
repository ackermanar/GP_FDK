# GP_FDK

Project goals are to analyze the potential advantages of using VIBE imaging (https://www.vibeia.com/) for phenotyping of fusarium damage kernels over traditional methods.

## Data

phenotypes = FDK_ALLDATA_hmpEntries_NoLA.csv
kinship matrix developed using A.mat = LD08_numeric.kin.txt

## Scripts

FDK_MultivariateGP.R - uses each phenotype as a covariate to predict DON content, testing effectiveness of all individual phenotypes as covariates. Hypothesis is that predicting DON using VIBE as a covariate will have more accurate predictions than DON alone, and DON + other FDK phenotypes.
