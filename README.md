# VGRegion-Hospital-Wastewater-Selection-Project
On selection toward antibiotic-resistant bacteria from hospital and municipal wastewater in the Västra Götaland Region of west Sweden


To run the scripts for the statistical analysis of resistance rates in e.coli and ARG carriage rates,*ecoli_resistance_rates_models.R* and *metagenomic_arg_models.R* , the following R packages are required along with an R installation ( R version 4.4.0 tested on windows x64-based laptop). Typical installation time of the packages is around 10 minutes. 

| Package   | Version   | Purpose                                         |
|-----------|-----------|-------------------------------------------------|
| tidyverse | 2.0.0     | Data wrangling                                  |
| glmmTMB   | 1.1.14    | Negative binomial mixed effects modelling       |
| emmeans   | 2.0.1    | Estimated marginal means and contrasts          |
| DHARMa    | 0.4.7    | Checking residuals and model fit          |

For the script *ecoli_resistance_rates_models.R*, the input dataset is the *Merged Counts Ecoli* sheet of the source data excel sheet, named:
And the expected outputs are found under output/

To run the script *metagenomic_arg_models.R*, the input dataset is supplementary information excel sheet, named:
And the expected outputs are found under output/

