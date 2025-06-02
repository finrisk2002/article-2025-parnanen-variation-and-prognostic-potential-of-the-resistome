
# General

This repository contains code for the FINRISK 2002 article Variation and
prognostic potential of the gut antibiotic resistome in the FINRISK 2002
cohort, 2025.
Preprint: Population variation and prognostic potential of gut antibiotic resistome
K. Pärnänen, M. Ruuskanen, G. Sommeria-Klein, V. Laitinen, P. Kantanen, G. Méric, 
C. Gazolla Volpiano, M. Inouye, R. Knight, V. Salomaa, A. S. Havulinna, T. Niiranen,
L. Lahti
medRxiv 2024.08.08.24311663; doi: https://doi.org/10.1101/2024.08.08.24311663

Source code:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15574151.svg)](https://doi.org/10.5281/zenodo.15574151)

Code contributed by Katariina Pärnänen, Matti Ruuskanen, Ville Laitinen
and Leo Lahti.

# Instructions

Open R

```         
R
```

Set working directory

Run `Setup.R` script from the that installs all necessary packages

For MAC users, install gcc via Homebrew and follow these instructions to
set up installations requiring c++/gcc/gfortran
<https://medium.com/biosyntax/following-up-library-dependency-when-compiling-r-packages-89f191b9f227>

sessionInfo.txt contains the exact package versions.

Run `Carpentry_TreeSE.R` script that creates all the necessary data
objects. Enough to run once. Note that this requires access to sensitive
phenotype data from FINRISK.

```         
source("Setup.R") # Enough to run once. Installation code for main packages, check also sessionInfo.txt to get exact versions

source("Carpentry_TreeSE.R") # Creates the TSE object from mapping data and sample data files.
source("Metaphlan4_TSE_ES.r") # Updates the TSE object to include Metaphlan4 and nmfs 
```

Then run the individual R files for specific analyses:

"Ready-made" figures. Each figure script is self-contained and runs
indepedently.

## Figures

**Figure 1**

-\> done manually. BioRender

**Figure 2**

`source("Figure2_finland.R")` \# -\> Fig2.pdf

**Figure 3**

`source("Figure3.R")` \# -\> Fig3.jpg

**Figure 4**

`source("Figure4.R")` \# -\> Fig4.pdf

## Supplementary material

### Data

**Data S1** `source("Figure2_finland.R")` -\> runs:
source("Linear_models_ARG_load_TreeSE.R") -\> generates: DataS1.xlsx

**Data S2** `source("DataS2.R")` \# DataS2.xlsx

**Data S3** `source("Figure2_finland.R")` -\> DataS3.xlsx

## Supplementary Figures

**Supplementary Figure S1** `source("FigS1.R")`

**Supplementary Figure S2** `source("Figure2_finland.R")` -\> runs
GLM_plots_TreeSE.R -\> FigS2.pdf

**Supplementary Figure S3** `source("Figure2_finland.R")` \# -\> runs:
source("Figure2_regional_trends.R") \# -\> generates:
FigureS3_regional_trends_fitted.pdf

**Supplementary Figure S4** `source("FigS4_ARGbarplot.R")` \# -\>
FigS4_ARGbarplot.pdf

**Supplementary Figure S5** `source("Figure3.R")` -\>
FigS5A_Enterosignature_sampleprofile.pdf -\>
FigS5B_Enterosignature_genusprofile.pdf

**Supplementary Figure S6** \# These figures (AB + C) have been combined
manually `source("FigureS6AB.R")` \# -\> FigureS6AB.pdf
source("Figure3.R") \# -\> FigureS6C.pdf

**Supplementary Figure S7** `source("Figure3.R")` \# -\> FigureS7.pdf,
`source("FigureS7_meta4.R")` -\> FigureS7_meta4.jpg

**Supplementary Figure S8** - survival figures and tables \# survival
tables (all-cause, cause-specific, sepsis)
`source("Figure_supplement_survival.R")` \# -\> FigS8.pdf

**Supplementary Figure S9** - Library size and ARG metrics
`source("FigureS9.R")` -\> FigS9.pdf

**Supplementary Figure S10** `source("FigureS10.R")`

### Supplementary Tables

**Supplementary Table S1** `source("Tables_S1_S9.R"`) \# -\>
Table_S1.xlsx

**Supplementary Table S2** `source("Figure2_finland.R")` \# This calls
source("Figure2_map.R") -\> Table_S2.xlsx

**Supplementary Table S3** `source("Figure2_finland.R")` \# This calls
source("Figure2_regional_trends.R") -\> Table_S3.xlsx

**Supplementary Table S4** `source("Table_S4.R")` \# -\>
Table_S4_ARG_ALUE_pairwise.xlsx

**Supplementary Table S5** - ES-ARG associations `source("Figure3.R")`
-\> Table_S5_ES_ARG_association.xlsx

**Supplementary Table S6** `source("Figure_supplement_survival.R")` \#
-\> Extended_data_table_SURVIVAL.xlsx

**Supplementary Table S7** `source("Figure_supplement_survival.R")` \#
-\> Extended_data_table_SURVIVAL.xlsx

**Supplementary Table S8** `source("Figure_supplement_survival.R"`) \#
-\> Extended_data_table_SURVIVAL.xlsx

**Supplementary Table S9** `source("Tables_S1_S9.R")` \# -\>
Table_S9.xlsx

**Supplementary Tables S10-S11**

`source("GLMs_TreeSE.r")` \# -\> TableS10.xlxs, TableS11.xlsx

**Supplementary Table S12**

`source("Additional_checks.r")` \# -\> Table12.xlsx, DATAS1_women.xlsx.
DATAS1_women.xlxs sheets copied manually as a part of DataS1.xlsx
