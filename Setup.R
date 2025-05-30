## Set working directory to FR_metagenomes/R

## Set up installations
if (!requireNamespace("BiocManager")) {  # check if the BiocManager package is available
  install.packages("BiocManager")  # install the BiocManager package if not available
}


# Define a character vector with the required package names
packages <- c(
  "ggplot2",
  "circlize",  
  "svglite",
  "factoextra",
  "phyloseq",
  "ggpubr",
  "microbiome",
  "vegan",
  "microbiome",
  "tidyverse",
  "doParallel",
  "mia",
  "miaViz",
  "multcomp",
  "dplyr",
  "randomForest",
  "DESeq2",
  "caret",
  "jtools",
  "FactoMineR",
  "polycor",
  "magrittr",
  "remotes",
  "randomForestSRC",
  "rmarkdown",
  "mboost",
  "pdp",
  "scater",
  "NMF"
)

# Define a function to check if the packages are installed, and if not, install them
is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask = F)
  }
  sapply(pkg, require, character.only  = TRUE)
}

# Check and install the required packages
is.installed(packages)
lapply(packages, require, character.only  = TRUE)

# --------------------------------------------------------------------------

# Custom installations

BiocManager::install("ComplexHeatmap")

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)


# Upgrade dev packages 
remotes::install_github("microbiome/mia")
remotes::install_github("alanocallaghan/scater") 
remotes::install_github("KarstensLab/microshades")

# -------------------------------------------------------------------------

source("R_functions.R")  # load the functions defined in the file "R_functions.R"


