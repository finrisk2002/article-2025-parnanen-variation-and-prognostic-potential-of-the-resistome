#install.packages("ape")
#install.packages("ggtree")
#install.packages("tidyverse")

# I separated data loading from tree construction, so we can more
# swiftly run the consctruction in scripts where data is already
# available.

source("loadData.R")
source("doTree.R")

