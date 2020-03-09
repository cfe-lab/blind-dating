# Run this with 'Rscript config.R' to make sure the correct packages are installed

install.packages(c(
"BiocManager",
"optparse",
"ape",
"ggplot2",
"dplyr",
"magrittr",
"tidytree",
"lubridate",
"chemCal"
))

BiocManager::install("ggtree")
