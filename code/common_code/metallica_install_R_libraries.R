#`---
#`  Title: "metallica_install_R_libraries.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 3 July 2023
#`  Description: Script to install all libraries used in the metallica repository
#`---

# base R
# data wrangling
install.packages(c("tidyr",
                   "dplyr",
                   "tidyverse",
                   "reshape2",
                   "readxl",
                   "DT",
                   "lubridate",
                   "doParallel",
                   "stringr"
                   ))

# bioinformatics

install.packages("gprofiler2"
                 )

# graphics

install.packages(c("ggplot2",
                   "ggpattern",
                   "plotly",
                   "patchwork",
                   "gridExtra",
                   "viridis",
                   "RColorBrewer",
                   "ggrepel",
                   "ggcorrplot",
                   "UpSetR",
                   "pheatmap",
                   "showtext",
                   "venn",
                   'R.utils'
                   ))

# install Bioconductor

install.packages("BiocManager")

# install Bioconductor libraries

BiocManager::install(pkgs = c("IRanges","S4Vectors","GenomeInfoDb","KEGGREST","org.Sc.sgd.db","piano","DEP",
                              "clusterProfiler", "PCAtools","ggkegg"))

# install devtools

install.packages("devtools")

# install libraries from github

library(devtools)

devtools::install_github("sprouffske/growthcurver")

# for circular plots in R
devtools::install_github("jokergoo/circlize")

devtools::install_github("symbioticMe/proBatch", build_vignettes = TRUE)

devtools::install_github("https://github.com/vdemichev/diann-rpackage")

install_github('ievaKer/aPEAR')
