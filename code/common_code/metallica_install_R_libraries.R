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
                   "DT"
                   ))


# graphics

install.packages(c("ggplot2",
                   "plotly",
                   "patchwork",
                   "gridExtra",
                   "viridis",
                   "RColorBrewer",
                   "ggrepel",
                   "ggcorrplot",
                   "UpSetR"
                   ))

# install Bioconductor

install.packages("BiocManager")

# install Bioconductor libraries



# install devtools

install.packages("devtools")

# install libraries from github

library(devtools)

devtools::install_github("sprouffske/growthcurver")

devtools::install_github("https://github.com/vdemichev/diann-rpackage")

