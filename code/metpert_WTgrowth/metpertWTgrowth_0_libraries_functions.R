#`---
#`  Title: "metpertWTgrowth_0_libraries_functions.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 8 April 2021 
#`  Description: Script to analyze all ICP-MS data from intracellular metallomics measurements
#`---
options(repos=structure(c(CRAN="http://cran.r-project.org")))

#install.packages('/camp/lab/ralserm/working/Simran Aulakh/metallica/packages/growthcurver_0.3.0.tar.gz', repos=NULL, type='source')

require(tidyverse)
require(readxl)
require(lubridate)
require(DT)
require(ggplot2)
require(plotly)
require(patchwork)
require(growthcurver)
require(reshape2)
require(gridExtra)
require(pheatmap)
require(viridis)

