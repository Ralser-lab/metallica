##########################################################
####          Script to plot PCA, UMAP, tSNE          ####
##########################################################

#`---
#`  Title: "metperWTproteomics_plot_PCA_UMAP_tSNE.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 14 Nov 2023
#`  Description: Script to plot PCA, UMAP and tSNE of complete matrix dataset
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

plot_dir = paste0(metpert_WTproteomics_dir,"/output/plots/completematrix/allmetals")

inp_data_dir = paste0(metpert_WTproteomics_dir,"/output/tables/completematrix/allmetals/allsamples")


prot_da_noNA <- read.csv(paste0(inp_data_dir,"/Protein Quantities wImpValues FullMatrix.csv"), stringsAsFactors = F)%>%
                dplyr::select(Protein.Ids,Log2.Protein.Quantity.,BioSpecID,Element)%>%
                unique()%>%
                group_by(BioSpecID, Element, Protein.Ids)%>%
                summarise(mean_log2_PQ = mean(Log2.Protein.Quantity.,na.rm=T))%>%
                ungroup()


library(PCAtools)
# Load the required libraries
library(ggplot2)
library(umap)


data_for_pca <- prot_da_noNA%>%
                reshape2::dcast(Protein.Ids ~ BioSpecID, value.var = "mean_log2_PQ")
rownames(data_for_pca) = data_for_pca$Protein.Ids
data_for_pca <- as.matrix(data_for_pca[,-1])

pca_res <- pca(data_for_pca,center = TRUE)
screeplot(pca_res, axisLabSize = 18, titleLabSize = 22)

biplot(p)
