
#`---
#`  Title: "prepare_metpertWTproteomics_dataset_for_metallica_app.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 12 April 2023
#`  Description: reshapes all metal perturbation WT S cerevisiae proteomics data for visualisatin in app 
#`---


#############################################
### source paths functions and libraries  ###
#############################################


# general
source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

############################################################
### Get Functions and variables from EDP Analysis script ###
############################################################

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))


###################################################
### load and reshape metpert WT proteomics data ###
###################################################


metpert_WTprot_data <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
                       mutate(metal = Element,
                              metal_concentration = Element.Concentration,
                              log2_FCvsAE = Log2FC_vs_AE,
                              protein_name = Genes,
                              protein_name_ORF = paste(protein_name,ORF,sep="|")
                            )%>%
                       dplyr::select(protein_name_ORF,metal_concentration,
                                     metal,BioSpecID,log2_FCvsAE,
                                     PValAdj_BH, LeastComplexModel,Significant)%>%
                       unique()%>%
                       na.omit()


write.csv(metpert_WTprot_data, paste0(metallica_app_dir,"/metallica_app_metpertWTproteomics.csv"),row.names = F)





