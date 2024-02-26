##################################################################################################################################################
####          Script to calculate how many metabolic, signalling and TF pathways and protein complexes are affected by metal perturbation     ####
##################################################################################################################################################

#`---
#`  Title: "metperWTproteomics_metpaths_TFs_sigpaths_complexes_affected.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 11 Jan 2024
#`  Description: calculate signalling, metabolic pathways, TFs and protein complexes affected
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific paths and scripts

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/metpaths_sigpaths_TFs_complexes_affected")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/metpaths_sigpaths_TFs_complexes_affected")
dir.create(output_tables_dir,recursive = T)

signalling_pathways <- GO_gset_BP%>%
                       filter(grepl("signaling",term))%>%
                       filter(!grepl("regulation of",term))

## how many do we measure proteins from 

# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF,Significant)%>%
  unique()

lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF,Significant)
tested_ORFs <- unique(c(lmfitres_extracell$ORF), lmfitres_intracell$ORF)
significant_ORFs <- unique(c(filter(lmfitres_extracell,Significant == 1)$ORF,
                             filter(lmfitres_intracell, Significant == 1)$ORF))
###########################
### Signalling Pathways ###
###########################

signalling_pathways <- signalling_pathways%>%
                       mutate(tested = ORF %in% tested_ORFs,
                              significant = ORF%in% significant_ORFs,
                              Gene_Name = as.character(lapply(ORF, convert_ORF2GeneName)))

singalling_pathway_summary <- signalling_pathways%>%
                              group_by(term)%>%
                              summarize(total = length(ORF),
                                     num_tested = sum(tested),
                                     num_sig = sum(significant),
                                     frac_meas = num_tested/total,
                                     frac_sig = num_sig/num_tested)
write.csv(singalling_pathway_summary,
          paste0(output_tables_dir,"/signalling_pathways_measured_significant_summary.csv"),row.names = F)

#############################
### Transcription Factors ###
#############################

TFs <- GO_gset_MF%>%
       filter(grepl("transcription factor",term))%>%
        mutate(tested = ORF %in% tested_ORFs,
               significant = ORF%in% significant_ORFs)


TFs_summary <- TFs%>%
                group_by(term)%>%
                summarize(total = length(ORF),
                          num_tested = sum(tested),
                          num_sig = sum(significant),
                          frac_meas = num_tested/total,
                          frac_sig = num_sig/num_tested)


##########################
### Protein Complexes  ###
##########################

prot_complexes <- GO_gset_CC%>%
  filter(grepl("complex",term))%>%
  filter(!grepl("protein-containing complex",term))%>%
  mutate(tested = ORF %in% tested_ORFs,
         significant = ORF%in% significant_ORFs)


prot_complexes_summary <- prot_complexes%>%
  group_by(term)%>%
  summarize(total = length(ORF),
            num_tested = sum(tested),
            num_sig = sum(significant),
            frac_meas = num_tested/total,
            frac_sig = num_sig/num_tested)

write.csv(prot_complexes_summary,
          paste0(output_tables_dir,"/protein_complexes__measured_significant_summary.csv"),row.names = F)


#####################
### KEGG pathways ###
#####################

KEGG_paths <- KEGG%>%
  mutate(tested = ORF %in% tested_ORFs,
         significant = ORF%in% significant_ORFs)


KEGG_paths_summary <- KEGG_paths%>%
  group_by(term)%>%
  summarize(total = length(ORF),
            num_tested = sum(tested),
            num_sig = sum(significant),
            frac_meas = num_tested/total,
            frac_sig = num_sig/num_tested)

write.csv(KEGG_paths_summary,
          paste0(output_tables_dir,"/KEGG_measured_significant_summary.csv"),row.names = F)

