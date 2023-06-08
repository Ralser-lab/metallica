
#`---
#`  Title: "prepare_KO_and_OE_metallomics_data_for_metallica_app.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 5 April 2023
#`  Description: add gene names to KO and OE metallomics data ( from Iacovacci et al 2021)
#`---

library(tidyr)
library(dplyr)
library(ggplot2)
library(gprofiler2)

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

#####################
### Read in data ###
#####################


KOmet <- read.csv(paste0(published_dataset_dir,"/metallomics/KOmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F)%>%
         mutate(KOgenename = as.character(lapply(ORF, convert_ORF2SingleGeneName)))%>%
         reshape2::melt(id.vars = c("ORF","KOgenename"),variable.name = "metal", value.name = "zscore_KOmetallomics")


OEmet <- read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F)%>%
         mutate(OEgenename = as.character(lapply(ORF, convert_ORF2SingleGeneName)))%>%
         reshape2::melt(id.vars = c("ORF","OEgenename"), variable.name = "metal", value.name = "zscore_OEmetallomics")

write.csv(KOmet, paste0(metallica_app_dir,"/metallica_app_KOmetallomics.csv"))

write.csv(OEmet, paste0(metallica_app_dir,"/metallica_app_OEmetallomics.csv"))


######################################
### summarise -- count "hits" data ###
######################################


KOmet_smry <- KOmet%>%
              filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
              mutate(significant = ifelse(abs(zscore_KOmetallomics) > 1.5, T, F))%>%
              group_by(metal)%>%
              summarise(num_significant_in_KO = sum(significant))
              


OEmet_smry <- OEmet%>%
              filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
              mutate(significant = ifelse(abs(zscore_OEmetallomics) > 1.5, T, F))%>%
              group_by(metal)%>%
              summarise(num_significant_in_OE = sum(significant))
            
pdf(paste0(published_dataset_dir,"/metallomics/KO_OE_metallomics_summary.pdf"),width = 8, height = 5 )

# KO metallomics summary
ggplot(KOmet_smry,
       aes(x = metal,
           y = num_significant_in_KO,
           fill = metal))+
  geom_bar(stat = "identity",width = 0.6)+
  theme_metallica()+
  scale_fill_manual(values = colkey_Ele)

# OE metallomics summary
ggplot(OEmet_smry,
                  aes(x = metal,
                      y = num_significant_in_OE,
                      fill = metal))+
  geom_bar(stat = "identity",width = 0.6)+
  theme_metallica()+
  scale_fill_manual(values = colkey_Ele)

dev.off()

#############################
### GSEA using gprofiler2 ###
#############################

## KO metallomics GSEA

KOm_forGSEA <- KOmet%>%
                filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
                mutate(significant = ifelse(abs(zscore_KOmetallomics) > 1.5, T, F))

# set background for enrichments 
background_ORFs_KOm <- unique(KOm_forGSEA$ORF)

# Run the gene set enrichment analysis
gsea_results_any_metal_KOm <- gost(query =  unique(filter(KOm_forGSEA,significant)$ORF),
                               organism = "scerevisiae",
                               multi_query = T, 
                               custom_bg = background_ORFs_KOm,
                               domain_scope = "custom", 
                               sources = c("GO:MF","GO:BP","KEGG","TF"))

top_terms <- gsea_results_any_metal_KOm$result %>%
             filter(p_values < 0.05)%>%
             pull(term_id)

p <- gostplot(gsea_results_any_metal_KOm, capped = TRUE, interactive = F)
go_enrichres_plot<- publish_gostplot(p, 
                                     highlight_terms = top_terms,
                                     width = 10, 
                                     height = 20,
                                     filename = paste0(plot_dir,"/KOmetallomics_gprofiler_enrich_result_anymetal",
                                                       "_overrep.pdf"))
dev.off()
## OE metallomics GSEA

OEm_forGSEA <- OEmet%>%
  filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  mutate(significant = ifelse(abs(zscore_OEmetallomics) > 1.5, T, F))

# set background for enrichments 
background_ORFs_OEm <- unique(OEm_forGSEA$ORF)

# Run the gene set enrichment analysis
gsea_results_any_metal_OEm <- gost(query =  unique(filter(OEm_forGSEA,significant)$ORF),
                                   organism = "scerevisiae",
                                   multi_query = T, 
                                   custom_bg = background_ORFs_KOm,
                                   domain_scope = "custom", 
                                   sources = c("GO:MF","GO:BP","KEGG","TF"))

top_terms <- gsea_results_any_metal_KOm$result %>%
  filter(p_values < 0.05)%>%
  pull(term_id)

p <- gostplot(gsea_results_any_metal_KOm, capped = TRUE, interactive = F)
go_enrichres_plot<- publish_gostplot(p, 
                                     highlight_terms = top_terms,
                                     width = 10, 
                                     height = 20,
                                     filename = paste0(plot_dir,"/OEmetallomics_gprofiler_enrich_result_anymetal",
                                                       "_overrep.pdf"))

dev.off()




