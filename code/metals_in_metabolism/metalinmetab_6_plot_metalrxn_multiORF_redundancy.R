#`---
#`  Title: "plot_metalrxn_multiORF_redundancy"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to analyze enzyme category wise changes in protein abundance + other datasets
#`---


library(tidyr)
library(RColorBrewer)
library(plotly)
library(readr)

#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))


## specific
source(paste0(code_dir,"/metpert_ecYeast8simulation/metpert_ecYeast8simulation_0_libraries_functions.R"))

plot_dir <- paste0(metpert_sim_vs_exp_comparison_dir,"/output/plots")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_sim_vs_exp_comparison_dir,"/output/tables")
dir.create(output_tables_dir,recursive = T)


#############################################################################
### Plot protein abundance of potential metal depedent redundant proteins ###
#############################################################################

## ENV CONC gradient DA

redundancy_multiORFs_env <- read.csv(paste0(metpert_sim_vs_exp_comparison_dir,
                                            "/metabolic_network_properties/redundancy_rxns_multipleORF_envDA_correlations.csv"),
                                     stringsAsFactors = F) %>%
  # create an identifier of unique ORF pairs and ORF pair - metal combinations
                            mutate(ORF_pair_ID = paste(pmin(Significant.ORF, Measured.Parallel.ORF), 
                                                                 pmax(Significant.ORF, Measured.Parallel.ORF), sep="_"),
                                   ORF_pair_ID_metal = paste(ORF_pair_ID, Element, sep = "_"))%>%
                            dplyr::select(Pearson.Correlations, Pearson.pvalues, Pearson.Direction, ORF_pair_ID_metal)%>%
                            unique()%>%
                            separate(ORF_pair_ID_metal, into = c("ORF1","ORF2","metal"), sep = "_",remove = F)

print(paste("number unique enzyme pairs and metal comb tested:",
            length(unique(redundancy_multiORFs_env$ORF_pair_ID_metal))))


print(paste("number unique enzyme pairs and metal comb significantly correlated:",
            length(unique(filter(redundancy_multiORFs_env,Pearson.Direction != "low")$ORF_pair_ID_metal))))

print(paste("number unique enzyme pairs and metal comb negatively correlated:",
            length(unique(filter(redundancy_multiORFs_env,Pearson.Direction == "negative")$ORF_pair_ID_metal))))


print(paste("number unique enzyme pairs and metal comb tested:",
            length(unique(paste(redundancy_multiORFs_env$ORF1,
                                redundancy_multiORFs_env$ORF2)))))

ORFs2keep <- unique(unlist(strsplit(c(redundancy_multiORFs_env$ORF1,
                                      redundancy_multiORFs_env$ORF2), ",")))

# Count observations for each Pearson Direction for each Element
direction_counts_unique <- redundancy_multiORFs_env %>%
                           dplyr::select(ORF_pair_ID_metal, Pearson.Direction, metal) %>%
                           group_by(metal, Pearson.Direction, ORF_pair_ID_metal) %>%
                           summarise(n = n()) %>%
                           filter(Pearson.Direction == "negative") %>%
                           summarise(total_neg = n())

# Count observations for each Pearson Direction for each Element
direction_counts_permetal_acrossmetalduplicates_Allowed <- redundancy_multiORFs_env %>%
                                                            dplyr::select(ORF_pair_ID_metal, Pearson.Direction, metal) %>%
                                                            group_by(metal, Pearson.Direction) %>%
                                                            summarise(n = n()) 

# Define custom colors
custom_colors <- c("low" = "gray", "negative" = "darkblue", "positive" = "darkred")
# Plot
pdf(paste0(plot_dir,"/redundancy_correlations_envDAdata.pdf"),width = 12,height = 7)
 ggplot(data = unique(redundancy_multiORFs_env[,c("ORF1","ORF2",
                                                  "metal","Pearson.Correlations",
                                                  "Pearson.pvalues", "Pearson.Direction")]), 
            aes(x = metal,
                y = Pearson.Correlations,
                group = metal,
                color = Pearson.Direction, size = -log10(Pearson.pvalues))) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5))+
  scale_color_manual(values = custom_colors) +
  labs(title = "multiORF pearson correlation env DA", 
       x = " ", 
       y = "Pearson Correlations", 
       size = "-log10(Pearson p-value)", 
       color = "Pearson Direction") +
   scale_size_continuous(range = c(0.5, 5))+
  theme_metallica() +
  coord_flip()
dev.off()

#################################################
### visualise ORFs with negative correlations ###
#################################################

negative_corr_sig_ORFs <- unique(filter(redundancy_multiORFs_env, Pearson.Direction == "negative")$ORF1)
  
## read in OQs along env metal cocn

FCvsAE_allproteins_env <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                  stringsAsFactors = F)%>%
                          dplyr::select(Genes, ORF,Element,Significant,Log2FC_vs_AE,Element.Concentration)%>%
                          unique()%>%
                          ## only keep ORFs we may want to plot 
                          filter(ORF %in% ORFs2keep)



dir.create(paste0(plot_dir,"/redundancy_negcorr_env"))
for(orf in 1:length(negative_corr_sig_ORFs)){
  
  corrdata2plot = unique(filter(redundancy_multiORFs_env, Pearson.Direction == "negative" &
                                  ORF1 == negative_corr_sig_ORFs[[orf]])[,c("ORF1",
                                                                                       "ORF2",
                                                                                       "metal","Pearson.Correlations",
                                                                                       "Pearson.pvalues","Pearson.Direction")])
  
  for (metal in unique(corrdata2plot$metal)) {
    
    specific_corrdata = corrdata2plot[corrdata2plot$metal == metal, ]
    orfs2plot = unique(c(negative_corr_sig_ORFs[[orf]], specific_corrdata$ORF2))
    
    df2plot = filter(FCvsAE_allproteins_env, ORF %in% orfs2plot & Element == metal) %>%
      mutate(sig_metal_anno = ifelse(ORF %in% GO_MF_all_metal_related_anno$ORF, T, F))
    
    title_parts = sapply(orfs2plot, function(orfs) {
      pearson_corr = specific_corrdata[specific_corrdata$ORF2 == orfs, "Pearson.Correlations"]
      pearson_pval = specific_corrdata[specific_corrdata$ORF2 == orfs, "Pearson.pvalues"]
      return(paste0(metal, " ", convert_ORF2GeneName(orfs), " pc: ", round(pearson_corr, 2), ", pr_pvalue: ", format(pearson_pval, scientific = TRUE)))
    })
    plot_title = paste(title_parts, collapse = " \n ")
    
    p = ggplot(df2plot, aes(x = log2(Element.Concentration),
                            y = Log2FC_vs_AE,
                            group = Genes,
                            colour = Genes,
                            shape = sig_metal_anno)) +
      geom_point(size = 3) +
      geom_smooth(method = 'loess', se = T, alpha = 0.2) +
      facet_wrap("Element") +
      scale_color_brewer(palette = "Set2") +
      scale_shape_manual(values = c(16, 17)) +
      labs(title = plot_title,
           x = "log2(env. metal concentration)",
           y = "log2(fold difference vs control",
           colour = "Gene", 
           shape = "Metal Anno") +
      theme_metallica() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            plot.title = element_text(size = 8))
    
    # Save each plot to a separate PDF
    pdf(paste0(plot_dir,"/redundancy_negcorr_env/redundancy_negativecorr_envDA_protabun_", 
               convert_ORF2GeneName(negative_corr_sig_ORFs[[orf]]), "_", metal, ".pdf"), width = 4, height =6)
    print(p)
    dev.off() 
  }
}


##############################################################
### add gene names and write correlation df for supp table ###
##############################################################

redundancy_multiORFs_env <- redundancy_multiORFs_env%>%
                            mutate(enzyme1_genename = as.character(lapply(ORF1, convert_ORF2GeneName)),
                                   enzyme2_genename = as.character(lapply(ORF2, convert_ORF2GeneName)))%>%
                            dplyr::select(enzyme1_genename,enzyme2_genename,
                                          ORF1, ORF2, metal,
                                          Pearson.Direction, Pearson.Correlations, Pearson.pvalues)%>%
                            unique()

write.csv(redundancy_multiORFs_env, paste0(output_tables_dir,"/redundancy_enzyme_pair_correlations_envDA.csv"),row.names = F)


#####################################
#####################################
#####################################
#### Cellular metal conc DA data ####
#####################################
#####################################
#####################################


## ENV CONC gradient DA

redundancy_multiORFs_cell <- read.csv(paste0(metpert_sim_vs_exp_comparison_dir,
                                            "/metabolic_network_properties/redundancy_rxns_multipleORF_cellDA_correlations.csv"),
                                     stringsAsFactors = F) %>%
  # create an identifier of unique ORF pairs and ORF pair - metal combinations
  mutate(ORF_pair_ID = paste(pmin(Significant.ORF, Measured.Parallel.ORF), 
                             pmax(Significant.ORF, Measured.Parallel.ORF), sep="_"),
         ORF_pair_ID_metal = paste(ORF_pair_ID, Element, sep = "_"))%>%
  dplyr::select(Pearson.Correlations, Pearson.pvalues, Pearson.Direction, ORF_pair_ID_metal)%>%
  unique()%>%
  separate(ORF_pair_ID_metal, into = c("ORF1","ORF2","metal"), sep = "_",remove = F)

print(paste("number unique enzyme pairs and metal comb tested - cell DA:",
            length(unique(redundancy_multiORFs_cell$ORF_pair_ID_metal))))


print(paste("number unique enzyme pairs and metal comb significantly correlated - cell DA:",
            length(unique(filter(redundancy_multiORFs_cell,Pearson.Direction != "low")$ORF_pair_ID_metal))))

print(paste("number unique enzyme pairs and metal comb negatively correlated - cell DA:",
            length(unique(filter(redundancy_multiORFs_cell,Pearson.Direction == "negative")$ORF_pair_ID_metal))))


print(paste("number unique enzyme pairs and metal comb tested - cell DA:",
            length(unique(paste(redundancy_multiORFs_cell$ORF1,
                                redundancy_multiORFs_cell$ORF2)))))

ORFs2keep <- unique(unlist(strsplit(c(redundancy_multiORFs_cell$ORF1,
                                      redundancy_multiORFs_cell$ORF2), ",")))
#############################################
### Read in cellular DA data and filter   ###
#############################################

FCvsAE_allproteins_cell <- read.csv(paste0(metpert_WTproteomics_dir,
                                           "/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                    stringsAsFactors = F) %>%
  dplyr::select(Genes, ORF,Element, median_relative_intracellular_concentration,median_log2_foldchangevsAE,Significant)%>%
  unique()%>%
  ## only keep ORFs we may want to plot 
  filter(ORF %in% ORFs2keep)

#####################################################
### visualise ORFs with negative correlations    ###
### for CELLULAR metal concentration gradient DA ###
#####################################################

dir.create(paste0(plot_dir,"/redundancy_negcorr_cell"))
for(orf in 1:length(negative_corr_sig_ORFs)){
  
  corrdata2plot = unique(filter(redundancy_multiORFs_cell, Pearson.Direction == "negative" &
                                  ORF1 == negative_corr_sig_ORFs[[orf]])[,c("ORF1",
                                                                            "ORF2",
                                                                            "metal","Pearson.Correlations",
                                                                            "Pearson.pvalues","Pearson.Direction")])
  
  for (metal in unique(corrdata2plot$metal)) {
    
    specific_corrdata = corrdata2plot[corrdata2plot$metal == metal, ]
    orfs2plot = unique(c(negative_corr_sig_ORFs[[orf]], specific_corrdata$ORF2))
    
    df2plot = filter(FCvsAE_allproteins_cell, ORF %in% orfs2plot & Element == metal) %>%
      mutate(sig_metal_anno = ifelse(ORF %in% GO_MF_all_metal_related_anno$ORF, T, F))
    
    title_parts = sapply(orfs2plot, function(orfs) {
      pearson_corr = specific_corrdata[specific_corrdata$ORF2 == orfs, "Pearson.Correlations"]
      pearson_pval = specific_corrdata[specific_corrdata$ORF2 == orfs, "Pearson.pvalues"]
      return(paste0(metal, " ", convert_ORF2GeneName(orfs), " pc: ", round(pearson_corr, 2), ", pr_pvalue: ", format(pearson_pval, scientific = TRUE)))
    })
    plot_title = paste(title_parts, collapse = " \n ")
    
    p = ggplot(df2plot, aes(x = log2(median_relative_intracellular_concentration),
                            y = median_log2_foldchangevsAE,
                            group = Genes,
                            colour = Genes,
                            shape = sig_metal_anno)) +
      geom_point(size = 3) +
      geom_smooth(method = 'loess', se = T, alpha = 0.2) +
      facet_wrap("Element") +
      scale_color_brewer(palette = "Set2") +
      scale_shape_manual(values = c(16, 17)) +
      labs(title = plot_title,
           x = "log2(cellular metal concentration)",
           y = "log2(fold difference vs control",
           colour = "Gene", 
           shape = "Metal Anno") +
      theme_metallica() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            plot.title = element_text(size = 8))
    
    # Save each plot to a separate PDF
    pdf(paste0(plot_dir,"/redundancy_negcorr_cell/redundancy_negativecorr_cellDA_protabun_", 
               convert_ORF2GeneName(negative_corr_sig_ORFs[[orf]]), "_", metal, ".pdf"), width = 4, height =6)
    print(p)
    dev.off() 
  }
}

##############################################################
### add gene names and write correlation df for supp table ###
##############################################################

redundancy_multiORFs_cell <- redundancy_multiORFs_cell%>%
  mutate(enzyme1_genename = as.character(lapply(ORF1, convert_ORF2GeneName)),
         enzyme2_genename = as.character(lapply(ORF2, convert_ORF2GeneName)))%>%
  dplyr::select(enzyme1_genename,enzyme2_genename,
                ORF1, ORF2, metal,
                Pearson.Direction, Pearson.Correlations, Pearson.pvalues)%>%
  unique()

write.csv(redundancy_multiORFs_cell, paste0(output_tables_dir,"/redundancy_enzyme_pair_correlations_cellDA.csv"),row.names = F)

