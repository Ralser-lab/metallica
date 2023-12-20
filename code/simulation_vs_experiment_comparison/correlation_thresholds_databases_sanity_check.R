### Correlation sanity check ###
#`---
#`  Title: "analyze_enzyme_changes_per_metal_v2.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to analyze enzyme category wise changes in protein abundance + other datasets
#`---


library(tidyr)
library(RColorBrewer)
library(plotly)
library(corrplot)

#########################
### Script specific functions ###
#################################


transform_anno_df <- function(anno_df) {
  if(ncol(anno_df)==2){
    colnames(anno_df)<- c("ORF","group")
  } else if (ncol(anno_df)==3){
    colnames(anno_df)<- c("ORF","group_db1","group_db2")
    
    anno_df <- anno_df %>%
      mutate(group = paste(group_db1, group_db2, sep = "_")) %>%
      dplyr::select(ORF, group) %>%
      unique()
  } else if (ncol(anno_df)==4){
    colnames(anno_df)<- c("ORF","group_db1","group_db2", "group_db3")
    
    anno_df <- anno_df %>%
      mutate(group = paste(group_db1, group_db2, group_db3, sep = "_")) %>%
      dplyr::select(ORF, group) %>%
      unique()
  } else if (ncol(anno_df)==5){
    colnames(anno_df)<- c("ORF","group_db1","group_db2", "group_db3","group_db4")
    
    anno_df <- anno_df %>%
      mutate(group = paste(group_db1, group_db2, group_db3, group_db4, sep = "_")) %>%
      dplyr::select(ORF, group) %>%
      unique()
  } else{
    stop("too many columns in annotation dataframe. Check it !!")
  }
  
  return(anno_df)
}

calculate_groupwise_corr <- function(PQ_df, metconc_type, anno_df, anno_name){
  
  # 
  PQ_df_annotated <- merge(PQ_df, anno_df, by = "ORF")
  
  groups_to_keep <- PQ_df_annotated%>%
    dplyr::select(ORF,group)%>%
    unique()%>%
    group_by(group)%>%
    summarise(num_ORF = n())%>%
    ungroup()%>%
    # filter out anything wiht < 5 observations
    filter(num_ORF > num_ORF_min & num_ORF < num_ORF_max )%>%
    pull(group)
  
  
  ## run correlations on protein quantity vs metal concentration in each group
  
  if(metconc_type =="env"){
    PQ_df_annotated <- filter(PQ_df_annotated, group %in% groups_to_keep)%>%
      drop_na(Log2FC_vs_AE, Element.Concentration) %>%
      group_by(group, Element) %>%
      summarise(
        correlation = cor.test(Log2FC_vs_AE, log2(Element.Concentration), method = "spearman")$estimate,
        p_value = cor.test(Log2FC_vs_AE, log2(Element.Concentration), method = "spearman")$p.value,
        num_ORF = length(unique(ORF)),
      )%>%
      ungroup()%>%
      mutate(p_value_adj_global = p.adjust(p_value),
             database = anno_name,
             metconc_type = metconc_type)%>%
      group_by(Element)%>%
      mutate(p_value_adj_metalwise = p.adjust(p_value))%>%
      ungroup()
  }else{
    PQ_df_annotated <- filter(PQ_df_annotated, group %in% groups_to_keep)%>%
      drop_na(median_log2_foldchangevsAE, median_relative_intracellular_concentration) %>%
      group_by(group, Element) %>%
      summarise(
        correlation = cor.test(median_log2_foldchangevsAE, log2(median_relative_intracellular_concentration), method = "spearman")$estimate,
        p_value = cor.test(median_log2_foldchangevsAE, log2(median_relative_intracellular_concentration), method = "spearman")$p.value,
        num_ORF = length(unique(ORF))
      )%>%
      ungroup()%>%
      mutate(p_value_adj_global = p.adjust(p_value, method = "BH"),
             database = anno_name,
             metconc_type = metconc_type)%>%
      group_by(Element)%>%
      mutate(p_value_adj_metalwise = p.adjust(p_value, method = "BH"))%>%
      ungroup()
  }
  
  return(PQ_df_annotated)
  
}
#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/enzyme_changes")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/enzyme_changes")
dir.create(output_tables_dir,recursive = T)


##################
## Input Files ###
##################


## Environmental metal conc

FCvsAE_allproteins_env<- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                  stringsAsFactors = F)%>%
  dplyr::select(Genes, ORF,Element,Significant,Log2FC_vs_AE,Element.Concentration)%>%
  unique()


FCvsAE_allproteins_cell<- read.csv(paste0(metpert_WTproteomics_dir,
                                          "/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                   stringsAsFactors = F)%>%
  dplyr::select(Genes, ORF,Element, median_relative_intracellular_concentration,median_log2_foldchangevsAE,Significant)%>%
  unique()


####################################################
### Overlap between metabolism related databases ###
####################################################

metabolism_related_dbs <- list(
  ## genome scale metabolic model
  unique(Sc_GEM$ORF),
  # Enzyme Classificatio numbers
  unique(EC_df$ORF),
  # metal binding or transport
  unique(c(GO_gset_MF_binding_metalwise$ORF, GO_gset_MF_metalbinding_unspecific$ORF,
           GO_gset_MF_transporter_metalwise$ORF, GO_gset_MF_metaltransporter_unspecific$ORF)),
  # KEGG
  unique(KEGG$ORF))

names(metabolism_related_dbs) <- c( "ScGEM",
                                    "EC",
                                    "metal_b+t",
                                    "KEGG")

PQ_dfs <- list(FCvsAE_allproteins_env, FCvsAE_allproteins_cell)
names(PQ_dfs) <- c("env","cell")


anno_dfs <- list(KEGG,
                 EC_df[,c("ORF","EC_level1_name")],
                 EC_df[,c("ORF","EC_level1_name","EC_level_2")],
                 GOslim_CC,
                 
                 merge(KEGG,EC_df[,c("ORF","EC_level1_name")], by = "ORF"),
                 merge(KEGG,EC_df[,c("ORF","EC_level1_name","EC_level_2")], by = "ORF"),
                 merge(KEGG, GOslim_CC, by = "ORF"),
                 
                 merge(EC_df[,c("ORF","EC_level1_name")],GOslim_CC, by = "ORF"),
                 merge(EC_df[,c("ORF","EC_level1_name","EC_level_2")],GOslim_CC, by = "ORF"),
                 
                 merge(merge(KEGG, GOslim_CC, by = "ORF"),EC_df[,c("ORF","EC_level1_name")], by = "ORF"),
                 merge(merge(KEGG, GOslim_CC, by = "ORF"),EC_df[,c("ORF","EC_level1_name","EC_level_2")], by = "ORF"))

# Apply the transformation function to each data frame in the list
anno_dfs <- lapply(anno_dfs, transform_anno_df)

names(anno_dfs) <- c("KEGG",
                     "EC",
                     "EC_level_1_2",
                     "CC",
                     
                     "KEGG-EC",
                     "KEGG-EC_level_1_2",
                     "KEGG-CC",
                     
                     "EC-CC",
                     "EC_level_1_2-CC",
                     
                     "KEGG-CC-EC",
                     "KEGG-CC-EC_level_1_2")

corr_check_results <- vector()

for(nmin in seq(3,21, by = 3)){
  
  num_ORF_min <- nmin
  
  for(nmax in seq(80,500,by = 20)){
  
    num_ORF_max <- nmax
    
      all_corr_df <- vector()
      
      for(adb in 1:length(anno_dfs)){
        
        for(pqd in 1:length(PQ_dfs)){
          
          all_corr_df <- rbind(all_corr_df, 
                               calculate_groupwise_corr(PQ_df = PQ_dfs[[pqd]],
                                                        metconc_type= names(PQ_dfs)[pqd],
                                                        anno_df = anno_dfs[[adb]],
                                                        anno_name = names(anno_dfs)[adb])
          )
          
        }
      }
      ### Summarize and visualise all correlation results
      
      
      all_corr_df <- all_corr_df%>%
        mutate(high_sig_spearmancorr = ifelse(p_value_adj_global < 0.05 & abs(correlation) > 0.3,T,F),
               database = factor(database, levels = c(
                 "KEGG-CC-EC_level_1_2",
                 "KEGG-EC_level_1_2",
                 
                 "KEGG-EC",
                 "KEGG-CC-EC",
                 "KEGG-CC",
                 "EC_level_1_2-CC",
                 "EC_level_1_2",
                 
                 "KEGG",
                 "EC-CC",
                 "EC",
                 "CC")),
               metconc_type = factor(metconc_type, levels = c("env","cell")))%>%
        filter(num_ORF > num_ORF_min & num_ORF < num_ORF_max )
      
      correlation_check_df <- unique(all_corr_df[,c("correlation","num_ORF","database")])%>%
        group_by(database)%>%
        summarise(correlation_spearman = cor.test(abs(correlation),num_ORF)$estimate,
                  correlation_pvalue = cor.test(abs(correlation),num_ORF)$p.value
        )%>%
        ungroup()
      
   
      num_terms_per_db <- all_corr_df%>%
                          dplyr::select(database,group)%>%
                          unique()%>%
                          group_by(database)%>%
                          summarise(num_groups = length(group))%>%
                          ungroup()%>%
                          summarise(median_num_groups = median(num_groups,na.rm= T))
      
      correlation_check_summary <- correlation_check_df %>%
                                   mutate(usable_db = correlation_pvalue >= 0.05 & abs(correlation_spearman) < 0.2)%>%
                                   summarise(num_db_low_corr = sum(usable_db),
                                             perc_usable = num_db_low_corr/length(database))
      
      ggplot(correlation_check_df,
             aes(x = correlation_spearman,
                 y = -log10(correlation_pvalue),
                 color = database))+
        geom_point(size = 5)+
        geom_text_repel(aes(label=database))+
        labs(color = "")+
        geom_hline(yintercept = -log10(0.05), color = "red")+
        theme_metallica()
      
     # pdf(paste0(plot_dir,"/spearman_correlations_min9_max160.pdf"),width =16,height = 13)
      ggplot(unique(all_corr_df[,c("correlation","num_ORF","database")]),
             aes( x = abs(correlation),
                  y = log2(num_ORF)))+
        geom_point()+
        theme_metallica()+
        facet_wrap("database",scales ="free")+
        geom_smooth(method = "lm", colour = "red")
      #dev.off()
      
      
      print("thresholds for usable databases:")
      print("spearmann pvalue >= 0.05 and abs(correlation_spearman) < 0.3")
      print("correlation calculated between number of ORFs in a group and correlation of the group with any metal concentration")
      print(paste("at num_ORF threshold:","min-",num_ORF_min,"max-",num_ORF_max,"only",correlation_check_summary$num_db_low_corr,"or",
                  100*round(correlation_check_summary$perc_usable,2),"% dbs are usable"))
      
      correlation_check_summary$min = num_ORF_min
      correlation_check_summary$max = num_ORF_max
    
      corr_check_results <- rbind(corr_check_results,correlation_check_summary)
  }
}   


corr_check_results
write.csv(corr_check_results,paste0(output_tables_dir,"/num_ORF_thresholds_for_corr_pearson_sanitycheck.csv",row.names = F))


####################################################
### Visualise the correlation check optimisation ###
####################################################

opt_min = min(filter(corr_check_results,perc_usable == max(corr_check_results$perc_usable))$min)
opt_max = max(filter(corr_check_results,perc_usable == max(corr_check_results$perc_usable))$max)


pdf(paste0(plot_dir,"/spearman_correlation_sanitycheck.pdf"),width = 8,height = 6)
ggplot(corr_check_results,
       aes(x = min,
           y = max,
           colour = perc_usable,
           fill = perc_usable))+
  geom_tile()+
  scale_color_viridis()+  
  scale_fill_viridis()+
  theme_metallica()+
  geom_vline(xintercept = opt_min, colour = "red")+
  geom_hline(yintercept = opt_max, colour = "red")+
  labs(title = paste("optimum thresholds are min",opt_min,"max",opt_max),
       fill = "% database usable",
       color = "% database usable")+
  theme(title = element_text(size = 12))

dev.off()

