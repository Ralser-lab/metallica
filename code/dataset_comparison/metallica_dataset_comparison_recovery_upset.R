#`---
#`  Title: "metallica_dataset_comparison.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 25 May 2023
#`  Description: scompares results of metpertWTproteomics, metdepKOgrowth, KOmetallomics, OEmetallomics and y5kmetalspecific
#`---

#############################################
### source paths functions and libraries  ###
#############################################

library(tidyr)
library(dplyr)

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))

plot_dir = paste0( datasetcomparison_dir, "/output/plots")
dir.create(plot_dir,recursive = T)
dir.create(paste0(plot_dir,"/HyperGSA"))
output_tables_dir = paste0( datasetcomparison_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)

##########################
### Input all datasets ###
##########################

## metal perturbation proteomics results  

## lms with extracellular metal concentration

metpertWTproteomics_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                 stringsAsFactors = F)[,c("Element","Genes","ORF","Significant")]
                                 

# make a list of all  proteins that were quantified in metpertWTproteomics
metpertWTproteomics_measured <- unique(metpertWTproteomics_extracell$ORF)

## lms with intracellular metal concentration

metpertWTproteomics_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                          stringsAsFactors = F)[,c("Element","Genes","ORF","Significant")]

## join the two metpertWTproteomics datasets


metpertWTproteomics <- rbind(metpertWTproteomics_extracell,metpertWTproteomics_intracell)%>%
                       data.frame()%>%
                       group_by(Element,Genes,ORF)%>%
                       mutate(Significant = ifelse(any(Significant), T, F),
                              dataset = "metpertWTproteomics",
                              dataset_type = "proteomics",
                              metal = Element)%>%
                       ungroup()%>%
                       unique()%>%
                       dplyr::select(ORF, metal, Significant, dataset, dataset_type)%>%
                       data.frame()
    
                  

###############################################################################
### Input metdepKOgrowth, KOmettalomics, OEmetallomice and y5kmetalspecific ###
###############################################################################

y5kmetalspecific <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),stringsAsFactors = F)%>%
                                     mutate(Significant = ifelse(p_value < 0.05 & abs(log2FC) > log2(1.5),T,F),
                                              metal = term)%>%
                                     group_by(metal,Protein.Group)%>%
                                     mutate(Significant = ifelse(any(Significant),T,F))%>%
                                     ungroup()%>%
                                     dplyr::select(Protein.Group,metal,Significant)%>%
                                     unique()%>%
                                     na.omit()%>%
                                     mutate(ORF = as.character(lapply(Protein.Group, convert_Uniprot2singleORF)),
                                            dataset = "y5kmetalspecific",
                                            dataset_type = "proteomics")%>%
                                     dplyr::select(ORF, metal, Significant, dataset, dataset_type)

y5kmetalspecific_measured <- unique(y5kmetalspecific$ORF)

## phenotypic ( growth ) screen on metal depletion conditions of the BY4741 + pHLUM KO library

metdepKOgrowth <- read.csv(paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),stringsAsFactors = F)%>%
                  mutate(Significant = sig_in_metal,
                         dataset = "metdepKOgrowth",
                         dataset_type = "KOgrowth")%>%
                  dplyr::select(ORF, metal, Significant, dataset, dataset_type)%>%
                  unique()
metdepKOgrowth_measured <- unique(metdepKOgrowth$ORF)
                 
## metallomics screen on metal depletion conditions of the BY4741 KO library

KOmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/KOmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
                 reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_KOmetallomics")%>%
                 filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
                 mutate(Significant = ifelse(abs(Zscore_KOmetallomics) > 1.959, T,F),
                        dataset = "KOmetallomics",
                        dataset_type = "metallomics")%>%
                 dplyr::select(ORF, metal, Significant, dataset, dataset_type)%>%
                 unique()


KOmetallomics_measured <- unique(KOmetallomics$ORF)
                 
## metallomics screen on metal depletion conditions of the OE library ( check the exact name ? whats the BY number)

OEmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
                 reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_OEmetallomics")%>%
                 filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
                 mutate(Significant = ifelse(abs(Zscore_OEmetallomics) > 1.959, T,F),
                        dataset = "OEmetallomics",
                        dataset_type = "metallomics")%>%
                 dplyr::select(ORF, metal, Significant, dataset, dataset_type)%>%
                 unique()

OEmetallomics_measured <- unique(OEmetallomics$ORF)

###################################################
### Comparing what was measured across datasets ###
###################################################

## all measured 

measured_across_ds <- list(metpertWTproteomics_measured,
                           y5kmetalspecific_measured,
                           metdepKOgrowth_measured,
                           KOmetallomics_measured,
                           OEmetallomics_measured)

names(measured_across_ds) <- c( "metpertWTproteomics",
                                "y5kmetalspecific",
                                "metdepKOgrowth",
                                "KOmetallomics",
                                "OEmetallomics")

pdf(paste0(plot_dir,"/venn_measured_acrooss_datasets.pdf"),width = 10, height = 10)
venn::venn(measured_across_ds,zcolor = "style",box = F, ilcs=2.5) 
dev.off()

###################################################
### Make venn of number shared between datasets ###
###################################################

# Get unique ORFs
all_ORFs <- unique(unlist(measured_across_ds))

# Create a data frame
ORF_df <- data.frame(ORF = all_ORFs, 
                     metpertWTproteomics = sapply(all_ORFs, function(x) x %in% measured_across_ds$metpertWTproteomics),
                     y5kmetalspecific = sapply(all_ORFs, function(x) x %in% measured_across_ds$y5kmetalspecific),
                     metdepKOgrowth = sapply(all_ORFs, function(x) x %in% measured_across_ds$metdepKOgrowth),
                     KOmetallomics = sapply(all_ORFs, function(x) x %in% measured_across_ds$KOmetallomics),
                     OEmetallomics = sapply(all_ORFs, function(x) x %in% measured_across_ds$OEmetallomics))

# Add a new column to the data frame that contains the sum of TRUE values
ORF_df$total <- rowSums(ORF_df[2:6])

# Create the summary data frame
summary_df <- as.data.frame(table(ORF_df$total))
colnames(summary_df) <- c("number_of_ds_measured_in", "number_ORFs")
 
pdf(paste0(plot_dir,"/num_datasets_measuredin_vs_num_ORFs.pdf"),width = 7,height = 4)
ggplot(summary_df,
       aes(x = number_of_ds_measured_in,
           y = number_ORFs))+
  geom_bar(stat = "identity", width = 0.5,
           colour = "white",
           fill = "#3A5FCDB3")+
  geom_text(aes(y = number_ORFs + 100,
                label = number_ORFs),
            size = 6)+
  theme_metallica()+
  labs(x = "number of datasets",
       y = "number of ORFs")
dev.off()

###########################################################
### Plot number of hits across all datasets in one plot ###
###########################################################

all_datasets <- rbind(metdepKOgrowth,KOmetallomics,OEmetallomics,y5kmetalspecific,metpertWTproteomics)%>%
                filter(metal != "AE")%>%
                data.frame()%>%
                mutate(dataset = factor(dataset,
                                       levels = c("metdepKOgrowth","KOmetallomics","OEmetallomics",
                                                  "y5kmetalspecific","metpertWTproteomics")))

write.csv(all_datasets, paste0(output_tables_dir,"/all_datasets_combined_significant_or_not.csv"),row.names = F)


################################################################################
### Count what % of uniprot annotation score > 1 proteins were sig in any ds ###
################################################################################

length(unique(all_datasets$ORF))

all_datasets_sig_in_any <- all_datasets %>%
  dplyr::select(ORF, Significant, dataset)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(sig_in_any = any(Significant))%>%
  ungroup()

sum(unique(all_datasets_sig_in_any[,c("ORF","sig_in_any")])$sig_in_any)

all_datasets_sig_in_any_unp_1more <- merge(all_datasets_sig_in_any, GenProt_SGD, by = "ORF")%>%
                           filter(Uniprot.Annotation.Score > 1)
  
sum(unique(all_datasets_sig_in_any_unp_1more[,c("ORF","sig_in_any")])$sig_in_any)

############################################################################
### count number of ORFs that are shared hits in 1, 2,3,4 and 5 datasets ###
############################################################################

# Calculate the number of unique dataset types for each ORF
ORF_counts <- filter(all_datasets, Significant == T) %>%
  group_by(ORF) %>%
  summarise(n_distinct_dataset = n_distinct(dataset)) %>%
  ungroup()

write.csv(ORF_counts, paste0(output_tables_dir,"/ORFs_num_dataset_significant_in.csv"),row.names = F)

ORF_counts<-ORF_counts%>%
  group_by(n_distinct_dataset) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  arrange(desc(n_distinct_dataset)) %>%
  mutate(cumulative_count = cumsum(count))


print("number of ORFs that were significant in x number of datasets")
ORF_counts
write.csv(ORF_counts,paste0(output_tables_dir,"/number_datasets_sig_in_x_datasets.csv"),row.names = F)
## metal wise

# Count the unique dataset types for each ORF_metal combination
ORF_metal_counts <- filter(all_datasets, Significant == T) %>%
  mutate(ORF_metal = paste(ORF, metal, sep = "_"))%>%
  group_by(ORF_metal,metal) %>%
  summarise(n_distinct_dataset= n_distinct(dataset))

write.csv(ORF_metal_counts, paste0(output_tables_dir,"/ORF_metal_numdataset_significant_in.csv"),row.names = F)

## number hits in > 1 dataset by metal

ORF_metal_counts_numhits_morethan1_ds_per_metal <- ORF_metal_counts%>%
                                                   filter(n_distinct_dataset >1)%>%
                                                   group_by(metal)%>%
                                                   summarise(num_ORFs_sig_in_morethan1ds = n())%>%
                                                   ungroup()
write.csv(ORF_metal_counts_numhits_morethan1_ds_per_metal,paste0(output_tables_dir,"/number_ORFs_sig_morethan2ds_permetal.csv"),row.names = F)

ORF_metal_counts<- ORF_metal_counts %>%
  group_by(n_distinct_dataset) %>%
  summarise(count = n())%>%
  ungroup()%>%
  arrange(desc(n_distinct_dataset)) %>%
  mutate(cumulative_count = cumsum(count))

print("number of ORF-metal combinations that were significant in x number of datasets")
ORF_metal_counts
write.csv(ORF_metal_counts,paste0(output_tables_dir,"/number_datasetmetal_pairs_sig_in_x_datasets.csv"),row.names = F)


###############################################
### Plot number of hits across all datasets ###
###############################################

num_hits_all_datasets <- all_datasets%>%
                         group_by(metal,dataset, dataset_type)%>%
                         summarise(num_significant = sum(Significant, na.rm = T))%>%
                         ungroup()%>%
                         data.frame()

pdf(paste0(plot_dir,"/num_hits_across_all_datasets.pdf"),width = 8,height = 5)
ggplot(filter(num_hits_all_datasets, num_significant > 0),
       aes(x = metal,
           y = num_significant,
           colour = dataset_type,
           fill = dataset_type)) +
 # geom_bar(stat = "identity", position = "dodge", width = 0.5, colour = "white")+
  geom_point(aes(shape = dataset),
             size = 4) +
  scale_y_log10() +
  scale_colour_manual(values = colkey_dataset_type) +
  scale_fill_manual(values = colkey_dataset_type) +
  scale_shape_manual(values = c(19,2,24, 22, 8)) +
  labs(x = "", y = "number signficant",colour = "", fill = "",shape = "")+
  theme_metallica()
dev.off()

##############################################
### Upset Plot of hits across all datasets ###
##############################################

all_datasets_hits <-  list(unique(filter(all_datasets,dataset == "metdepKOgrowth" & Significant)$ORF),
                           unique(filter(all_datasets,dataset == "KOmetallomics" & Significant)$ORF),
                           unique(filter(all_datasets,dataset == "OEmetallomics" & Significant)$ORF),
                           unique(filter(all_datasets,dataset == "y5kmetalspecific" & Significant)$ORF),
                           unique(filter(all_datasets,dataset == "metpertWTproteomics" & Significant)$ORF)
                           )
names(all_datasets_hits) <- c( 
                    "metdepKOgrowth",
                    "KOmetallomics",
                    "OEmetallomics",
                    "y5kmetalspecific",
                    "metpertWTproteomics")

pdf(paste0(plot_dir,"/upset_all_significant_metalORF_across_datasets.pdf"),width = 9, height = 6)

upset(fromList(all_datasets_hits), 
      keep.order = T,
      order.by = "freq",
      cutoff = 100,
      nintersects = NA,
      main.bar.color ="#CD919E"  ,
      sets.bar.color = "#BCD2EE",
      text.scale = 2,
      point.size = 3,
      matrix.color = "#2b0b57",
      mb.ratio = c(0.5, 0.5)

)
dev.off()

#############################################
### metal wise upset plots for validation ###
#############################################

metals = unique(all_datasets$metal)

pdf(paste0(plot_dir,"/upset_all_significant_across_datasets_metalwise.pdf"),width = 8, height = 4)

for(m in 1:length(metals)){
  
  all_datasets_ele = filter( all_datasets, metal == metals[m])
  
  all_datasets_ele_hits <-  list(unique(filter(all_datasets_ele,dataset == "metdepKOgrowth" & Significant)$ORF),
                             unique(filter(all_datasets_ele,dataset == "KOmetallomics" & Significant)$ORF),
                             unique(filter(all_datasets_ele,dataset == "OEmetallomics" & Significant)$ORF),
                             unique(filter(all_datasets_ele,dataset == "y5kmetalspecific" & Significant)$ORF),
                             unique(filter(all_datasets_ele,dataset == "metpertWTproteomics" & Significant)$ORF)
  )
  names(all_datasets_ele_hits) <- c( 
    "metdepKOgrowth",
    "KOmetallomics",
    "OEmetallomics",
    "y5kmetalspecific",
    "metpertWTproteomics")
  
print(
  upset(fromList(all_datasets_ele_hits), nsets = 6, order.by = "freq",
        cutoff = 6,
        matrix.color = "black", 
        main.bar.color = colkey_Ele[metals[m]],
        sets.bar.color =  colkey_Ele[metals[m]])
)
  
}
dev.off()


######################################################################
### cdf plots all dataset hits ### --  keeping hits of all studies ###
######################################################################

get_plot_recovery_summary <- function(all_datasets, specific_annotation_df, unspecific_annotation_df, type_of_anno){
  
  # note what % of specific_annotation_df were measured in at least one screen 
  meas_in_any_ds <- unique(all_datasets$ORF)
  print(paste("total ORF intesect",type_of_anno,round(
    length(intersect(unique(specific_annotation_df$ORF), meas_in_any_ds)) / length(unique(specific_annotation_df$ORF)),2)))
  
  recov_in_any_ds <- unique(filter(all_datasets, Significant)$ORF)
  print(paste("total ORF intesect",type_of_anno,round(
    length(intersect(unique(specific_annotation_df$ORF), recov_in_any_ds)) / length(unique(specific_annotation_df$ORF)),2)))
  
  meas_in_any_ds_ORFmetcombo <- unique(paste(all_datasets$ORF, all_datasets$metal))
  spec_anno_ORFmetcombo <- paste(specific_annotation_df$ORF, specific_annotation_df$term)
  
  print(paste("ORF-met combo intersect",type_of_anno,round(length(intersect(meas_in_any_ds_ORFmetcombo,spec_anno_ORFmetcombo))/
                                               length(unique(spec_anno_ORFmetcombo)),2)))
  
  print(paste("unspecific ORF intersect",type_of_anno,round(length(intersect(meas_in_any_ds,unspecific_annotation_df$ORF))/
                length(unique(unspecific_annotation_df$ORF)),2)))
  
  all_datasets_merge_anno <- merge(filter(specific_annotation_df,
                                          ORF %in% meas_in_any_ds), all_datasets, 
                                   by.x = c("ORF","term") ,
                                   by.y = c("ORF","metal"))
  
  recovery_anno <- all_datasets_merge_anno %>%
    group_by(ORF, term) %>%
    mutate(recovered = ifelse(is.na(Significant), F, Significant),
           recovered_by_any = any(recovered),
           recovered_by_any = ifelse(is.na(recovered_by_any), F, recovered_by_any))
  
  # all dataset recovery
  recovery_smry <- recovery_anno %>%
    dplyr::select(ORF, term, recovered_by_any) %>%
    unique() %>%
    group_by(term) %>%
    summarize(total_annotype = length(ORF),
              num_recovered = sum(recovered_by_any),
              frac_recovered = num_recovered/total_annotype)
  
  # ds wise recovery 
  recovery_smry_dswise <- recovery_anno %>%
    mutate( dataset = factor(dataset,
                             levels = c("metdepKOgrowth","KOmetallomics","OEmetallomics",
                                        "y5kmetalspecific","metpertWTproteomics"))) %>%
    group_by(term, dataset) %>%
    summarize(total_annotype = length(ORF),
              num_recovered = sum(recovered),
              frac_recovered = num_recovered/total_annotype) %>%
    ungroup()
  
  pdf(paste0(plot_dir,"/", type_of_anno,"_recovery_all_datasets.pdf"),width = 10,height = 6)
  print(
    ggplot(recovery_smry_dswise,
         aes(x = dataset))+
    geom_bar(aes(y = total_annotype,
                 colour = term),
             stat = "identity",
             position = "dodge", width = 0.5,
             alpha = 0)+
    geom_bar(aes(y = num_recovered,
                 colour = term,
                 fill = term),
             stat = "identity",
             position = "dodge", width = 0.5,
             alpha = 0.8)+
    scale_fill_manual(values = colkey_Ele)+
    scale_colour_manual(values = colkey_Ele)+
    theme_metallica()+
    labs(x = "",
         y = paste0("fraction of ", type_of_anno, " recovered"))
  )
  dev.off()
  
  
  ## count how many known metal binders are recovered cumulatively by each dataset for each metal 
  
  colnames(recovery_anno)[which(colnames(recovery_anno)== "term")] <- "metal"
  recovery_anno$ORF_metal = paste(recovery_anno$ORF, recovery_anno$metal)
  
  cuml_rec_metdepKOgrowth <- filter(recovery_anno,dataset == "metdepKOgrowth" & recovered)%>%
    group_by(metal)%>%
    mutate(num_cuml_recovered = length(ORF))
  
  cuml_rec_KOmetallomics <- filter(recovery_anno,dataset == "KOmetallomics" & recovered &
                                     !ORF_metal %in% unique(cuml_rec_metdepKOgrowth$ORF_metal))%>%
    group_by(metal)%>%
    mutate(num_cuml_recovered = length(ORF))
  
  
  cuml_rec_OEmetallomcis <- filter(recovery_anno,dataset == "OEmetallomics" & recovered & 
                                     !ORF_metal %in% c(unique(cuml_rec_metdepKOgrowth$ORF_metal,
                                                              cuml_rec_KOmetallomics$ORF_metal)))%>%
    group_by(metal)%>%
    mutate(num_cuml_recovered = length(ORF))
  
  
  cuml_rec_y5kmetalspecific <- filter(recovery_anno,dataset == "y5kmetalspecific" & recovered &
                                        !ORF_metal %in% unique(c( cuml_rec_metdepKOgrowth$ORF_metal,
                                                                  cuml_rec_KOmetallomics$ORF_metal,
                                                                  cuml_rec_OEmetallomcis$ORF_metal)))%>%
    group_by(metal)%>%
    mutate(num_cuml_recovered = length(ORF))%>%
    ungroup()
  
  cuml_rec_metpertWTproteomics <- filter(recovery_anno,dataset == "metpertWTproteomics" & recovered &
                                           !ORF_metal %in% unique(c(cuml_rec_metdepKOgrowth$ORF_metal,
                                                                    cuml_rec_KOmetallomics$ORF_metal,
                                                                    cuml_rec_OEmetallomcis$ORF_metal,
                                                                    cuml_rec_y5kmetalspecific$ORF_metal)))%>%
    group_by(metal)%>%
    mutate(num_cuml_recovered = length(ORF))%>%
    ungroup()
  
  for_ecdf <- merge(rbind(  unique(cuml_rec_metdepKOgrowth[,c("metal","dataset","num_cuml_recovered")]),
                            unique(cuml_rec_KOmetallomics[,c("metal","dataset","num_cuml_recovered")]),
                            unique(cuml_rec_OEmetallomcis[,c("metal","dataset","num_cuml_recovered")]),
                            unique(cuml_rec_metpertWTproteomics[,c("metal","dataset","num_cuml_recovered")]),
                            unique(cuml_rec_y5kmetalspecific[,c("metal","dataset","num_cuml_recovered")])
                            
  ), recovery_smry[,c("term","total_annotype")], by.x = "metal", by.y = "term")%>%
    mutate(frac_cuml_recovered = num_cuml_recovered/total_annotype)%>%
    mutate(dataset = gsub("Significant_","", dataset),
           dataset = factor(dataset, levels = c("metdepKOgrowth","KOmetallomics","OEmetallomics",
                                                "y5kmetalspecific","metpertWTproteomics")))%>%
    group_by(metal)%>%
    arrange(dataset)%>%
    mutate(cumsum = cumsum(frac_cuml_recovered))
  
  pdf(paste0(plot_dir,"/",type_of_anno,"_cumfracrecovery_all_datasets.pdf"),width = 7,height = 7)
  print(ggplot(for_ecdf,
         aes(x = dataset,
             y = cumsum,
             colour = factor(metal),
             group = metal))+
    geom_point( alpha = 0.9, size = 3)+
    geom_line(linewidth = 1.5, alpha = 0.8)+
    scale_colour_manual(values = colkey_Ele)+
    theme_metallica()+
    theme(axis.text.x = element_text(angle = 90))+
    labs(x = "",
         y = paste0("cumulative % recovery of",type_of_anno),
         colour = "")
  )
   dev.off()
  
  
  ## Plot number of unspecific metal binders or metal transporters that were recovered in each 
  
  # filter dataframe for ORFs that have an unspecific annotation only - i.e. no metal annotation at all
  
  unspecific_anno_df <- unspecific_annotation_df%>%
    filter(!ORF %in% unique(specific_annotation_df$ORF))
  
  recovery_unspecific  <- merge(filter(unspecific_anno_df,
                                       ORF %in% meas_in_any_ds), all_datasets, 
                                by= c("ORF"))%>%
    group_by(metal,dataset)%>%
    dplyr::summarise(total_measured = length(ORF),
                     num_recovered = sum(Significant, na.rm = T))%>%
    ungroup()%>%
    filter(num_recovered >0)%>%
    mutate(xlab = dataset,
           xlab = ifelse(grepl("metdepKOgrowth",xlab),"kog",
                         ifelse(grepl("KOmetallomics",xlab),"kom",
                                ifelse(grepl("OEmetallomics",xlab),"oem",
                                       ifelse(grepl("metpertWTproteomics",xlab),"wtp",
                                              ifelse(grepl("y5kmetalspecific",xlab),"y5k",xlab))))),
           xlab = paste(xlab,metal))
  
  unsp_ORF_meas = unique(merge(filter(unspecific_anno_df,
                               ORF %in% meas_in_any_ds), all_datasets, 
                        by= c("ORF"))$ORF)
  
  unsp_ORF_sig = unique(merge(filter(unspecific_anno_df,
                                      ORF %in% meas_in_any_ds), all_datasets, 
                               by= c("ORF")))%>%
                filter(Significant)%>%
                pull(ORF)%>%
                unique()
  
  unsp_ORF_sig_1metalonly = unique(merge(filter(unspecific_anno_df,
                                                ORF %in% meas_in_any_ds), all_datasets, 
                                         by= c("ORF")))%>%
                             filter(Significant)%>%
                             dplyr::select(ORF,metal)%>%
                             unique()%>%
                             group_by(ORF)%>%
                             summarize(num_metals_sig_in = n())%>%
                             ungroup()%>%
                             filter(num_metals_sig_in <2)
  
  print(paste("unspecific", type_of_anno, "measured across datasets:",length(unsp_ORF_meas)))
  print(paste("unspecific", type_of_anno, "signficant across datasets:",length(unsp_ORF_sig)))
  
  print(paste("unspecific", type_of_anno,"recovery total ORFs:",round(length(unsp_ORF_sig)/length(unsp_ORF_meas),2)))
  print(paste("unspecific", type_of_anno, "only significant in 1 metal:",nrow(unsp_ORF_sig_1metalonly)))
  write.csv(unsp_ORF_sig_1metalonly,paste0(output_tables_dir,"unspecific",type_of_anno,"sig_only_for_1metal.csv"),row.names = F)
  
  pdf(paste0(plot_dir,"/unspecific_", type_of_anno,"_recovery_all_datasets.pdf"),width = 10,height = 8.15)
  print(ggplot(recovery_unspecific, aes(x = xlab, y = num_recovered, color = metal)) +
    geom_segment(aes(xend = xlab, yend = 0, colour = metal), linewidth = 0.2) +
    geom_point(position = position_dodge(0.9), size = 3) +
    scale_color_manual(values = colkey_Ele) +
    theme_metallica() +
    labs(x = "", y = paste0("Number unspecific", type_of_anno, " recovered"))+
    theme(axis.text.x = element_text(angle = 90))
  )
  dev.off()
  
}


## calculate and plot recovery of metal binding proteins

get_plot_recovery_summary(all_datasets, 
                          specific_annotation_df = GO_gset_MF_binding_metalwise,
                          unspecific_annotation_df = GO_gset_MF_metalbinding_unspecific,
                          type_of_anno = "metalbinding")

## calculate and plot recovery of metal specific metal transporters # only GO annotations

#get_plot_recovery_summary(all_datasets, 
 #                         GO_gset_MF_transporter_metalwise, 
  #                        GO_gset_MF_metaltransporter_unspecific, "metaltransporter_GO")

## calculate and plot recovery of metal specific metal transporters #  GO and philpott review annotations

#get_plot_recovery_summary(all_datasets, rbind(GO_gset_MF_transporter_metalwise,philpott_metal_transporter_df),
 #                         GO_gset_MF_metaltransporter_unspecific,"metaltransporter_GO_philpott")

############################################################################
### plot % of metal binders and % metal transporters sig in all datasets ###
############################################################################

## all metals together 


all_ds_metal_binding_allmet <- merge(all_datasets,rbind(GO_gset_MF_binding_metalwise,
                                                       GO_gset_MF_metalbinding_unspecific
                                                        ),
                                    by.x = c("ORF"),
                                    by.y = c("ORF"))%>%
  group_by(ORF,dataset)%>%
  mutate(Significant = ifelse(any(Significant),T,F))%>%
  ungroup()%>%
  dplyr::select(dataset,Significant,ORF)%>%
  unique()%>%
  group_by(dataset)%>%
  mutate(total = length(ORF))%>%
  ungroup()%>%
  group_by(dataset, total)%>%
  reframe(num_sig = sum(Significant,na.rm=T),
          perc_sig = round(100*(num_sig / total),2))%>%
  ungroup()%>%
  unique()%>%
  mutate(type = "metal binder")



all_ds_metal_transport_allmet <- merge(all_datasets,rbind(GO_gset_MF_transporter_metalwise,
                                                          philpott_metal_transporter_df),
                                     by.x = c("ORF"),
                                     by.y = c("ORF"))%>%
  group_by(ORF,dataset)%>%
  mutate(Significant = ifelse(any(Significant),T,F))%>%
  ungroup()%>%
  dplyr::select(dataset,Significant,ORF)%>%
  unique()%>%
  group_by(dataset)%>%
  mutate(total = length(ORF))%>%
  ungroup()%>%
  group_by(dataset, total)%>%
  reframe(num_sig = sum(Significant,na.rm=T),
          perc_sig = round(100*(num_sig / total),2))%>%
  ungroup()%>%
  unique()%>%
  mutate(type = "metal transporter")


all_ds_metal_other_allmet <- merge(all_datasets,rbind(GO_gset_MG_other_metalrelated),
                                       by.x = c("ORF"),
                                       by.y = c("ORF"))%>%
  group_by(ORF,dataset)%>%
  mutate(Significant = ifelse(any(Significant),T,F))%>%
  ungroup()%>%
  dplyr::select(dataset,Significant,ORF)%>%
  unique()%>%
  group_by(dataset)%>%
  mutate(total = length(ORF))%>%
  ungroup()%>%
  group_by(dataset, total)%>%
  reframe(num_sig = sum(Significant,na.rm=T),
          perc_sig = round(100*(num_sig / total),2))%>%
  ungroup()%>%
  unique()%>%
  mutate(type = "other metal related")


all_ds_bindtransp <- rbind(all_ds_metal_binding_allmet, all_ds_metal_transport_allmet, all_ds_metal_other_allmet)

pdf(paste0(plot_dir,"/metal_binder_metal_transporter_recovery_by_dataset.pdf"),width = 7, height = 5)
ggplot(all_ds_bindtransp,
       aes(x = dataset,
           colour = type,
           fill = type))+
  geom_bar(aes(y = total), fill = NA, stat = "identity", position = "dodge2", width = 0.5)+
  geom_bar(aes(y = num_sig), stat = "identity", position = "dodge2", width = 0.5)+
  geom_text(aes(y = num_sig+30, label = paste0(round(perc_sig,0),"%")), size = 4, color = "black")+
  theme_metallica()+
  facet_wrap("type",scales = "free_x")+
  scale_fill_brewer(palette = "Pastel1")+
  scale_color_brewer(palette = "Pastel1")+
  coord_flip()+
  labs(y = "number of ORFs",
       x = "",
       colour = "",
       fill = "")+
  theme(legend.position = "none")

dev.off()




all_ds_metal_binding_metsp <- merge(all_datasets,GO_gset_MF_binding_metalwise,
                                    by.x = c("ORF","metal"),
                                    by.y = c("ORF","term"))%>%
                              group_by(metal,dataset)%>%
                              mutate(total = length(ORF),
                                     type = "metal binders")%>%
                              ungroup()%>%
                              group_by(metal,dataset, total)%>%
                              reframe(num_sig = sum(Significant,na.rm=T),
                                        perc_sig = round(100*(num_sig / total),2))%>%
                              ungroup()%>%
                              unique()
                             
ggplot(all_ds_metal_binding_metsp,
       aes(x = dataset,
           color = metal,
           fill = metal))+
  geom_bar(aes(y = perc_sig),  width = 0.6, stat = "identity", position = "dodge2")+
  geom_text(aes(y = perc_sig +0.5, label = paste0(num_sig,"/",total)),size = 3)+
  scale_color_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  theme(axis.text.x= element_text(angle = 90,hjust = 0.9))+
  labs(x = "",
       y = "% specific metal binders significant")



## specific metal transporters 

all_ds_metal_transp_sp <- merge(all_datasets,rbind(GO_gset_MF_transporter_metalwise,
                                                          philpott_metal_transporter_df),
                                    by.x = c("ORF","metal"),
                                    by.y = c("ORF","term"))%>%
  group_by(metal,dataset)%>%
  mutate(total = length(ORF),
         type = "metal transporters")%>%
  ungroup()%>%
  group_by(metal,dataset, total)%>%
  reframe(num_sig = sum(Significant,na.rm=T),
          perc_sig = round(100*(num_sig / total),2))%>%
  ungroup()%>%
  unique()


ggplot(all_ds_metal_transp_sp,
       aes(x = dataset,
           color = metal,
           fill = metal))+
  geom_bar(aes(y = perc_sig),  width = 0.6, stat = "identity", position = "dodge2")+
  geom_text(aes(y = perc_sig +0.5, label = paste0(num_sig,"/",total)),size = 3)+
  scale_color_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  theme(axis.text.x= element_text(angle = 90,hjust = 0.9))+
  labs(x = "",
       y = "% specific metal transporters significant")

write.csv(rbind(all_ds_metal_binding_metsp,all_ds_metal_transp_sp),
          paste0(output_tables_dir,"/all_dataset_recovert_metal_binders_transporters.csv"), row.names = F)

#################################################################
### Run elewise, datasetwise HyperGSA on hits in each dataset ###
#################################################################

## function for running hyperGSA per element per dataset

run_metalwisedbwise_hyperGSA <- function(hits_df,ds_name){
  
  HyperGSA_res_df = data.frame()
  
  metals = unique(hits_df$metal)
  
  for(m in 1:length(metals)){
    
    metal_df_to_hypGSA <- na.omit(filter(hits_df, metal == metals[m])[,c("ORF","Significant")])
    
    if(sum(metal_df_to_hypGSA$Significant) > 5){
      
    hres <- run_HyperGSA(data.frame(metal_df_to_hypGSA),
                              EnrichPV_thresh=0.05)
    
    if(length(hres)>0){
      
      hres$dataset = ds_name
      hres$metal = metals[m]
      
      HyperGSA_res_df <- rbind(HyperGSA_res_df,hres)
    }
    }
  }
  return(HyperGSA_res_df)
}
# filter alldatasets dataframe for each datasetthen run hyperGSA on each dataframe
hyperGSA_elewise_res_all_datasets <- vector()

ds_names <- as.character(unique(all_datasets$dataset))

for(ds in 1:length(ds_names)){
  
  
  for_hypgsa <- filter(all_datasets, dataset == ds_names[ds])[,c("metal","ORF","Significant")]
  
  hyperGSA_res <- run_metalwisedbwise_hyperGSA(for_hypgsa,ds_names[ds])
  hyperGSA_elewise_res_all_datasets <- rbind(hyperGSA_elewise_res_all_datasets, hyperGSA_res)
  
}

write.csv(hyperGSA_elewise_res_all_datasets, paste0(output_tables_dir,"/HyperGSAresults_alldatasets_metalwise.csv"),row.names = F)

###########################################################################
### Make 3 column Sankey plot -- dataset to metal to gset term enriched ###
###########################################################################

enriched_gsnames <- unique(hyperGSA_elewise_res_all_datasets$Gset.Type)

for(egs in 1:length(enriched_gsnames)){
  sn <- plot_threecolumn_dataset2metal2genesetterm_Sankey(hyperGSA_elewise_res_all_datasets,enriched_gsnames[egs])
  htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/HyperGSA/HyperGSA_dswisemetalwise_Sankey",enriched_gsnames[egs],".html"))
  plotly::export(sn, file = paste0(plot_dir,"/HyperGSA/HyperGSA_dswisemetalwise_Sankey",enriched_gsnames[egs],".pdf"))
}


###########################################################################







