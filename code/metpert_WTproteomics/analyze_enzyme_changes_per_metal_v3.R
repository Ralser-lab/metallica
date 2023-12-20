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
    filter(num_ORF > 9 & num_ORF < 160 )%>%
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

# Create a function to merge with metalwise_binding_transport_db and calculate percentages
calculate_metal_dependence_per_group_db <- function(anno_df,metal_anno_db) {
  
  # Merge with metalwise_binding_transport_db using "ORF" column
  merged_df <- merge(anno_df, metal_anno_db, by = "ORF", all.x = T)
  
  result <- merged_df %>%
    group_by(ORF) %>%
    mutate(metal_anno_present = any(!is.na(term)))%>%
    dplyr::select(-term)%>%
    unique()%>%
    group_by(group)%>%
    summarise(total = length(ORF),
              num_metal_related = sum(metal_anno_present),
              frac_metal_related = num_metal_related/total)
  
  return(result)
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

####################################################
### Overlap between metabolism related databases ###
####################################################

metal_related_ORFs_df  <- unique(rbind(GO_gset_MF_binding_metalwise, GO_gset_MF_metalbinding_unspecific,
                                       GO_gset_MF_transporter_metalwise, GO_gset_MF_metaltransporter_unspecific))

metabolism_related_dbs <- list(
  ## genome scale metabolic model
  unique(Sc_GEM$ORF),
  # Enzyme Classificatio numbers
  unique(EC_df$ORF),
  # metal binding or transport
  #unique(c(GO_gset_MF_binding_metalwise$ORF, GO_gset_MF_metalbinding_unspecific$ORF,
  #        GO_gset_MF_transporter_metalwise$ORF, GO_gset_MF_metaltransporter_unspecific$ORF)),
  # KEGG
  unique(KEGG$ORF))

names(metabolism_related_dbs) <- c( "ScGEM",
                                    "EC",
                                    #   "metal_b+t",
                                    "KEGG")

pdf(paste0(plot_dir,"/venn_ScGEM_EC_KEGG_overlap.pdf"), width = 10,height = 10)
venn::venn(metabolism_related_dbs,zcolor = "style",box = F, ilcs=2.5) 
dev.off()


pdf(paste0(plot_dir,"/metabolism_related_datasets_overlap.pdf"),width = 9, height = 6)

upset(fromList(metabolism_related_dbs), 
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

########################################################################################
### Fraction of each metabolism related database that has a metal binding annotation ###
########################################################################################

## scGEM

ScGEM_metalfrac <- merge(Sc_GEM , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  dplyr::select(-term.x)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "ScGEM")

## KEGG

KEGG_metalfrac <- merge(KEGG , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  dplyr::select(-term.x)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "KEGG")

KEGGnum_metalfrac <- merge(KEGG , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(term.x)%>%
  mutate(metal_anno = ifelse(any(!is.na(term.y)),1,0))%>%
  ungroup()%>%
  dplyr::select(metal_anno,term.x)%>%
  unique()%>%
  summarise(total = length(term.x),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "KEGGnum")


## EC_df

EC_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name")] , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  dplyr::select(-EC_level1_name)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term),1,0))%>%
  ungroup()%>%
  dplyr::select(-term)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "EC")


ECnum_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name","EC_level_2")] , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(EC_level1_name,EC_level_2)%>%
  mutate(metal_anno = ifelse(any(!is.na(term)),1,0))%>%
  ungroup()%>%
  dplyr::select(metal_anno,EC_level1_name,EC_level_2)%>%
  unique()%>%
  summarise(total = length(EC_level1_name),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "ECnum_level_1_2")



metabolism_metal_frac <- rbind(ScGEM_metalfrac, KEGG_metalfrac,EC_metalfrac)

pdf(paste0(plot_dir,"/percentage_metdatabaseORFs_metalrelated.pdf"),width = 7,height = 6)
ggplot(metabolism_metal_frac,
       aes(x = reorder(database,-frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total, color = database), width = 0.5,fill = NA, linewidth = 0.5)+
  geom_bar(stat = "identity",
           aes(y = num_metanno, color = database, fill = database), width = 0.5,)+
  geom_text(aes(y = num_metanno + 50, label = paste0(round(frac_metanno,2)*100,"%")), size =4 )+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "", y = "ORFs")
dev.off()

metabolism_metal_frac_grouped <- rbind(KEGGnum_metalfrac,ECnum_metalfrac)

pdf(paste0(plot_dir,"/percentage_metdatabaseORFs_metalrelated_pathwaywise.pdf"),width = 4,height = 6)

ggplot(metabolism_metal_frac_grouped,
       aes(x = reorder(database,-frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total, color = database), width = 0.5,fill = NA, linewidth = 0.5)+
  geom_bar(stat = "identity",
           aes(y = num_metanno, color = database, fill = database), width = 0.5,)+
  geom_text(aes(y = num_metanno + 2, label = paste0(round(frac_metanno,2)*100,"%")),size = 4)+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "")+
  theme(legend.position = "none")
dev.off()


EC_level1_2_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name","EC_level_2")] , metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(ORF,EC_level1_name, EC_level_2)%>%
  mutate(metal_anno = ifelse(!is.na(term),1,0))%>%
  ungroup()%>%
  dplyr::select(-term)%>%
  unique()%>%
  group_by(EC_level1_name,EC_level_2)%>%
  summarise(total_ORFs = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total_ORFs)%>%
  ungroup()


pdf(paste0(plot_dir,"/EC_level_1_2_metaldependence.pdf"),width = 6,height = 11)

ggplot(EC_level1_2_metalfrac,
       aes(x = reorder(paste(EC_level1_name,EC_level_2),frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total_ORFs), width = 0.5,fill = NA, linewidth = 0.5, color = brewer.pal(3, name = "Set2")[[1]])+
  geom_bar(stat = "identity",
           aes(y = num_metanno), width = 0.5, fill = brewer.pal(3, name = "Set2")[[1]])+
  geom_text(size = 3,aes(y = total_ORFs + 15, label = paste0(round(frac_metanno,2)*100,"%")))+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "")+
  theme(axis.text.y= element_text(size = 8))+
  coord_flip()
dev.off()



KEGG_metalfrac <- merge(KEGG, metal_related_ORFs_df, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(ORF,term.x)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  group_by(term.x)%>%
  summarise(total_ORFs = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total_ORFs)%>%
  filter(total_ORFs > 4 & total_ORFs < 200)


pdf(paste0(plot_dir,"/KEGG_pathway_metaldependence.pdf"),width = 10,height = 20)

ggplot(KEGG_metalfrac,
       aes(x = reorder(term.x,frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total_ORFs), width = 0.5,fill = NA, linewidth = 0.5, color = brewer.pal(3, name = "Set2")[[2]])+
  geom_bar(stat = "identity",
           aes(y = num_metanno), width = 0.5, fill = brewer.pal(3, name = "Set2")[[2]])+
  geom_text(aes(y = total_ORFs + 5, label = paste0(round(frac_metanno,2)*100,"%")))+
  theme_metallica()+
  labs(fill = "", color = "")+
  coord_flip()+
  theme(axis.text.y= element_text(size = 10))
dev.off()

#####################################################################################################
### What % of the metal binding and non-metal binding proteins respond in flux, DAenv and DAint ? ###
#####################################################################################################


########################
### Input fluxes ###
#######################


GEMsim_0pt9grwth_fluxes <- read.csv(paste0(metpert_ecYeast8simulation_dir,"/output/tables/cofactorYeast8_flux_change_summary.csv"),stringsAsFactors = F)%>%
  
  ## filter out metal transport fluxes -- these were the ones that were minimized or maximized to achieve simulation 
  filter( !Reaction.Name %in% c("Ca(2+) exchange","Ca(2+) transport", ## Ca  
                                "Uptake of Cu(+)","Cu2(+) transport", "Cu2(+) exchange",
                                "iron (II) transport","	iron(3+) exchange","iron(2+) exchange",
                                "Mg(2+) exchange",
                                "Mn(2+) exchange",
                                "potassium exchange","potassium transport",
                                "sodium exchange",
                                "Zn(2+) exchange"))%>%
  ## filter out very high flux fold changes 
  filter(abs(log2(fold_change_flux)) < log2 (5) )


GEMsim_0pt9grwth_fluxes_separated<- GEMsim_0pt9grwth_fluxes %>%
  separate_rows(grRules, sep = " and | or ")


# Counting the number of NA values in the 'flux' column
na_count_flux <- sum(GEMsim_0pt9grwth_fluxes_separated$flux ==0)

# Counting the number of NA values in the 'fold_change_flux' column
na_count_fold_change_flux <- sum(is.na(GEMsim_0pt9grwth_fluxes_separated$fold_change_flux))
# Printing the counts
print(paste("Number of NA values in 'flux' column: ", na_count_flux))
print(paste("Number of NA values in 'fold_change_flux' column: ", na_count_fold_change_flux))



# Summarize fold_change_flux by the median per ORF and metal_perturbation
GEMsim_0pt9grwth_fluxes_separated_summary <- GEMsim_0pt9grwth_fluxes_separated%>%
  group_by(grRules, metal_perturbation) %>%
  summarize(median_fc_flux = median(fold_change_flux, na.rm = TRUE),
            frac_flux_change_vs_cntrl = (median(flux,na.rm=T)-median(control_condition_flux,na.rm=T))/median(control_condition_flux,na.rm=T) ) %>%
  ungroup()

colnames(GEMsim_0pt9grwth_fluxes_separated_summary)[which(colnames(GEMsim_0pt9grwth_fluxes_separated_summary)=="grRules")] <- "ORF"

### What % of metal related and unknown are changing in GEM ? ###

ScGEM_flux_metrel_changes <- merge(GEMsim_0pt9grwth_fluxes_separated_summary, 
                          metal_related_ORFs_df, by = "ORF",all.x =T)%>%
  mutate(metal_related = ifelse(!is.na(term),"metal related","unknown"))%>%
  group_by(ORF)%>%
  mutate(changing_in_any = ifelse(any(abs(log2(median_fc_flux)) > log2(1.5)),T,F ))%>%
  dplyr::select(ORF, metal_related, changing_in_any)%>%
  unique()%>%
  group_by(metal_related)%>%
  summarize(
    total_ORFs = length(ORF),
    num_ORF_changing_flux = sum(changing_in_any),
    frac_ORF_pathway_changingflux = num_ORF_changing_flux/total_ORFs)%>%
  ungroup()%>%
  filter(total_ORFs > 4)

### merge with KEGG ###

KEGG_flux_metrel_changes <- merge(merge(GEMsim_0pt9grwth_fluxes_separated_summary,
                                        KEGG, by = "ORF"), 
                                  metal_related_ORFs_df, by = "ORF",all.x =T)%>%
  mutate(metal_related = ifelse(!is.na(term.y),"metal related","unknown"))%>%
  group_by(ORF)%>%
  mutate(changing_in_any = ifelse(any(abs(log2(median_fc_flux)) > log2(1.5)),T,F ))%>%
  dplyr::select(ORF, metal_related, changing_in_any)%>%
  unique()%>%
  group_by(metal_related)%>%
  summarize(
    total_ORFs = length(ORF),
    num_ORF_changing_flux = sum(changing_in_any),
    frac_ORF_pathway_changingflux = num_ORF_changing_flux/total_ORFs)%>%
  ungroup()%>%
  filter(total_ORFs > 4)



### merge with KEGG ###

