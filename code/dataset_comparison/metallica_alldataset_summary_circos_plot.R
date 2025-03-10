
#`---
#`  Title: "metallica_alldataset_summary_circos_plot.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 25 May 2023
#`  Description: make circos plot to compare all hits
#`---

#############################################
### source paths functions and libraries  ###
#############################################

library(tidyr)
library(dplyr)
library(circlize)

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

###############################
### Input all dataset table ###
###############################

all_datasets <- read.csv(paste0(output_tables_dir,"/all_datasets_combined_significant_or_not.csv"),stringsAsFactors = F)

all_datasets <- merge(all_datasets, unique(GenProt_SGD[,c("ORF","Uniprot.Annotation.Score")]), by = "ORF")
all_datasets_for_poorlycharacterised <- all_datasets

all_datasets <- merge(all_datasets, GO_gset_MF_binding_metalwise, by = "ORF", all.x=T)
colnames(all_datasets)[which(colnames(all_datasets)=="term")] <- "metal_binder_specific"

all_datasets <- merge(all_datasets, GO_gset_MF_transporter_metalwise, by = "ORF", all.x = T)
colnames(all_datasets)[which(colnames(all_datasets) =="term")] <- "metal_transporter_specific"

all_datasets <- merge(all_datasets, GO_gset_MF_metalbinding_unspecific, by = "ORF", all.x=T)
colnames(all_datasets)[which(colnames(all_datasets)=="term")] <- "metal_binder_unspecific"



all_datasets <- all_datasets%>%
  mutate(metal_binding_anno = ifelse(is.na(metal_binder_specific) & is.na(metal_binder_unspecific),F,
                                     ifelse(!is.na(metal_binder_specific),"sp","unsp")),
         metal_transporter_anno = ifelse(is.na(metal_transporter_specific) ,F,
                                      ifelse(!is.na(metal_transporter_specific),"sp","unsp")),
         sig_in_metpertWTproteomics = ifelse(!is.na(Significant),
                                             ifelse(dataset =="metpertWTproteomics" & Significant, T, F),F))%>%
  dplyr::select(ORF,metal, dataset, Significant, Uniprot.Annotation.Score, metal_binding_anno, metal_transporter_anno,sig_in_metpertWTproteomics)%>%
  unique()

num_ds_sig_df <- all_datasets%>%
                 dplyr::select(ORF,dataset,metal, Significant)%>%
                 group_by(dataset,ORF)%>%
                 summarise(sig_in_ds = any(Significant))%>%
                 ungroup()%>%
                 group_by(ORF)%>%
                 summarise(num_ds_sig_in = sum(sig_in_ds))%>%
                 unique()

num_metals_sig_df <- all_datasets %>%
                      dplyr::select(ORF, dataset, metal, Significant) %>%
                      group_by(metal, ORF) %>%
                      summarise(sig_in_metal = any(Significant)) %>%
                      ungroup() %>%
                      group_by(ORF) %>%
                      summarise(num_metals_sig_in = sum(sig_in_metal)) %>%
                      unique()

# Now let's merge these two new dataframes with the original one to get the final result
all_datasets <- merge(unique(all_datasets[,c("ORF","Uniprot.Annotation.Score",
                                             "metal_binding_anno","metal_transporter_anno","sig_in_metpertWTproteomics")]),
                      num_ds_sig_df, by = "ORF")
all_datasets <- merge(all_datasets, num_metals_sig_df, by = "ORF")%>%

  # arrange datasets by decreasing Uniprot annotation score
                  arrange(desc(as.numeric(Uniprot.Annotation.Score))) %>% 
                  mutate(id = 1:n())


df_uniprot_anno <- all_datasets %>% 
  dplyr::select(ORF, Uniprot.Annotation.Score, id) %>%
  unique()%>%
  mutate(
    color = viridis_pal(option = "B",
                              begin = 0.2, 
                              end = 0.9)(max(as.numeric(Uniprot.Annotation.Score) + 1))[as.numeric(Uniprot.Annotation.Score) + 1]
  )%>%
  na.omit()

df_metalbindng_anno <- all_datasets %>% 
  dplyr::select(ORF, metal_binding_anno, id) %>%
  unique()%>%
  mutate(
    color = ifelse(metal_binding_anno == "sp","#287d8eb3",
                   ifelse(metal_binding_anno == "unsp",NA,NA)))%>%
  na.omit()
#287d8eb3
df_metaltransporter_anno <- all_datasets %>% 
                            dplyr::select(ORF, metal_transporter_anno, id) %>%
                            unique()%>%
                            mutate(
                              color = ifelse(metal_transporter_anno == "sp","#9b2964",
                                             ifelse(metal_transporter_anno == "unsp",NA,NA)))%>%
                            na.omit()
#f8870eb3

df_metperWTproteomics <- all_datasets%>%
                        dplyr::select(ORF, sig_in_metpertWTproteomics, id) %>%
                        unique()%>%
                        mutate(
                          color = ifelse(sig_in_metpertWTproteomics,"#453781" ,NA))%>%
                        na.omit()
                         

df_numsigin_ds_met <- all_datasets%>%
             dplyr::select(ORF, num_ds_sig_in,num_metals_sig_in, id)%>%
             unique()%>%
             mutate(color = viridis_pal(end = 0.9)(max(as.numeric(num_metals_sig_in) + 1))[as.numeric(num_metals_sig_in) + 1],
                    height = num_ds_sig_in/5
             )%>%
  na.omit()
     

# Plotting
pdf(paste0(plot_dir,"/circos_plot_all_dataset_comparison_5tracks.pdf"),width = 10,height=10)
circos.clear()
circos.par(cell.padding = c(0,0,0,0), track.margin = c(0,0), start.degree = 90,
           gap.degree = 0)
circos.initialize(factors = "a", xlim = c(1, max(df_uniprot_anno$id)))

# number of datasets significant in 
circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.15,
                       panel.fun = function(x, y) {
                         circos.segments(x0 = df_numsigin_ds_met$id, 
                                         x1 = df_numsigin_ds_met$id, 
                                         y0 = 0, 
                                         y1 = df_numsigin_ds_met$height, 
                                         col = df_numsigin_ds_met$color, 
                                         border = NA,
                                         lwd = 0.35)
                       })

# number of datasets significant in 
circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.15,
                       panel.fun = function(x, y) {
                         circos.segments(x0 = df_metperWTproteomics$id, 
                                         x1 = df_metperWTproteomics$id, 
                                         y0 = 0, 
                                         y1 = 1, 
                                         col = df_metperWTproteomics$color, 
                                         border = NA,
                                         lwd = 0.35)
                       })

# metal transporter annotation
circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.1,
                       panel.fun = function(x, y) {
                         circos.segments(x0 = df_metaltransporter_anno$id, 
                                         x1 = df_metaltransporter_anno$id, 
                                         y0 = 0.05, 
                                         y1 = 0.95, 
                                         col = df_metaltransporter_anno$color,
                                         border = NA,
                                         lwd = 0.35)
                       })

# metal binding annotation
circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.1,
                       panel.fun = function(x, y) {
                         circos.segments(x0 = df_metalbindng_anno$id, 
                                         x1 = df_metalbindng_anno$id, 
                                         y0 = 0.05, 
                                         y1 = 0.95, 
                                         col = df_metalbindng_anno$color, 
                                         border = NA,
                                         lwd = 0.35)
                       })

# Uniprot annotation score
circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.05,
                       panel.fun = function(x, y) {
                         circos.rect(xleft = df_uniprot_anno$id, 
                                     xright = df_uniprot_anno$id + 1, ybottom = 0, 
                                     ytop = 1, 
                                     col = df_uniprot_anno$color, border = NA)
                       })


dev.off()

## plot legend for circos plot

# Define the colors and labels for the legend
legend_colors <- c(
  "#287d8eb3", "#453781",   # Metal binding annotation
  "#f8870eb3", "#9b2964",   # Metal transporter annotation
  viridis_pal(option = "B", begin = 0.2, end = 0.9)(5), # Uniprot annotation score
  viridis_pal(end = 0.9)(10)[1:9] # Significant in datasets / metals
)

legend_labels <- c(
  "Metal binding annotation (unsp)",
  "Metal binding annotation (sp)",
  "Metal transporter annotation (unsp)",
  "Metal transporter annotation (sp)",
  paste("Uniprot annotation score", 1:5),
  paste("Significant in metals", 1:9)
)

# Create a data frame
df_legend <- data.frame(labels = legend_labels, colors = legend_colors)

pdf(paste0(plot_dir,"/circos_plot_legend.pdf"),width = 6,height = 10)
# Create a bar plot
ggplot(df_legend, aes(x = labels, y = 1, fill = colors)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  coord_flip()+
  theme_metallica()+
  labs(x = "", y = "", fill = "")
dev.off()


###################################################################
### Summarise poorly characterised proteins across all datasets ###
###################################################################

all_datasets_for_poorlycharacterised <- all_datasets_for_poorlycharacterised %>%
                                        mutate(poorly_characterised = ifelse(as.numeric(Uniprot.Annotation.Score) < 3,
                                                                             TRUE, FALSE),
                                               dataset = factor(dataset, levels = c("metpertWTproteomics",
                                                                                    "y5kmetalspecific",
                                                                                    "metdepKOgrowth",
                                                                                    "KOmetallomics",
                                                                                    "OEmetallomics")))%>%
                                        arrange(dataset)

## Look at overlaps between poorly characterised proteins that are signficant in each dataset

poorly_characterised_df <- all_datasets_for_poorlycharacterised %>%
  filter(poorly_characterised == TRUE  & Significant == TRUE)

poorly_characterised_list <- lapply(unique(poorly_characterised_df$dataset), function(x) {
  unique(poorly_characterised_df$ORF[poorly_characterised_df$dataset == x])
})

names(poorly_characterised_list) <- unique(poorly_characterised_df$dataset)

pdf(paste0(plot_dir,"/venn_poorly_characterised_across_datasets.pdf"), width = 10, height = 10)
venn::venn(poorly_characterised_list, zcolor = "style", box = FALSE, ilcs=2.5)
dev.off()


## summarise numbers of poorly characterised proteins

print(paste("total significant poorly characterised proteins in any dataset:",
      length(unique(poorly_characterised_df$ORF))))


# make table with counts of datasets where each ORF is significant
nds_count_table_ORF <- poorly_characterised_df %>%
  group_by(ORF) %>%
  summarise(n_datasets = n_distinct(dataset)) %>%
  ungroup()

write.csv(nds_count_table_ORF, 
          paste0(output_tables_dir,"/poorlycharacterised_ORFs_numds_sig_in.csv",row.names = F))

# Count how many proteins are significant in exactly 1, 2, 3, 4, and 5 datasets
nds_count_table_ORF_summary <- nds_count_table_ORF %>%
  group_by(n_datasets) %>%
  summarise(n_proteins = n()) %>%
  ungroup()%>%
# Calculate cumulative number of proteins
  arrange(desc(n_datasets)) %>%
  mutate(cumulative_n_proteins = cumsum(n_proteins))
write.csv(nds_count_table_ORF_summary, 
          paste0(output_tables_dir,"/poorlycharacterised_numORFs_sig_in_numds.csv"),row.names = F)


# make table with counts of datasets where each ORF-metal combination is significant
nds_count_table_metal_ORF <- poorly_characterised_df %>%
  group_by(ORF, metal) %>%
  summarise(n_datasets = n_distinct(dataset)) %>%
  ungroup()
write.csv(nds_count_table_metal_ORF, 
          paste0(output_tables_dir,"/poorlycharacterised_ORFmetalpairs_numds_sig_in.csv"),row.names = F)

# Count how many ORF-metal combinations are significant in exactly 1, 2, 3, 4, and 5 datasets
nds_count_summary_metal_ORF <- nds_count_table_metal_ORF %>%
  group_by(n_datasets) %>%
  summarise(n_combinations = n()) %>%
  ungroup()%>%
  # Calculate cumulative number of combinations
  arrange(desc(n_datasets)) %>%
  mutate(cumulative_n_combinations = cumsum(n_combinations))
write.csv(nds_count_summary_metal_ORF, 
          paste0(output_tables_dir,"/poorlycharacterised_numORFmetalpairss_sig_in_numds.csv"),row.names = F)



### count number of poorly characterised ORFs that are significant in each metal condition in each dataset
num_poorlychar_sigORF_permetalds <- poorly_characterised_df %>%
  mutate(dataset = factor(dataset, levels = c("metdepKOgrowth",
                                              "KOmetallomics",
                                              "OEmetallomics",
                                              "y5kmetalspecific",
                                              "metpertWTproteomics")))%>%
  group_by(metal, dataset, ORF) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(metal, dataset) %>%
  summarise(total_ORFs = n_distinct(ORF)) %>%
  ungroup()

### Write table with poorly characterised ORF-metal pairs that are significant in 3 or more datasets
poorly_charac_in3ormoreds <- poorly_characterised_df %>%
                            mutate(ORF_metal = paste(ORF, metal))%>%
                            filter(ORF_metal %in%  paste(filter(nds_count_table_metal_ORF, n_datasets >2)$ORF,
                                                filter(nds_count_table_metal_ORF, n_datasets >2)$metal) )
write.csv(poorly_charac_in3ormoreds, 
          paste0(output_tables_dir,"/poorlycharacterised_ORFmetal_sigin3ormore.csv"),row.names = F)


pdf(paste0(plot_dir,"/poorlycharacterised_ORFs_significant_datasets.pdf"), width = 12,
    height = 6)
ggplot(num_poorlychar_sigORF_permetalds, aes(x = metal, y = dataset, fill = total_ORFs)) +
  geom_tile() +
  geom_text(aes(label = total_ORFs), color = "white", size = 5) +
  scale_fill_viridis() +
  theme_metallica() +
  labs(x = " ", y = " ", fill = "poorly\ncharacterised\nsignificant\nORFs")
dev.off()



