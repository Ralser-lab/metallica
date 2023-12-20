
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


all_datasets <- merge(all_datasets, GO_gset_MF_binding_metalwise, by = "ORF", all.x=T)
colnames(all_datasets)[which(colnames(all_datasets)=="term")] <- "metal_binder"

all_datasets <- merge(all_datasets, GO_gset_MF_transporter_metalwise, by = "ORF", all.x = T)
colnames(all_datasets)[which(colnames(all_datasets) =="term")] <- "metal_transporter"

all_datasets <- all_datasets%>%
                mutate(metal_binder = ifelse(is.na(metal_binder),F,T),
                       metal_transporter = ifelse(is.na(metal_transporter),F,T))%>%
                group_by(ORF,dataset)%>%
                mutate(Significant = any(Significant))%>%
                ungroup()%>%
                dplyr::select(ORF,Significant,dataset,metal_binder,metal_transporter,Uniprot.Annotation.Score)%>%
                unique()%>%
                reshape2::dcast(ORF+metal_binder+metal_transporter+Uniprot.Annotation.Score ~ dataset, value.var = "Significant")%>%
                arrange(desc(as.numeric(Uniprot.Annotation.Score))) %>% 
                mutate(id = 1:n())%>%
                reshape2::melt(id.vars = c("ORF","Uniprot.Annotation.Score","id"), variable.name = "dataset_type")
               
                                

# prepare annotation dataset

df_annotation <- all_datasets %>% 
  dplyr::select(ORF, Uniprot.Annotation.Score, id) %>%
  unique()%>%
  mutate(
    id = 1:n(),
    color = viridis_pal(option = "B",
                              begin = 0.2, 
                              end = 0.9)(max(as.numeric(Uniprot.Annotation.Score) + 1))[as.numeric(Uniprot.Annotation.Score) + 1]
  )
# prepare list of dataframes with correct colour

ds_tracks <- c("metpertWTproteomics", "y5kmetalspecific",
              "OEmetallomics","KOmetallomics",
              "metdepKOgrowth",  "metal_binder","metal_transporter")
ds_list <- list()

for(i in 1:length(ds_tracks)){
  
  if(ds_tracks[i] == "metal_binder"){
    
    df =  all_datasets %>%
      filter(dataset_type == ds_tracks[i]) %>%
      dplyr::select(ORF, id, value) %>%
      mutate(color = ifelse(value, "black", NA))
    
  }else if(ds_tracks[i]=="metal_transporter"){
    df =  all_datasets %>%
      filter(dataset_type == ds_tracks[i]) %>%
      dplyr::select(ORF, id, value) %>%
      mutate(color = ifelse(value,"maroon", NA))
  }
  else{
  df =  all_datasets %>%
    filter(dataset_type == ds_tracks[i]) %>%
    dplyr::select(ORF, id, value) %>%
    mutate(color = ifelse(value, colkey_dataset_circos[as.character(ds_tracks[i])],
                          ifelse(is.na(value),NA,NA)))
  }
  ds_list[[i]] = df
}

names(ds_list) <- ds_tracks

pdf(paste0(plot_dir,"/circos_plot_all_dataset.pdf"),width = 10,height= 10)
# Plotting
circos.clear()
circos.par(cell.padding = c(0,0,0,0), track.margin = c(0,0), start.degree = 90, gap.degree = 0)
circos.initialize(factors = "a", xlim = c(1, max(df_annotation$id)))


# Datasets
for(i in 1:length(ds_list)) {
  
  dataset <- ds_list[i]
  df <- na.omit(ds_list[[i]])

  circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.1,
                         panel.fun = function(x, y) {
                           circos.segments(df$id, 0.05, df$id, 0.95,
                                           straight = T,
                                           col = df$col,
                                           lwd = 0.3)
                           
                         })
}

# Uniprot Annotation Score
df_annotation <- na.omit(df_annotation) # remove rows with NA values

circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.05,
                       panel.fun = function(x, y) {
                         circos.rect(xleft = df_annotation$id, 
                                     xright = df_annotation$id + 1, ybottom = 0, 
                                     ytop = 1, 
                                     col = df_annotation$color, border = NA)
                       })
dev.off()
