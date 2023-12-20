
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
  dplyr::select(ORF,dataset,Significant,metal_binder,metal_transporter,Uniprot.Annotation.Score)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(num_sig = sum(Significant,na.rm=T))%>%
  ungroup()%>%
  dplyr::select(-dataset)%>%
  unique()%>%
  arrange(desc(as.numeric(Uniprot.Annotation.Score))) %>% 
  mutate(id = 1:n())



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

# prepare num ds significant in 
df_numsig <- all_datasets %>% 
  dplyr::select(ORF, num_sig, id) %>%
  unique()%>%
  mutate(
    id = 1:n(),
    color = viridis_pal(option = "A",
                        begin = 0.2, 
                        end = 0.9)(max(as.numeric(num_sig) + 1))[as.numeric(num_sig) + 1]
  )
    
df_metal_binder =  all_datasets %>%
      dplyr::select(ORF, id, metal_binder) %>%
      mutate(color = ifelse(metal_binder, "black", NA))

df_metal_transporter=  all_datasets %>%
      dplyr::select(ORF, id, metal_transporter) %>%
      mutate(color = ifelse(metal_transporter,"maroon", NA))


# Plotting
circos.clear()
circos.par(cell.padding = c(0,0,0,0), track.margin = c(0,0), start.degree = 90, gap.degree = 0)
circos.initialize(factors = "a", xlim = c(1, max(df_annotation$id)))


# number of datasets ORF was a hit in 

circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.1,
                       panel.fun = function(x, y) {
                         circos.segments(df_numsig$id, 0.05, df_numsig$id, 0.95,
                                         straight = T,
                                         col = df_numsig$color,
                                         lwd = 0.3)
                         
                       })

# number 

circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.1,
                         panel.fun = function(x, y) {
                           circos.segments(df_metal_transporter$id, 0.05, df_metal_transporter$id, 0.95,
                                           straight = T,
                                           col = df_metal_transporter$color,
                                           lwd = 0.3)
                           
                         })

# Uniprot Annotation Score
df_annotation <- na.omit(df_annotation) # remove rows with NA values

circos.trackPlotRegion(factors = "a", ylim = c(0,1), track.height = 0.05,
                       panel.fun = function(x, y) {
                         circos.rect(xleft = df_annotation$id, 
                                     xright = df_annotation$id + 1, ybottom = 0, 
                                     ytop = 1, 
                                     col = df_annotation$color, border = NA)
                       })

