#####################################################################
####          Script to plot results of ensemble clustering      ####
#####################################################################

#`---
#`  Title: "metperWTproteomics_visualise_ensembleclustering_summary.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 12 July 2023
#`  Description: Script to visualise summary of ensemble clustering using Sankey Plots
#`---

#############################################
### source paths functions and libraries  ###
#############################################
library(plotly)
library(visNetwork)
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/ensemble_clustering")
dir.create(paste0(plot_dir,"/summary/"), recursive = T)

input_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering/clustering_results/summary")

all_ensemble_clustering_annotated <- read.csv(paste0(input_tables_dir,"/ensemble_clustering_allresults_annotated.csv"),stringsAsFactors = F)%>%
                                     dplyr::select(-Uniprot_ID, -Gene)%>%
                                     mutate(poorly_characterised = ifelse(Uniprot.Annotation.Score < 3, T, F))
         
##########################################
### prepare dataframe for Sankey plots ###
##########################################

gstypes = c("GOslim_BP","GOslim_CC","GOslim_MF","KEGG","GO_MF_metalspecific_metalbinding")


for(gs in 1:length(gstypes)){

all_ensemble_clustering_annotated_gsettype <- all_ensemble_clustering_annotated%>%
                                              filter( (Gset_type_metalwise == gstypes[gs] | is.na(Gset_type_metalwise)) & 
                                                     (Gset_type_allmetal == gstypes[gs] | is.na(Gset_type_allmetal)) &
                                                      !is.na(cluster_allmetal))%>%
                                              mutate(Gset_term_enriched_metalwise = ifelse(is.na(Gset_term_enriched_metalwise),"no enrichment",
                                                                                           Gset_term_enriched_metalwise),
                                                     Gset_term_enriched_allmetal = ifelse(is.na(Gset_term_enriched_allmetal),"no enrichment",
                                                                                          Gset_term_enriched_allmetal))

## prepare df metalwise clustering enrichments --> metal wise clusters

df_mwenrch2mwclust <- all_ensemble_clustering_annotated_gsettype%>%
                      dplyr::select(ORF,
                                    metal,
                                    Gset_term_enriched_metalwise,
                                    cluster_metalwise,
                                    poorly_characterised)%>%
                      unique()%>%
                      group_by(metal, cluster_metalwise,
                               poorly_characterised,Gset_term_enriched_metalwise)%>%
                      summarise(value = length(ORF))%>%
                      ungroup()%>%
                      mutate(source  = paste("mw_enrich",Gset_term_enriched_metalwise,sep="_"),
                             target = paste("mw_cluster",metal, cluster_metalwise,sep="_"),
                             link_colour = ifelse(poorly_characterised, "black",
                                    sapply(metal, function(x) paste0(colkey_Ele[[x]], "60"))))%>%
                      dplyr::select(source, target, link_colour, value)

                                    
df_mclust2allmetclust <- all_ensemble_clustering_annotated_gsettype %>%
                        mutate(cluster_metalwise = ifelse(is.na(cluster_metalwise), "no mw cluster", cluster_metalwise),
                               cluster_allmetal = ifelse(is.na(cluster_allmetal), "no allm cluster", cluster_allmetal)) %>%
                        dplyr::select(ORF,metal,cluster_metalwise, cluster_allmetal, poorly_characterised) %>%
                        unique() %>%
                        group_by(metal,cluster_metalwise, cluster_allmetal, poorly_characterised) %>%
                        summarise(value = length(ORF)) %>%
                        ungroup() %>%
                        arrange(metal,cluster_metalwise)%>%
                        mutate(source  = paste("mw_cluster",metal, cluster_metalwise,sep="_"),
                               target = paste("allm_cluster",as.character(cluster_allmetal),sep="_"),
                               link_colour = ifelse(poorly_characterised, "black",
                                                    sapply(metal, function(x) paste0(colkey_Ele[[x]], "20")))) %>%
                        dplyr::select(source, target, link_colour, value)

  
df_allmetal2gset <- all_ensemble_clustering_annotated_gsettype %>%
                    dplyr::select(ORF, cluster_allmetal, Gset_term_enriched_allmetal, poorly_characterised) %>%
                    unique() %>%
                    group_by(cluster_allmetal, Gset_term_enriched_allmetal) %>%
                    ungroup() %>%
                    group_by(cluster_allmetal,poorly_characterised, Gset_term_enriched_allmetal) %>%
                    mutate(value = length(Gset_term_enriched_allmetal)) %>%
                    ungroup() %>%
                    mutate(source  = paste("allm_cluster",as.character(cluster_allmetal),sep= "_"),
                           target = paste("allm_enrichm",Gset_term_enriched_allmetal,sep="_"),
                            link_colour = ifelse(poorly_characterised, "black",
                                                    "beige"))%>%
                    dplyr::select(source, target, link_colour, value)%>%
                    unique()

                    
# Combine the three dataframes
df_for_sankey <- rbind(df_mwenrch2mwclust, df_mclust2allmetclust, df_allmetal2gset)%>%
                 filter(! (grepl("no enrichment", source)|
                          grepl("no enrichment", target)))

# Generate all nodes
buffer <- 0.001 # adjust this value to get the desired spacing between y-coordinates
sankey_nodes <- data.frame(name = c(as.character(df_for_sankey$source), 
                                    as.character(df_for_sankey$target)) %>% 
                             unique()) %>%
  # Add node positions 
  mutate(x_coord = ifelse(grepl("mw_enrich", name), 0.1,
                          ifelse(grepl("mw_cluster", name), 0.4,
                                 ifelse(grepl("allm_cluster", name), 0.6, 0.9)))) %>%
  mutate(n = 1) %>%
  # y coord
  group_by(x_coord) %>%
  mutate(y_cood = length(x_coord),
         counter = cumsum(n)) %>%
  ungroup() %>%
  mutate(y_cood = (counter + buffer) / (y_cood + 2 * buffer)) %>% # add buffer to both the counter and y_cood to maintain the relative positioning
  dplyr::select(-n, -counter)

df_for_sankey$IDsource = match(df_for_sankey$source, sankey_nodes$name) - 1 
df_for_sankey$IDtarget = match(df_for_sankey$target, sankey_nodes$name) - 1

# Create a function to assign colors based on node names
convert_nodenames_to_hexcolour <- function(x) {
  colour = ifelse(grepl("Ca_",x),colkey_Ele[["Ca"]],
                  ifelse(grepl("Cu_",x),colkey_Ele[["Cu"]],
                         ifelse(grepl("Fe_",x),colkey_Ele[["Fe"]],
                                ifelse(grepl("K_",x),colkey_Ele[["K"]],
                                       ifelse(grepl("Mg_",x),colkey_Ele[["Mg"]],
                                              ifelse(grepl("Mn_",x),colkey_Ele[["Mn"]],
                                                     ifelse(grepl("Mo_",x),colkey_Ele[["Mo"]],
                                                            ifelse(grepl("Na_",x),colkey_Ele[["Na"]],
                                                                   ifelse(grepl("Zn_",x),colkey_Ele[["Zn"]],
                                                                          "beige")))))))))
  return(colour)
  
}

node_color_vector <- as.character(lapply(sankey_nodes$name, convert_nodenames_to_hexcolour))

# Convert node names to readable format
sankey_nodes <- sankey_nodes %>%
  mutate(name = gsub("mw_enrich_", "", name),
         name = gsub("mw_cluster_", "", name),
         name = gsub("allm_cluster_", "", name),
         name = gsub("allm_enrichm_", "", name))

fig <- plotly::plot_ly(
  type = "sankey",
  node = list(
    label = sankey_nodes$name,
    x = sankey_nodes$x_coord,
    y = sankey_nodes$y_cood,
    color = node_color_vector,
    pad = 20,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  link = list(
    source = df_for_sankey$IDsource,
    target = df_for_sankey$IDtarget,
    value =  df_for_sankey$value,
    color = df_for_sankey$link_colour
  )
)

fig <- fig %>% layout(
  title = paste(gstypes[gs]), # replace GeneSetType with the title you want
  font = list(
    size = 15,
    family = "Times"
  )
)

plotly::export(fig, 
               file = paste0(plot_dir,"/summary/sankey_ensembleclustering",gstypes[gs],".pdf"))

}
  