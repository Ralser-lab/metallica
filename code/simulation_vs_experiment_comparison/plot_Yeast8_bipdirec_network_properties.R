#`---
#`  Title: "plot_Yeast8_bipdirec_network_properties.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 Nov 2022
#`  Description: Script to visualise network properties in context of metal anno of yeast8 bipartite directed network graph
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



centrality_metrics = read.csv(paste0(metpert_sim_vs_exp_comparison_dir,
                                     "/metabolic_network_properties/Yeast8_bip_direc_revasdouble_network_centrality_metrics.csv"),
                              stringsAsFactors = F)%>%
                      reshape2::melt(id.vars = c("Centrality.Measure", "Significant","p.value"))%>%
                      mutate(type = ifelse(grepl("Non.Metal",variable),"no metal","metal"),
                             variable = ifelse(grepl("Mean",variable), "Mean",
                                             ifelse(grepl("SEM",variable),"SEM",
                                                    ifelse(grepl("n[.]",variable),"N",NA))))%>%
                      reshape2::dcast(Centrality.Measure + Significant + p.value + type ~ variable, value.var = "value")


pdf(paste0(plot_dir,"/Yeast8_bipdirec_network_centrality_metrics_metnomet.pdf"),width = 13,height = 7)
ggplot(centrality_metrics,
       aes(x = type, 
           y = Mean,
           group = type,
           color = type,
           fill = type))+
  facet_wrap("Centrality.Measure",scales = "free",ncol = 5)+
  geom_bar(stat = "identity", width = 0.6,
           alpha = 0.5)+
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.25, colour = "black", linewidth = 0.2)+
  geom_text(aes( y = Mean /2,
                label = ifelse(p.value < 0.001, "***",
                                                          ifelse(p.value < 0.01, "**",
                                                                 ifelse(p.value < 0.05,"*", "ns")))),
            colour = "black",size = 3)+
  geom_text(aes( y = Mean/4,
               label = paste0("n=",`N`)),
          colour = "black", size = 3)+
  scale_colour_manual(values = c("#440154","#21918c"))+
  scale_fill_manual(values = c("#440154","#21918c"))+
  theme_metallica()+
    labs(x = "", fill = "", colour = "")+
  theme(legend.position = "bottom")
dev.off()



##################
### metal wise ###
##################


centrality_metrics_metalwise <- read.csv(paste0(metpert_sim_vs_exp_comparison_dir,
                                                 "/metabolic_network_properties/Yeast8_bipdir_network_centrality_metrics_metalwise.csv"),
                                          stringsAsFactors = FALSE) %>%
  group_by(Measure) %>%
  mutate(Rank = rank(-Mean)) %>%
  ungroup()  # to remove the grouping

max_rank <- max(centrality_metrics_metalwise$Rank, na.rm = TRUE)

pdf(paste0(plot_dir,"/centrality_metrics_summary.pdf"), width = 10,height = 7)
ggplot(centrality_metrics_metalwise,
       aes(x = reorder(Metal,-Rank),
           y = Measure,
           fill = Rank))+
  geom_tile(colour = "white")+
  geom_text(aes(label = N))+
  scale_fill_viridis_c(breaks = c(1, 5, max_rank)) +  
  
  theme_metallica()

dev.off()








