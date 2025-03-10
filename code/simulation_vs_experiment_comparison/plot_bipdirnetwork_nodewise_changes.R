#`---
#`  Title: "plot_bipdirnetwork_nodewise_changes.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 10 Nov 2022
#`  Description: Script to visualise flux and protein levels changes at each node distance
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


########################################################################################
### read in node-wise flux and enzyme abundance changes data -- calculated in python ###
########################################################################################

## shortest distance to any metal node file 

node_dist_sig_df <- read.csv(paste0(metpert_sim_vs_exp_comparison_dir,"/metabolic_network_properties/node_distance_wise_flux_env_cell_significant_grouped.csv"),
                             stringsAsFactors = F)%>%
                    reshape2::melt(id.vars = c("distance_group","total_reactions"),
                                   value.name = "num_sig",
                                   variable.name = "type")%>%
                    mutate(frac_sig = num_sig/total_reactions,
                           distance_group = factor(distance_group,
                                                   levels = c("0","2","4","6",">6")),
                           type = gsub("_significant_count","",type),
                           type = factor(type, levels = c("flux",
                                                          "env",
                                                          "cell")))

pdf(paste0(plot_dir,"/percent_fluxenvcell_ScGEM_nodedistwise.pdf"), width = 8,height = 6)

ggplot(node_dist_sig_df)+
geom_bar(aes(x = distance_group,
             y = total_reactions,
             colour = type),
         fill = NA, stat = "identity", position = "dodge",
         width = 0.5,
         linewidth = 0.25)+
  geom_bar(aes(x = distance_group,
               y = num_sig,
               colour = type,
               fill = type), 
           alpha = 0.6,
           linewidth = 0.25,
           stat = "identity", position = "dodge2",
           width = 0.5)+
  geom_text(aes(x = distance_group,
                y = num_sig+10,
                group = type,
                colour = type,
                label = paste0(100*round(frac_sig,2),"%") ),
            size = 5)+
  theme_metallica()+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)+
  labs(x= "distance nearest metal annotated reaction",
       y = "number of reactions",
       fill = "",
       colour = "")


## bar height proportional to % of nodes

ggplot(node_dist_sig_df)+
  geom_bar(aes(x = distance_group,
               y = frac_sig,
               fill = type),
           colour = "white",
            stat = "identity", position = "dodge",
           width = 0.5,
           linewidth = 0.5)+
  geom_text(aes(x = distance_group,
                y = frac_sig + 0.02,
                group = type,
                colour = type,
                label = total_reactions ),
            size = 5)+
  theme_metallica()+
  ylim(0, 0.52)+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)+
  labs(x= "distance nearest metal annotated reaction",
       y = "fraction of significant reactions",
       fill = "",
       colour = "")

dev.off()


###################################################
### node distance and sig calculation per metal ###
###################################################

# metal wise number of nodes


metal_wise_nodedistance_sig <-read.csv(paste0(metpert_sim_vs_exp_comparison_dir,
                                              "/metabolic_network_properties/node_distance_wise_metal_wise_flux_env_cell_signficant_grouped.csv"),
                                       stringsAsFactors = F)%>%
  reshape2::melt(id.vars = c("distance_group","total_reactions","Metal"),
                 value.name = "num_sig",
                 variable.name = "type")%>%
  mutate(frac_sig = num_sig/total_reactions,
         distance_group = factor(distance_group, levels = c("0.0","2.0","4.0","6.0",">6")),
         Metal = gsub("unspecific","unsp",Metal),
         Metal = factor(Metal, levels = c("Ca","Cu","Fe","K","Mg","Mn","Na","Zn","unsp")),
         type = gsub("_significant_count","",type),
         type = factor(type, levels = c("flux","env","cell")))

pdf(paste0(plot_dir,"/percent_fluxenvcell_Sc_GEM_nodedistancewise_metalwise.pdf"),width = 15,height =15)
ggplot(metal_wise_nodedistance_sig,
       aes(x = Metal))+
  geom_bar(aes(
               y = total_reactions,
               colour = type),
           fill = NA, stat = "identity", position = "dodge",
           width = 0.5,
           linewidth = 0.25)+
  geom_bar(aes(x = Metal,
               y = num_sig,
               colour = type,
               fill = type), 
           alpha = 0.6,
           linewidth = 0.25,
           stat = "identity", position = "dodge2",
           width = 0.5)+
  facet_wrap(c("distance_group"),scales = "free", ncol = 2)+  
  theme_metallica()+
  labs(x = "")+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)
dev.off()
