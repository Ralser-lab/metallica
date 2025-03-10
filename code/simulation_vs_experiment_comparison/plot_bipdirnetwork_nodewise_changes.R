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
  #geom_text(aes(x = Metal,
   #            y = num_sig+10,
    #          group = type,
     #        colour = type,
      #      label = paste0(100*round(frac_sig,2),"%") ),
      # size = 5)+
  facet_wrap(c("distance_group"),scales = "free", ncol = 2)+  
  theme_metallica()+
  labs(x = "")+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)
dev.off()

##################################################################
##################################################################
##################################################################

### Example of metal node to non-metal node signal propagation ###

##################################################################
##################################################################
##################################################################

## read in df created using igraph in jupyter notebook "Yeast8_dir-bip-graph_create_network_calc_centrality_distances2metalnodes"


metalnode_2_signode <- read.csv(paste0(metpert_sim_vs_exp_comparison_dir,"/metabolic_network_properties/nearest_metal_nodes_to_significant_nodes.csv"),
                                stringsAsFactors = F)
                       ## Filter for distance = 4 -- we just want to look for one nice examples
metalnode_2_signode_dist4 <- metalnode_2_signode%>%
                              filter(Distance_from_Nearest_Metal_ID_Node == 4)%>%
                              unique()

metalnode_2_signode_dist2 <- metalnode_2_signode%>%
  filter(Distance_from_Nearest_Metal_ID_Node == 2)%>%
  unique()

metalnode_2_signode_dist0 <- metalnode_2_signode%>%
  filter(Distance_from_Nearest_Metal_ID_Node == 1)%>%
  unique()

# Create an empty dataframe to store the results
result_df <- data.frame(
  Significant_Node_ID = character(),
  Nearest_Metal_Node_ID = character(),
  Distance = integer(),
  ORFs_Significant_Node_Dist4 = character(),
  ORFs_Significant_Node_Dist2 = character(),
  ORFs_Significant_Node_Dist0 = character(),
  Reaction_IDs_Significant_Node_Dist4 = character(),
  Reaction_IDs_Significant_Node_Dist2 = character(),
  Reaction_IDs_Significant_Node_Dist0 = character(),
  stringsAsFactors = FALSE
)

# Iterate through each unique Significant Node ID in the distance 4 dataframe
for (sign_node_id in unique(metalnode_2_signode_dist4$Significant_Node_ID)) {
  
  # Find the corresponding Nearest Metal Node ID for the current Significant Node ID
  nearest_metal_id <- metalnode_2_signode_dist4$Nearest_Metal_Node_ID[metalnode_2_signode_dist4$Significant_Node_ID == sign_node_id]
  
  # Filter the distance 2 dataframe to find Significant Node IDs for the same Nearest Metal Node ID
  dist2_significant_nodes <- metalnode_2_signode_dist2$Significant_Node_ID[metalnode_2_signode_dist2$Nearest_Metal_Node_ID == nearest_metal_id]
  
  # Filter the distance 0 dataframe to find Significant Node IDs for the same Nearest Metal Node ID
  dist0_significant_nodes <- metalnode_2_signode_dist0$Significant_Node_ID[metalnode_2_signode_dist0$Nearest_Metal_Node_ID == nearest_metal_id]
  
  # Get the ORFs and Reaction IDs annotated to the Significant Node ID at distance 4
  orfs_dist4 <- metalnode_2_signode_dist4$ORFs_in_Significant_Node[metalnode_2_signode_dist4$Significant_Node_ID == sign_node_id]
  reaction_ids_dist4 <- metalnode_2_signode_dist4$Reaction_IDs[metalnode_2_signode_dist4$Significant_Node_ID == sign_node_id]
  
  # Get the ORFs and Reaction IDs annotated to the Significant Node IDs at distance 2 and distance 0
  orfs_dist2 <- metalnode_2_signode_dist2$ORFs_in_Significant_Node[metalnode_2_signode_dist2$Significant_Node_ID %in% dist2_significant_nodes]
  reaction_ids_dist2 <- metalnode_2_signode_dist2$Significant_Node_ID[metalnode_2_signode_dist2$Significant_Node_ID %in% dist2_significant_nodes]
  orfs_dist0 <- metalnode_2_signode_dist0$ORFs_in_Significant_Node[metalnode_2_signode_dist0$Significant_Node_ID %in% dist0_significant_nodes]
  reaction_ids_dist0 <- metalnode_2_signode_dist0$Significant_Node_ID[metalnode_2_signode_dist0$Significant_Node_ID %in% dist0_significant_nodes]
  
  # Add the results to the dataframe
  result_df <- rbind(result_df, data.frame(
    Significant_Node_ID = sign_node_id,
    Nearest_Metal_Node_ID = nearest_metal_id,
    Distance = 4,  # Distance 4
    ORFs_Significant_Node_Dist4 = paste(orfs_dist4, collapse = ", "),  # Combine ORFs as a string
    ORFs_Significant_Node_Dist2 = paste(orfs_dist2, collapse = ", "),  # Combine ORFs as a string
    ORFs_Significant_Node_Dist0 = paste(orfs_dist0, collapse = ", "),  # Combine ORFs as a string
    Reaction_IDs_Significant_Node_Dist4 = paste(reaction_ids_dist4, collapse = ", "),  # Combine Reaction IDs as a string
    Reaction_IDs_Significant_Node_Dist2 = paste(reaction_ids_dist2, collapse = ", "),  # Combine Reaction IDs as a string
    Reaction_IDs_Significant_Node_Dist0 = paste(reaction_ids_dist0, collapse = ", "),  # Combine Reaction IDs as a string
    stringsAsFactors = FALSE
  ))
}

# Reset row names of the result dataframe
row.names(result_df) <- NULL

# Print the result dataframe
print(result_df)

