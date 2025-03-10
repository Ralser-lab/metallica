
#`---
#`  Title: "summarise_y5kpubphenomics_data.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 26 April 2023
#`  Description: explored the phenomics dataset by Baryhs.. et al 2022 and summarises how many metal related screens and phenotypes are available
#`---

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#############################################
### source paths functions and libraries  ###
#############################################


# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
#source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

plot_dir <- paste0(published_dataset_dir,"/phenomics/output/plots")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(published_dataset_dir,"/phenomics/output/tables")
dir.create(outputtables_dir, recursive = T)

phenomics_data <- read.delim(paste0(published_dataset_dir,"/phenomics/yp_screens_haphom_20221025.txt"))


# List of conditions
condition_list <- c( "calcium", "copper", "iron", "potassium", "magnesium", 
                    "manganese", "molybdenum", "sodium", "zinc", "chelator", "batho","desferrioxamine", "BPS","EDTA" ,"EGTA")

phenomics_data <- phenomics_data %>%
  mutate(metal_related = ifelse(str_detect(conditionset, paste(condition_list, collapse="|")), T, F),
         metal = case_when(
           str_detect(conditionset, "chelator") ~ "Chelator",
           str_detect(conditionset, "batho") ~ "Chelator",
           str_detect(conditionset, "desferrioxamine") ~ "Chelator",
           str_detect(conditionset, "EDTA") ~ "Chelator",
           str_detect(conditionset, "EGTA") ~ "Chelator",
           str_detect(conditionset, "BPS") ~ "Chelator",
           str_detect(conditionset, "calcium") ~ "Ca",
           str_detect(conditionset, "copper") ~ "Cu",
           str_detect(conditionset, "iron") ~ "Fe",
           str_detect(conditionset, "potassium") ~ "K",
           str_detect(conditionset, "magnesium") ~ "Mg",
           str_detect(conditionset, "manganese") ~ "Mn",
           str_detect(conditionset, "molybdenum") ~ "Mo",
           str_detect(conditionset, "sodium") ~ "Na",
           str_detect(conditionset, "zinc") ~ "Zn",
           TRUE ~ "None"
         ))%>%
  mutate(metal_related = ifelse(conditionset %in% c("sodium acetate [125 mM]","sodium acetate [150 mM]",
                                                    "sodium chloride [0.5 M], microgravity, time [21 gen]","sodium chloride [0.5 M], microgravity, time [14 gen]",
                                                    "sodium chloride [0.4 M], temperature [39 C], time [72 h]",
                                                    "sodium chloride [0.4 M], temperature [39 C], time [48 h]",
                                    "sodium arsenite","sodium arsenite [0.075-1 mM]","sodium arsenite [0.2 mM]","sodium arsenite [0.4 mM]","sodium arsenite [0.5 mM], time [24 h]",
                                    "sodium arsenite [0.5 mM], time [48 h]","sodium arsenite [1 mM], time [24 h]","sodium arsenite [1 mM], time [48 h]","sodium arsenite [1.25 mM]",
                                    "sodium arsenite [1.5 mM], time [24 h]","sodium arsenite [1.5 mM], time [48 h]","sodium arsenite [150 uM], time [15 gen]","sodium arsenite [150 uM], time [5 gen]",
                                    "sodium arsenite [300 uM], time [15 gen]","sodium arsenite [300 uM], time [5 gen]","sodium arsenite [400 uM]","sodium arsenite [625 uM], time [o/n + 20 gen]",
                                    "sodium arsenite [75 uM], time [15 gen]","sodium arsenite [75 uM], time [5 gen]","sodium chloride [0.6 M], temperature [39 C], time [48 h]","sodium chloride [0.6 M], temperature [39 C], time [72 h]",
                                    "sodium hydroxide [8 mM]", "sodium hydroxide [10 mM]","sodium metabisulfite [12 mM]","sodium pyrithione [8.35 uM]","sodium selenide [1 uM], time [16 h]","sodium selenide [1 uM], time [27 h]",
                                    "sodium selenide [2 uM], time [16 h]","sodium selenide [2 uM], time [27 h]","sodium selenite [0.05-0.5 mM]","sodium tetradecyl sulfate [13.37 uM]",
                                    "2-(N-morpholino)ethanesulfonic acid [20 mM], EGTA [20 mM], sodium chloride [4 M]","cantharidin disodium [100 uM], time [o/n + 20 gen]",
                                    "FK506 [0.05 ug/ml], sodium chloride [0.6 M], time [o/n + 20 gen]","FK506 [0.1 ug/ml], sodium chloride [0.6 M], time [o/n + 20 gen]","hydrogen peroxide [1-1.2 mM], sodium chloride [0.7 M]"),
                F, metal_related),
         metal = ifelse(metal_related,metal,"None"))


pdf(paste0(plot_dir,"/pie_y5kphenomics_metal_related.pdf"),width = 6,height=4)
# Bar chart to show counts of TRUE vs FALSE in metal_related
ggplot(phenomics_data, aes(x = metal_related, fill = metal_related)) +
  geom_bar() +
  labs(x = "Metal Related", y = "Count") +
  scale_fill_manual(values = c("FALSE" = "darkred",
                               "TRUE" = "darkgreen")) +
  coord_polar(theta = "y")+
  theme_minimal()
dev.off()


print(paste("number of metal related screens:", sum(phenomics_data$metal_related))) 
print(paste("total number of screens:", length(unique(phenomics_data$name)))) 


print(paste("number of metal related conditions:", length(unique(filter(phenomics_data, metal_related)$conditionset))))
print(paste("total number of conditions:", length(unique(phenomics_data$conditionset)))) 

write.csv(metal_related_screens_data, paste0(outputtables_dir,"/phenomics_metalrelated_BEFORE_manual_curation.csv"),row.names = F)

# Create summary data
metal_related_screens_data <- phenomics_data %>%
                filter(metal_related == T & collection != "hap ?") %>%
                 mutate(medium_clean = ifelse(medium == "YPD", "YPD",
                               ifelse( medium == "SC", "SC",
                                       ifelse(medium == "CSM","SC",
                                              ifelse(medium == "LZM", "LZM",
                                                     ifelse(grepl("YPD [+]",medium)|grepl("YPD[],]",medium),"YPD modified",
                                                            ifelse(grepl("SC [+]", medium)| grepl("SC,", medium) | grepl("SC[(]",medium), "SC modified",
                                                                   ifelse(grepl("SD",medium), "SD",
                                                                          ifelse(grepl("YP",medium),"other YP",
                                                                                 ifelse(grepl("CSM [+]",medium),"SC modified",
                                                                                        ifelse(grepl("YTD",medium), "YTD","other")))))))))))%>%
                group_by( metal) %>%
                mutate(num_screens_for_metal =  length(name)) %>%
                ungroup() %>%
                group_by( metal,medium_clean) %>%
                mutate(num_screens_metalmedium = length(name))%>%
                ungroup()%>%
                group_by(collection, medium_clean, metal) %>%
                mutate(num_screens_metalmediumcollection = length(name))%>%
                ungroup()%>%
                mutate(
                        metal = factor(metal,levels = c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn","Chelator")),
                        medium_clean = factor(medium_clean, levels = c("YPD","YPD modified","SC", "SC modified","YTD","SD", "LZM","other YP","other")))


# Visualize the data
pdf(paste0(plot_dir, "/metal_related_screens_in_phenomics_collection.pdf"),width = 7, height = 5)
ggplot(unique(metal_related_screens_data[,c("metal","medium_clean","num_screens_metalmedium")]), 
       aes(x = metal,
           y = num_screens_metalmedium,
           group = medium_clean,
           fill = medium_clean))+
  geom_bar(stat  = "identity", width = 0.6, linewidth = 0.1, colour = "black")+
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4, "Paired"),
                               RColorBrewer::brewer.pal(8, "Set2")[c(2,3,4,6,8)]) )+
  theme_metallica()
dev.off()

