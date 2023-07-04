
#`---
#`  Title: "summarise_metdepKOgrowthscreen_data.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 18 April 2023
#`  Description: summarises results of the growth screen of y5k KO library on metal depletion media
#`---


#############################################
### source paths functions and libraries  ###
#############################################


# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
#source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

plot_dir = paste0( metdep_KOgrowth_dir, "/output/plots")
dir.create(plot_dir,recursive = T)
output_tables_dir = paste0( metdep_KOgrowth_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)

#################################################################
### Get metdepKOgrowth specific functions and libraries from  ###
#################################################################

source(paste0(code_dir,'/metdep_KOgrowth/metdepKOgrowth_0_libraries_functions.R'))

#############################################
### Read in processed metdepKOgrowth data (output from pyphe ) ###
#############################################

fix_colnames <-  function(x){
  return(unlist(strsplit(x,split="[.]"))[[1]] )
}

metdepKOgrowth_data <- data.table::fread(paste0(metdep_KOgrowth_dir,"/pyphe/pyphe-interpret_reps.csv"),sep = ",",stringsAsFactors = F)

conditions = unlist(lapply(X = colnames(metdepKOgrowth_data), FUN = fix_colnames))

colnames(metdepKOgrowth_data) <- paste(conditions, metdepKOgrowth_data[1,])
metdepKOgrowth_data <- metdepKOgrowth_data[-1,]
metdepKOgrowth_data$ORF <- metdepKOgrowth_data[,1]

metdepKOgrowth_data <-metdepKOgrowth_data[-c(1:2),]%>%
                      reshape2::melt(id.vars = c("ORF"),value.name = "normalised_colony_size")%>%
                      separate(variable, into = c("condition", "replicate"),sep = " ")%>%
                      filter(!condition %in% c("B","Ca_EGTA","EGTA","EDTA","NE") & !ORF %in% c("DNG","CONT","his3_grid","his3_extraGrid"))


metdepKOgrowth_statres <- read.csv(paste0(metdep_KOgrowth_dir,"/pyphe/pyphe-interpret_summaryStats.csv"),stringsAsFactors = F,
                                   row.names = 1)

conditions = unlist(lapply(X = colnames(metdepKOgrowth_statres), FUN = fix_colnames))

colnames(metdepKOgrowth_statres) <- paste(conditions, metdepKOgrowth_statres[1,])

metdepKOgrowth_statres$ORF <- rownames(metdepKOgrowth_statres)

metdepKOgrowth_statres <- metdepKOgrowth_statres[-c(1:2),]%>%
                          reshape2::melt(id.vars = c("ORF"))%>%
                          separate(variable, into = c("condition", "stat_metric"),sep = " ")%>%
                          filter(stat_metric %in% c("mean_effect_size_log2", "mean_fitness_log2","observation_count", "p_Welch_BH","p_Welch") &
                                 !condition %in% c("B","Ca_EGTA","EGTA","EDTA","NE") & !ORF %in% c("DNG","CONT","his3_grid","his3_extraGrid"))%>%
                          reshape2::dcast(ORF+condition ~ stat_metric, value.var = "value")%>%
                          mutate(mean_effect_size_log2 = as.numeric(mean_effect_size_log2),
                                 mean_fitness_log2 = as.numeric(mean_fitness_log2),
                                 observation_count = as.numeric(observation_count),
                                 p_Welch_BH = as.numeric(p_Welch_BH),
                                 significant = ifelse(p_Welch_BH < 0.10 & abs(mean_effect_size_log2) > log2(1.2),T,F ),
                                 metal = gsub("_DiP","", condition),
                                 metal = gsub("5000","", metal),
                                 metal = gsub("1000","", metal),
                                 metal = gsub("500","", metal),
                                 metal = gsub("100","", metal),
                                 metal = gsub("50","", metal),
                                 metal = gsub("20","", metal))%>%
                          group_by(metal, ORF)%>%
                          mutate(sig_in_metal = ifelse(any(significant), T, F))%>%
                          ungroup()%>%
                          group_by(ORF)%>%
                          mutate(sig_in_any_metal = ifelse(any(significant), T, F))%>%
                          ungroup()
                          

write.csv(metdepKOgrowth_statres, paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),row.names =F)
ggplot(metdepKOgrowth_statres,
       aes(x = 2^(mean_fitness_log2),
           y = -log10(p_Welch_BH),
           colour = sig_in_metal))+
  facet_wrap("condition",scales = "free")+
  geom_point()+
  scale_color_manual(values = c("lightgray","maroon"))+
  theme_metallica()


ggplot(metdepKOgrowth_statres,
       aes(x = mean_fitness_log2,
           fill = condition))+
  facet_wrap("condition",scales = "free")+
  geom_density()+
  theme_metallica()
                     

# Create a summary data frame for sig_in_metal
summary_sig_in_metal <- metdepKOgrowth_statres %>%
                        dplyr::select(metal,ORF,sig_in_metal)%>%
                        unique()%>%
                        group_by(metal) %>%
                        summarise(sig_in_metal = sum(sig_in_metal, na.rm=T))%>%
                        filter(metal != "AE")

write.csv(summary_sig_in_metal,paste0(output_tables_dir,"/metdepKOgrowth_numhits_summary.csv"),row.names = F)

# Plot for sig_in_metal
pdf(paste0(plot_dir,"/metdepKOgrowth_numhits_summary.pdf"),width = 7, heigh = 5)
ggplot(summary_sig_in_metal, aes(x=metal, y=sig_in_metal, fill=metal)) + 
  geom_bar(stat='identity', width = 0.6) +
  scale_fill_manual(values = colkey_Ele) +
  theme_metallica() +
  labs(title="Number of ORFs Significant in Each Metal", x="Metal", y="Number of ORFs") 
dev.off()
# Create a summary data frame for sig_in_any_metal
summary_sig_in_any_metal <- metdepKOgrowth_statres %>%
                            dplyr::select(ORF,sig_in_any_metal)%>%
                            unique()%>%
                            summarise(total_hits = sum(sig_in_any_metal, na.rm = T))

print(paste0("number of KOs with significantly altered growth in any metal : ",summary_sig_in_any_metal))

# Generate the list of significant ORFs
sig_ORFs_any_metal <- metdepKOgrowth_statres %>%
                      filter(sig_in_any_metal == T) %>%
                      pull(ORF)

############################
### GSEA with gprofiler2 ###
############################


# set background for enrichments 
background_ORFs <- as.character(unique(metdepKOgrowth_statres$ORF))

# Run the gene set enrichment analysis
gsea_results_any_metal <- gost(query = sig_ORFs_any_metal,
                           organism = "scerevisiae",
                           multi_query = T, 
                           custom_bg = background_ORFs,
                           domain_scope = "custom", 
                           sources = c("GO:MF","GO:BP","KEGG","TF"))

top_terms <- gsea_results_any_metal$result %>%
             filter(p_values < 0.05)%>%
             pull(term_id)

p <- gostplot(gsea_results_any_metal, capped = TRUE, interactive = F)
go_enrichres_plot<- publish_gostplot(p, 
                                     highlight_terms = top_terms,
                                     width = 10, 
                                     height = 15,
                                     filename = paste0(plot_dir,"/metdepKOgrowth_gprofiler_enrich_result_anymetal",
                                                     "_overrep.pdf"))

metals <- unique(metdepKOgrowth_statres$metal)

for (m in 1:length(metals)) {
  
  sig_ORFs_metal <- metdepKOgrowth_statres %>%
    filter(metal == metals[m] & sig_in_metal == T) %>%
    pull(ORF)%>%
    unique()
  
  if(length(sig_ORFs_metal) > 3) {
    
    # Run the gene set enrichment analysis
    gsea_results_metal <- gost(query = sig_ORFs_metal,
                               organism = "scerevisiae",
                               multi_query = T, 
                               custom_bg = background_ORFs,
                               domain_scope = "custom", 
                               sources = c("GO:MF","GO:BP","KEGG","TF"))
    
    if(!is.null(gsea_results_metal)){
      
      # Filter for the top 15 terms
      top_terms <- gsea_results_metal$result %>%
                    arrange(p_values) %>%
                    head(15)%>%
                    pull(term_id)
                  
      # Plot the results with only the top 15 terms
      p <- gostplot(gsea_results_metal, capped = TRUE, interactive = F)
      go_enrichres_plot<- publish_gostplot(p, 
                                           highlight_terms = top_terms,
                                           width = 10, 
                                           height = 8,
                                           filename = paste0(plot_dir,"/metdepKOgrowth_gprofiler_enrich_result_",
                                                             metals[m],"_overrep.pdf"))
   dev.off()
     }
  }
}


############################
### Stouffer Enrichments ###
############################



############################
### Export for metallica ###
############################

## Summarise by metal 


metdepKOgrowth_metalwise_hits<- metdepKOgrowth_data%>%
                                group_by(Metal,ORF)%>%
                                mutate(Significant_Welchs_CondVar_metalwise = ifelse(any(Significant_Welchs_CondVar),T,F),
                                       AbsMax_logrelcolsize = max(abs(log2(Fitness.Relative.to.AE))),
                                       Log2_RelColSize = ifelse(abs(log2(Fitness.Relative.to.AE)) == AbsMax_logrelcolsize,
                                                                log2(Fitness.Relative.to.AE),NA )) %>%
                                ungroup()%>%
                                dplyr::select(Metal,ORF,Ttest_Welchs_CondVar, Significant_Welchs_CondVar_metalwise,Log2_RelColSize)%>%
                                unique()%>%
                                na.omit()

metdepKOgrowth_metalwise_smry <- metdepKOgrowth_metalwise_hits%>%
                                  group_by(Metal)%>%
                                  summarise(num_measured = length(ORF),
                                            num_significant = sum(Significant_Welchs_CondVar_metalwise))%>%
                                  ungroup()

metdepKOgrowth_metalwise_hits <- metdepKOgrowth_metalwise_hits%>%
                                # filter(Significant_Welchs_CondVar_metalwise)%>%
                                 mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)))


write.csv(metdepKOgrowth_metalwise_hits, paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_allhits.csv"),row.names = F)

### write to disk for metallica app


for_metallica_app <- metdepKOgrowth_metalwise_hits%>%
                     mutate(ORF_Gene = paste(ORF,Gene_Name, sep = "|"),
                            label = ifelse(Significant_Welchs_CondVar_metalwise, Gene_Name,NA),
                            colour = ifelse(Significant_Welchs_CondVar_metalwise, Metal, NA))

write.csv(for_metallica_app, paste0(metallica_app_dir,"/metallica_app_metdepKOgrowth.csv"), row.names = F)

pdf(paste0(metdep_KOgrowth_dir,"/output/plots/metdepKOgrowth_numhits_summary.pdf"),width = 8,height = 5)
ggplot()+
  geom_bar(data = metdepKOgrowth_metalwise_smry,
           aes(x = Metal,
               y = num_significant,
               fill = Metal),
           stat = "identity",
           width = 0.6)+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()
dev.off()
