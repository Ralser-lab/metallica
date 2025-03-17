
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
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
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

# how many were cultivated 
metdepKOgrowth_cultivated <- read.csv(paste0(metdep_KOgrowth_dir,"/pyphe/20190122_screenArrangementLayout.csv"),stringsAsFactors = F)%>%
                              filter(!ID %in% c("DNG","CONT","his3_grid","his3_extraGrid",""))%>%
                              pull(ID)%>%
                              unique()

paste("total KOs cultivated:",length(metdepKOgrowth_cultivated))

# All post QC colony sizes 

metdepKOgrowth_data <- read.csv(paste0(metdep_KOgrowth_dir,"/pyphe/pyphe-analyse-report_postQC.csv"), 
                                stringsAsFactors = F)[,c("Colony_size_corr_checked","ID","Condition")]%>%
                       filter(!ID %in% c("DNG","CONT","his3_grid","his3_extraGrid","") &
                               !Condition %in% c("B","Ca_EGTA","EGTA","EDTA","NE"))%>%
                       unique()%>%
                       na.omit()


paste("number of colony sizes:",length(unique(metdepKOgrowth_data$Colony_size_corr_checked)))
paste("number of unique KOs in colony size measurements :",length(unique(metdepKOgrowth_data$ID)))
paste("number of unique conditions :",length(unique(metdepKOgrowth_data$Condition)))  

# Post QC and post stat - statistical results

metdepKOgrowth_statres <- read.csv(paste0(metdep_KOgrowth_dir,"/pyphe/pyphe-interpret_summaryStats.csv"),stringsAsFactors = F,
                                   row.names = 1)


fix_colnames <-  function(x){  return(unlist(strsplit(x,split="[.]"))[[1]] )}

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
                          ungroup()%>%
                          na.omit()
                          
write.csv(metdepKOgrowth_statres, paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),row.names =F)

print(paste("stat res number of unique KOs in stat results :",length(unique(metdepKOgrowth_statres$ORF))))
print(paste("stat res number of combined conditions :",length(unique(metdepKOgrowth_statres$metal))))  

# total number of significant genetic interactions

paste("total number of metal - gene interactions indentified:",
      sum(unique(filter(metdepKOgrowth_statres,condition != "AE")[,c("ORF","metal","sig_in_metal")])$sig_in_metal))
paste0("total unique KOs participating in metal- gene interactions :",
       sum(unique(filter(metdepKOgrowth_statres,condition != "AE")[,c("ORF","sig_in_any_metal")])$sig_in_any_metal))

##########################################################
### GSEa of KO ORFs which have metal-gene interactions ### all metals taken together
##########################################################

# set background for enrichments 
background_ORFs <- as.character(unique(metdepKOgrowth_statres$ORF))

metal_gene_interaction_ORFKOs <- unique(filter(metdepKOgrowth_statres,condition != "AE" & sig_in_any_metal)$ORF)

# Run the gene set enrichment analysis

# overrepresentation
gsea_results_any_metal <- gost(query = metal_gene_interaction_ORFKOs,
                               organism = "scerevisiae",
                               multi_query = F, 
                               custom_bg = background_ORFs,
                               exclude_iea = T,
                               domain_scope = "custom", 
                               sources = c("GO:MF","GO:BP","GO:CC","KEGG","TF"))

top_terms <- gsea_results_any_metal$result %>%
              filter(p_value < 0.01)%>%
              filter( (term_size > 400  & source == "GO:BP") |
                      (term_size < 100 & source == "KEGG") |
                      (term_size > 400 & source == "GO:MF") |
                      (term_size < 1500 | term_size > 200 & source == "KEGG")
                    )%>%
              pull(term_id)

p <- gostplot(gsea_results_any_metal, capped = TRUE, interactive = T)
go_enrichres_plot<- publish_gostplot(p, 
                                     highlight_terms = top_terms,
                                     width = 10, 
                                     height = 15,
                                     filename = paste0(plot_dir,"/metdepKOgrowth_gprofiler_enrich_result_anymetal",
                                                       "_overrep.pdf"))

dev.off()

# underrepresentation
gsea_results_any_metal <- gost(query = metal_gene_interaction_ORFKOs,
                               organism = "scerevisiae",
                               multi_query = F, 
                               custom_bg = background_ORFs,
                               measure_underrepresentation = T,
                               domain_scope = "custom", 
                               sources = c("GO:MF","GO:BP","GO:CC","KEGG","TF"))

top_terms <- gsea_results_any_metal$result %>%
             filter(p_value < 0.05)%>%
             pull(term_id)

p <- gostplot(gsea_results_any_metal, capped = TRUE, interactive = F)
go_enrichres_plot<- publish_gostplot(p, 
                                     highlight_terms = top_terms,
                                     width = 10, 
                                     height = 15,
                                     filename = paste0(plot_dir,"/metdepKOgrowth_gprofiler_enrich_result_anymetal",
                                                       "_underrep.pdf"))

dev.off()


################################
### summarise hits per metal ###
################################

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
  labs(title=" ", x=" ", y="Number of deleted ORFs with growth difference") 
dev.off()
# Create a summary data frame for sig_in_any_metal
summary_sig_in_any_metal <- metdepKOgrowth_statres %>%
                            dplyr::select(ORF,sig_in_any_metal)%>%
                            unique()%>%
                            summarise(total_hits = sum(sig_in_any_metal, na.rm = T))

print(paste0("number of KOs with significantly altered growth in any metal : ",summary_sig_in_any_metal))


############################
### GSEA with gprofiler2 ### - per metal 
############################

metals <- unique(metdepKOgrowth_statres$metal)
GSEA_metdepKOgrowth_metalwise_results <- vector()

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
      top_terms = gsea_results_metal$result %>%
                    arrange(p_values) %>%
                    head(15)%>%
                    pull(term_id)
    
      df2b = data.frame(gsea_results_metal$result[,c("term_name","term_id","source","term_size","query_sizes","intersection_sizes","p_values")])%>%
             mutate(p_values = unlist(p_values),
                    term_size = unlist(term_size),
                    query_sizes = unlist(query_sizes),
                    intersection_sizes = unlist(intersection_sizes)
                    )
      df2b$metal = metals[m]
      GSEA_metdepKOgrowth_metalwise_results <- rbind(GSEA_metdepKOgrowth_metalwise_results,df2b  )  
      # Plot the results with only the top 15 terms
      p = gostplot(gsea_results_metal, capped = TRUE, interactive = F)
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

write.csv(GSEA_metdepKOgrowth_metalwise_results, paste0(output_tables_dir,
                                                        "/metdepKOgrowth_metalwise_enrichments_GSA_gprofiler2.csv"),
          row.names = F)

############################
### Export for metallica ###
############################

ORF2GeneName_map <- data.frame(ORF = unique(metdepKOgrowth_statres$ORF),
                               Gene_Name = as.character(lapply(unique(metdepKOgrowth_statres$ORF),
                                                               convert_ORF2SingleGeneName)))
## Summarise by metal 
  
for_metallica_app <- merge(unique(metdepKOgrowth_statres[,c("ORF","condition","p_Welch_BH","sig_in_metal","metal")]),
                           metdepKOgrowth_data,
                                  by.x = c("ORF","condition"),
                                  by.y = c("ID","Condition"))

for_metallica_app <- merge(for_metallica_app,ORF2GeneName_map, by = "ORF")%>%
                     mutate(GORF_Gene = paste(ORF,Gene_Name, sep = "|"),
                            label = ifelse(sig_in_metal, Gene_Name,NA),
                            colour = ifelse(sig_in_metal, metal, NA))

write.csv(for_metallica_app, paste0(metallica_app_dir,"/data/metallica_app_metdepKOgrowth.csv"), row.names = F)

