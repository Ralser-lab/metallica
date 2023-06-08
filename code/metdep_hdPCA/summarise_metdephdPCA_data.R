
#`---
#`  Title: "summarise_metdephdPCA_data.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 20 April 2023
#`  Description: summarises results of the homodimer protein complementation assay on metal depletion media


#############################################
### source paths functions and libraries  ###
#############################################


# general
source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
#source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

#################################################################
### Get metdepKOgrowth specific functions and libraries from  ###
#################################################################

source(paste0(code_dir,'/metdep_hdPCA/metdephdPCA_0_libraries_functions.R'))


dir.create(paste0(metdep_hdPCA_dir,"/output/plots/"),recursive = T)
dir.create(paste0(metdep_hdPCA_dir,"/output/tables/"),recursive = T)

#############################################
### Read in processed metdepKOgrowth data ###
#############################################

metdep_hdPCA_data <- read.csv(paste0(metdep_hdPCA_dir,"/mainscreen/output/tables/hdPCA CyclicLoessNorm colony sizes w Statistical Test Results w AE.csv"),
                              stringsAsFactors = F)%>%
  dplyr::select(Condition,ID.Condition, ORF, Rel.Mean.Log2.CS ,
                Ttest_Welchs, Ttest_Welchs.FDR)%>%
  na.omit()%>%
  filter(!Condition %in% c("B","Ca_EGTA","EGTA","EDTA"))%>%
  mutate(Metal = gsub("_DiP","", Condition),
         Metal = gsub("5000","", Metal),
         Metal = gsub("1000","", Metal),
         Metal = gsub("500","", Metal),
         Metal = gsub("100","", Metal),
         Metal = gsub("20","", Metal),
         
         Significant_WelchsFDR = ifelse(Ttest_Welchs < 0.01 & abs(Rel.Mean.Log2.CS) > log2(1.2) ,T,F))

## Plot the pvalue vs fitness relative to AE

pdf(paste0(metdep_hdPCA_dir,"/output/plots/metdephdPCA_volcanoplots.pdf"),width = 15, height = 14)
ggplot(metdep_hdPCA_data,
       aes(x = Rel.Mean.Log2.CS,
           y = -log10(Ttest_Welchs.FDR),
           colour = Significant_WelchsFDR))+
  geom_point(alpha = 0.5)+
  facet_wrap("Condition")+
  theme_metallica()+
  scale_colour_manual(values = c("lightgray","darkred"))+
  theme(legend.position = "none")
dev.off()

## Summarise by metal 


metdephdPCA_metalwise_hits <- metdep_hdPCA_data%>%
  group_by(Metal,ORF)%>%
  mutate(Significant_WelchsFDR_metalwise = ifelse(any(Significant_WelchsFDR),T,F),
         AbsMax_log2relmeancolsize = abs(max(Rel.Mean.Log2.CS)),
         Max_log2relmeancolsize = ifelse(abs(Rel.Mean.Log2.CS) == AbsMax_log2relmeancolsize, Rel.Mean.Log2.CS, NA) )%>%
  ungroup()%>%
  dplyr::select(Metal,ORF, Significant_WelchsFDR_metalwise,Max_log2relmeancolsize)%>%
  unique()%>%
  na.omit()

metdephdPCA_metalwise_smry <- metdephdPCA_metalwise_hits%>%
  group_by(Metal)%>%
  summarise(num_measured = length(ORF),
            num_significant = sum(Significant_WelchsFDR_metalwise))%>%
  ungroup()

metdephdPCA_metalwise_hits <- metdephdPCA_metalwise_hits%>%
  # filter(Significant_WelchsFDR_metalwise)%>%
  mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)))


write.csv(metdephdPCA_metalwise_hits, paste0(metdep_hdPCA_dir,"/output/tables/metdephdPCA_allhits.csv"),
          row.names = F)

### write to disk for metallica app


for_metallica_app <- metdephdPCA_metalwise_hits %>%
  mutate(label = ifelse(Significant_WelchsFDR_metalwise, Gene_Name,NA),
         colour = ifelse(Significant_WelchsFDR_metalwise, Metal, NA))

write.csv(for_metallica_app, paste0(metallica_app_dir,"/metallica_app_metdephdPCA.csv"), row.names = F)

pdf(paste0(metdep_hdPCA_dir,"/output/plots/metdephdPCA_numhits_summary.pdf"),width = 10,height = 6)
ggplot()+
  geom_bar(data = metdephdPCA_metalwise_smry,
           aes(x = Metal,
               y = num_significant,
               fill = Metal),
           stat = "identity")+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()
dev.off()

############################
### GSEA using gprofiler ### 
############################

eles <- unique(metdephdPCA_metalwise_hits$Metal)

for(e in 1:length(eles)){
  
  
  genes_for_gprofiler = unique(filter(metdephdPCA_metalwise_hits,Significant_WelchsFDR_metalwise == T & Metal == eles[e])$ORF)
  
  if(length(genes_for_gprofiler) > 3) {
    
    background_genes <- as.character(unique( filter(metdep_hdPCA_data,Metal == eles[e])$ORF))
    
    ## overrepresentation
    go_enrich_res <- gost(query = genes_for_gprofiler,organism = "scerevisiae", 
                          multi_query = T, custom_bg = background_genes,
                          domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))
    if(!is.null(go_enrich_res)){
      
      #   go_enrichres_plot <- gostplot(go_enrich_res, capped = TRUE, interactive = F)
      
      go_enrich_res_df <- go_enrich_res$result
      
      # publish_gostplot(go_enrichres_plot, highlight_terms = go_enrich_res_df$term_id, 
      #                 width = 10, height = ifelse(eles[e] %in% c("Mg","K","Zn"),20, 8), 
      #                filename = paste0(metdep_KOgrowth_dir,"/output/plots/metdepKOgrowth_gprofiler_enrich_result_",eles[e] ,"_overrep.pdf") )
    }
  }
}
