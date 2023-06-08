#####################################################################
####          Script to plot results of ensemble clustering      ####
#####################################################################

#`---
#`  Title: "metperWTproteomics_summarise_ensembleclustering.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 June 2021 
#`  Description: Script to plot results of ensemble clustering by Oliver Lemke (in python)
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))

# specific

source(paste0(code_dir,"/metpertWTproteomics/metpertWTproteomics_0_libraries_functions.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/ensemble_clustering")
dir.create(plots_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering")
dir.create(output_tables_dir,recursive = T)

####################################################################
### Read in Ensemble Clustering res using full proteome profiles ###
####################################################################

full_EC <- read.csv(paste0(metpert_WTproteomics_dir,"/outputs/ensemble_clustering/Clustered_full_solo.csv"),stringsAsFactors = F)
full_EC <- merge(full_EC, GenProt_SGD, by="ORF")

clus_annostat <- full_EC %>%
  group_by(Cluster, Uniprot.Annotation.Score)%>%
  summarize(Num_ORFs = length(ORF))%>%
  ungroup()

pubmedanno_AllLit <- read.csv(paste0(projdir,"/Databases/Number of pubmed publications per ORF All Literature.csv"),
                              stringsAsFactors = F)%>%
  filter(!grepl("ARS",ORF))
colnames(pubmedanno_AllLit)[3] <- "Num_Public_AllLit"

pubmedanno_PrimaryLit <- read.csv(paste0(projdir,"/Databases/Number of pubmed publications per ORF Primary Literature.csv"),
                                  stringsAsFactors = F)%>%
  filter(!grepl("ARS",ORF))
colnames(pubmedanno_PrimaryLit)[3] <- "Num_Public_PrimaryLit"

quantiles_AL <- quantile(pubmedanno_AllLit$Num_Public_AllLit)

clus_ALanno <- merge(full_EC, pubmedanno_AllLit,by="ORF")%>%
  mutate(Quantile = ifelse( Num_Public_AllLit < quantiles_AL[[2]], paste("<",quantiles_AL[[2]]),
                            ifelse(Num_Public_AllLit < quantiles_AL[[3]] , paste(quantiles_AL[[2]],"-",quantiles_AL[[3]]),
                                   ifelse(Num_Public_AllLit < quantiles_AL[[4]],paste(quantiles_AL[[3]] ,"-",quantiles_AL[[4]]),
                                          paste(quantiles_AL[[4]],"-",quantiles_AL[[5]]))))
  )
write.csv(clus_ALanno,"ORF ensemble clustering lit anno.csv")

clus_ALanno <- clus_ALanno%>%
  group_by(Cluster, Quantile)%>%
  mutate(Num_ORFs = length(ORF))%>%
  ungroup()

clus_ALanno$Quantile = factor(clus_ALanno$Quantile,
                              levels = c("< 8","8 - 26","26 - 65","65 - 1212"))


quantiles_PL <- quantile(pubmedanno_PrimaryLit$Num_Public_PrimaryLit)

clus_PLanno <- merge(full_EC, pubmedanno_PrimaryLit,by="ORF")%>%
  mutate(Quantile = ifelse( Num_Public_PrimaryLit < quantiles_PL[[2]], paste("<",quantiles_PL[[2]]),
                            ifelse(Num_Public_PrimaryLit < quantiles_PL[[3]] , paste(quantiles_PL[[2]],"-",quantiles_PL[[3]]),
                                   ifelse(Num_Public_PrimaryLit < quantiles_PL[[4]],paste(quantiles_PL[[3]] ,"-",quantiles_PL[[4]]),
                                          paste(quantiles_PL[[4]],"-",quantiles_PL[[5]]))))
  )%>%
  group_by(Cluster, Quantile)%>%
  mutate(Num_ORFs = length(ORF))%>%
  ungroup()

clus_PLanno$Quantile = factor(clus_PLanno$Quantile,
                              levels = c("< 3","3 - 8","8 - 20","20 - 435"))

write.csv(clus_PLanno,"ORF ensemble clustering primarylit anno.csv")

pdf(paste0(ensclust_dir,"/Poorly characterised proteins in clusters.pdf"),width=10,height=5)
ggplot(clus_annostat,
       aes(x=Cluster,
           y=Num_ORFs,
           fill=factor(Uniprot.Annotation.Score)))+
  geom_bar(stat="identity", position = "stack", colour="white", width = 0.6)+
  scale_fill_viridis_d(begin =0.2,end=0.9, option = "B")+
  theme_SKA()+
  theme(legend.position = "bottom")+
  labs(fill = "UniProt Annotation Score", y = "Number of proteins")

ggplot(filter(clus_annostat,Uniprot.Annotation.Score <=2),
       aes(x=Cluster,
           y=Num_ORFs,
           fill=factor(Uniprot.Annotation.Score)))+
  geom_bar(stat="identity", position = "stack", colour="white", width = 0.6)+
  scale_fill_viridis_d(begin =0.2,end=0.47, option = "B")+
  theme_SKA()+
  theme(legend.position = "bottom")+
  labs(fill = "UniProt Annotation Score", y = "Number of proteins")

## All Lit Quantiles


ggplot(unique(clus_ALanno[,c("Cluster","Num_ORFs","Quantile")]),
       aes(x=Cluster,
           y=Num_ORFs,
           fill=Quantile))+
  geom_bar(stat="identity", position = "stack", colour="black",
           size=0.25, width = 0.5)+
  scale_fill_brewer(palette="BuPu",direction=-1)+
  theme_SKA()+
  theme(legend.position = "bottom")

ggplot(unique(clus_PLanno[,c("Cluster","Num_ORFs","Quantile")]),
       aes(x=Cluster,
           y=Num_ORFs,
           fill=Quantile))+
  geom_bar(stat="identity", position = "stack", colour="black",
           size=0.25)+
  scale_fill_brewer(palette="PuRd",direction=-1)+
  theme_SKA()+
  theme(legend.position = "bottom")


dev.off()


#######################################################
### Read in protein quantities to plot cluster-wise ###
#######################################################

#PQ_wNA <- fread(paste0(FinalOutCSVdir,"/Protein Quantities w NAs FullMatrix.csv"),stringsAsFactors = F)%>%
#      group_by(Protein.Ids,BioSpecID,Element,Element.Concentration)%>%
#  summarize(Mean_Log2_FCvsAE = mean(Log2_FCvsAE,na.rm=T))%>%
#     mutate(GeneName = as.character(lapply(Protein.Ids,convert_Uniprot2SingleGeneName)))


PQs <- fread(paste0(FinalOutCSVdir,("/AllData matrix log2FCvsAE GeneName.csv")),stringsAsFactors = F)%>%
  melt(id.vars="V1")
colnames(PQs) <- c("Gene","BioSpecID","Log2_FCvsAE")

get_ele<-function(e){unlist(strsplit(e," "))[[1]]}
get_eleconc<-function(e){unlist(strsplit(e," "))[[2]]}

PQs <- PQs %>%
  data.frame()%>%
  mutate(BioSpecID = as.character(BioSpecID),   
         Element = as.character(lapply(BioSpecID,get_ele)),
         Element_Conc = factor(as.character(lapply(BioSpecID,get_eleconc)),levels = c("0","0.01","0.02","0.05","0.1","0.2","0.5","2","5","10","20","50","100")))%>%
  filter(Element_Conc != "0")%>%
  group_by(BioSpecID,Gene)%>%
  mutate(Mean_Log2_FCvsAE = mean(Log2_FCvsAE,na.rm=T))%>%
  ungroup()


full_EC <- full_EC %>%
  mutate(GeneName = as.character(lapply(ORF,convert_ORF2SingleGeneName)))

dir.create(paste0(ensclust_dir,"/Full/MeanProfiles"),recursive = T)
setwd(paste0(ensclust_dir,"/Full/MeanProfiles"))

go_enrich_res_fullclustering_df <- vector()

for(i in 1:length(unique(full_EC$Cluster))){
  
  Clus_df <- filter(PQs, Gene %in% filter(full_EC,Cluster == i)$GeneName)%>%
    group_by(BioSpecID,Element,Element_Conc)%>%
    summarize(Mean_FCvsAE_Cluster = mean(Log2_FCvsAE,na.rm=T))%>%
    ungroup()
  
  Clus_df_PQs <- filter(PQs, Gene %in% filter(full_EC,Cluster == i)$GeneName)
  
  pdf(paste0("Mean log2FCvsAE of all proteins in cluster",i,".pdf"),width=10,height=5)
  print( ggplot(Clus_df,
                aes(x=Element_Conc,
                    y=Mean_FCvsAE_Cluster,
                    colour = Element,
                    group=Element))+
           geom_point(size=3)+
           geom_line(size=1,alpha=0.8)+
           scale_colour_manual(values = colkey_Ele)+
           theme_SKA()+
           labs(title = paste0("Cluster",i))
         
  )
  
  print( ggplot(Clus_df_PQs,
                aes(x=Element_Conc,
                    y=Mean_Log2_FCvsAE,
                    colour = Element,
                    group=Element))+
           geom_smooth(alpha = 0.1,method = "loess")+
           scale_colour_manual(values = colkey_Ele)+
           theme_SKA()+
           labs(title = paste0("Cluster ",i))
  )
  dev.off()
  
  Clus_df <- filter(PQs, Gene %in% filter(full_EC,Cluster == i)$GeneName)
  nump = length(unique(Clus_df$Gene))
  
  pdf(paste0("Mean log2FCvsAE per protein cluster",i,".pdf"),width = round(sqrt(nump),0)*5.6,height =round(sqrt(nump),0)*4)
  print(
    ggplot(Clus_df,
           aes(x=log2(as.numeric(as.character(Element_Conc))),
               y=Mean_Log2_FCvsAE,
               colour = Element,
               group=Element))+
      geom_point(size=3)+
      geom_line(size=1,alpha=0.7)+
      facet_wrap("Gene",ncol=round(sqrt(nump),0),scales="free")+
      scale_colour_manual(values = colkey_Ele)+
      theme_SKA()
  )
  
  # geom smooth
  print(
    ggplot(Clus_df_PQs,
           aes(x=log2(as.numeric(as.character(Element_Conc))),
               y=Mean_Log2_FCvsAE,
               colour = Element,
               group=Element))+
      geom_smooth(alpha = 0.1,method = "loess")+
      facet_wrap("Gene",ncol=round(sqrt(nump),0),scales="free")+
      scale_colour_manual(values = colkey_Ele)+
      theme_SKA()
  )
  dev.off()
  
  ## GO enrichment
  
  background_proteome <- unique(full_EC$GeneName)
  
  ## overrepresentation
  go_enrich_res <- gost(query = unique(Clus_df_PQs$Gene),organism = "scerevisiae", 
                        multi_query = TRUE, custom_bg = background_proteome,
                        domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))
  
  go_enrich_res_df <- go_enrich_res$result
  if(!is.null(go_enrich_res_df)){
    go_enrich_res_df$Cluster = i
    go_enrich_res_fullclustering_df <-rbind(go_enrich_res_fullclustering_df,go_enrich_res_df)
  }
}

## Plot Enriched GO terms in all clusters

go_enrich_res_fullclustering_df_toplot <- go_enrich_res_fullclustering_df%>%
  dplyr::select(source, term_name,p_values, term_size, Cluster)%>%
  unique()%>%
  mutate(p_values = as.numeric(as.character(p_values)))

write.csv(go_enrich_res_fullclustering_df_toplot,paste0(ensclust_dir,"/gprofiler_enrichments_per_cluster.csv"))

#################################
### Summarise unchar proteins ###
#################################

adjpvthresh = 0.05
fcthresh = 1.5

DA_prots <- read.csv(paste0(FinalOutCSVdir,"/lmfit DE res with PQ and SigNotSig AdjPVthresh ",
                            adjpvthresh," FCthresh ",fcthresh,".csv"))%>%
  data.frame()%>%
  filter(Significant == 1)%>%
  dplyr::select(Genes, Element, ORF, Significant)%>%
  unique()

DA_prots_poorly_char <- merge(DA_prots,GenProt_SGD, by = "ORF")%>%
  dplyr::select(ORF,Genes, Element, Uniprot.ID, Uniprot.Annotation.Score)%>%
  filter( Uniprot.Annotation.Score <2)

full_EC <- full_EC %>%
  mutate(DiffAbund = ifelse(ORF %in% unique(DA_prots$ORF), T, F))


full_EC_poorly_char <- full_EC%>%
  filter(Uniprot.Annotation.Score  <=2)

write.csv(full_EC_poorly_char,
          paste0(ensclust_dir,"/poorly characterised proteins in clusters.csv"))


###############################
### ELEMENT WISE CLUSTERING ###
###############################

sig_eles <- c("Ca","Cu","Fe","Zn")
sig_genes_df <- read.csv(paste0(FinalOutCSVdir,"/lmfit DE res with PQ and SigNotSig AdjPVthresh ",
                                adjpvthresh," FCthresh ",fcthresh,".csv"), stringsAsFactors = F)%>%
  filter(Significant == 1)%>%
  dplyr::select(Genes,Element)%>%
  unique()


ele_wise_cluster <- vector()
ele_wise_cluster_enrichment <- vector()

for(e in 1:length(sig_eles)){
  
  
  ele_clust_dir <- paste0(ensclust_dir,"/Elements/",sig_eles[e],"/Output")
  
  num_ele_clust <- length(dir(ele_clust_dir)[grep(".txt", dir(ele_clust_dir))])
  
  
  for(c in 1:num_ele_clust){
    
    cl_uniprots = read.table(paste0(ele_clust_dir,"/C",c,".txt"), stringsAsFactors = F)
    cl_uniprots$Cluster = c
    cl_uniprots$Element_Series = sig_eles[e]
    colnames(cl_uniprots)[[1]] = "Uniprot_ID"
    
    ele_wise_cluster <- rbind(ele_wise_cluster,cl_uniprots)
    
    
    background_proteome <- filter(sig_genes_df,Element == sig_eles[e] & Genes %in% unique(filter(PQs, Element == sig_eles[e])$Gene))$Genes
    
    ## overrepresentation
    go_enrich_res <- gost(query = as.character(lapply(cl_uniprots$Uniprot_ID,convert_Uniprot2SingleGeneName)),organism = "scerevisiae", 
                          multi_query = TRUE, custom_bg = background_proteome,
                          domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))
    
    go_enrich_res_df <- go_enrich_res$result
    if(!is.null(go_enrich_res_df)){
      go_enrich_res_df$Cluster = c
      go_enrich_res_df$Element = sig_eles[e]
      ele_wise_cluster_enrichment <-rbind(ele_wise_cluster_enrichment,go_enrich_res_df)
    }
    
  }
  
}


ele_wise_cluster <- ele_wise_cluster %>%
  mutate(ORF = as.character(lapply(Uniprot_ID, convert_Uniprot2singleORF)),
         Gene = as.character(lapply(Uniprot_ID,convert_Uniprot2SingleGeneName)))

ele_wise_cluster<-  merge(ele_wise_cluster,GenProt_SGD[,c("Uniprot.Annotation.Score","ORF")], by = "ORF")

Conc_Conversion <- read.csv(paste0(FinalOutCSVdir,"/lmfit DE res with PQ and SigNotSig AdjPVthresh ",
                                   adjpvthresh," FCthresh ",fcthresh,".csv"))%>%
  dplyr::select(BioSpecID,Element,Element.Concentration)%>%
  unique()

PQ_df_elewise_clusters <- PQs%>% ## USING IMPUTED VALUE DATAFRAME FOR PLOTTING RESULTS OF CLUSTERING
  filter(Element %in% sig_eles)%>%
  mutate(ORF = as.character(lapply(Gene, convert_GeneName2ORF)))%>%
  dplyr::select(Gene, ORF, Element, BioSpecID, Log2_FCvsAE)%>%
  unique()

PQ_df_elewise_clusters <- merge(PQ_df_elewise_clusters, Conc_Conversion, by = c("BioSpecID","Element"))


PQ_df_elewise_clusters <- merge(PQ_df_elewise_clusters, ele_wise_cluster[,c("ORF","Cluster","Element_Series","Uniprot.Annotation.Score")], by = "ORF")%>%
  filter(Element == Element_Series)


PQ_df_elewise_clusters <-PQ_df_elewise_clusters%>%
  group_by(Element,Gene)%>%
  mutate(label = ifelse(Element.Concentration == min(Element.Concentration),
                        ifelse(Uniprot.Annotation.Score < 3,Gene,NA),NA),
         label = ifelse(Element == "Zn" & !is.na(label),
                        ifelse(BioSpecID == "Zn 0.01",label,NA),label))%>%
  ungroup()



## Plot all element wise clusters 

pdf(paste0(ensclust_dir,"/Elements/all_elementwise_clusters_noscalelim.pdf"),width = 30,height = 39)
ggplot(PQ_df_elewise_clusters,
       aes(x = log2(as.numeric(as.character(Element.Concentration))),
           y = Log2_FCvsAE   ))+
  facet_wrap(c("Element_Series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.25,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+
  geom_hline(yintercept = 0, size = 0.1)+
  geom_vline(xintercept = 0, size = 0.1)+
  geom_text_repel(aes(label = label, colour = Uniprot.Annotation.Score),size = 3)+
  theme_SKA()
dev.off()



pdf(paste0(ensclust_dir,"/Elements/all_elementwise_clusters_scalelim.pdf"),width = 30,height = 39)
ggplot(PQ_df_elewise_clusters,
       aes(x = log2(as.numeric(as.character(Element.Concentration))),
           y = Log2_FCvsAE   ))+
  facet_wrap(c("Element_Series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.5,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+
  geom_hline(yintercept = 0, size = 0.2)+
  geom_vline(xintercept = 0, size = 0.2)+
  ylim(-2,2)+
  geom_text_repel(aes(label = label, colour = Uniprot.Annotation.Score),size = 2.5)+
  theme_SKA()+
  labs(x = "log2(metal concentration)", y = "log2(fold difference vs allele)")
dev.off()


#####################################
### Summarise metal wise clusters ###
#####################################

PQ_df_elewise_clusters_summary <- PQ_df_elewise_clusters %>%
  group_by(Element, Cluster, Uniprot.Annotation.Score)%>%
  dplyr::summarise(num_protein = length(Gene))

pdf(paste0(ensclust_dir,"/Elements/metal_wise_clustering_summary.pdf"),width = 16, height = 3)
ggplot(PQ_df_elewise_clusters_summary,
       aes(x = Cluster, 
           y = num_protein,
           fill = Uniprot.Annotation.Score))+
  geom_bar(stat = "identity", position = "stack", size = 0, width = 0.6)+
  facet_wrap("Element", ncol = 4)+
  scale_fill_viridis_d(begin =0.2,end=0.9, option = "B")+
  theme_SKA()
dev.off()

PQ_df_elewise_clusters_summary_unchar <- PQ_df_elewise_clusters %>%
  filter(Uniprot.Annotation.Score < 3)%>%
  dplyr::select(Gene,Element, Cluster, Uniprot.Annotation.Score)%>%
  unique()

pdf(paste0(ensclust_dir,"/Elements/poorly_characterised_proteins_metal_wise_cluster_membership.pdf"),width = 15, height = 20)
ggplot(PQ_df_elewise_clusters_summary_unchar,
       aes(x = Gene,
           y = factor(Cluster),
           fill = Element,
           colour = Element))+
  geom_tile()+
  scale_fill_manual(values = colkey_Ele)+
  scale_colour_manual(values = colkey_Ele)+
  coord_flip()+
  theme_SKA()+
  labs("cluster number")
dev.off()

#########################################################
### Merge dataset with GO MF metal and MetalPDB annos ###
#########################################################

GO_gset_MF_metal_binding <- GO_gset_MF %>%
  filter(grepl("binding",term))%>%
  filter(grepl("calcium",term) |
           grepl("copper",term) |
           grepl("iron",term) |
           grepl("magnesium",term) |
           grepl("manganese",term) |
           grepl("molybdenum", term) |
           grepl("potassium", term) |
           grepl("sodium", term) |
           grepl("zinc", term) |
           grepl("metal",term))%>%
  filter(!grepl("calcium-dependent",term))%>%
  pull(ORF)%>%
  unique()

GO_gset_MF_binding_metalwise <- GO_gset_MF %>%
  filter(grepl("binding",term))%>%
  filter(grepl("calcium",term) |
           grepl("copper",term) |
           grepl("iron",term) |
           grepl("magnesium",term) |
           grepl("manganese",term) |
           grepl("molybdenum", term) |
           grepl("potassium", term) |
           grepl("sodium", term) |
           grepl("zinc", term) )%>%
  filter(!grepl("calcium-dependent",term))%>%
  mutate(term = ifelse(grepl("calcium",term),"Ca",
                       ifelse(grepl("copper", term),"Cu",
                              ifelse(grepl("iron",term),"Fe",
                                     ifelse(grepl("magnesium",term),"Mg",
                                            ifelse(grepl("manganese",term),"Mn",
                                                   ifelse(grepl("molybdenum",term),"Mo",
                                                          ifelse(grepl("potassium",term),"K",
                                                                 ifelse(grepl("sodium",term),"Na",
                                                                        ifelse(grepl("zinc",term),"Zn",NA))))))))))%>%
  dplyr::select(ORF,term)%>%
  unique()

MetPDB_biosig <- filter(MetPDB,term %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  pull(ORF)%>%
  unique()

ele_wise_clusters_wanno <- ele_wise_cluster %>%
  mutate(GO_MetPDB_metal_binding = ifelse(ORF %in% c(GO_gset_MF_metal_binding,MetPDB_biosig),T,F))


met_specific_metal_binders_meas <- unique(rbind(GO_gset_MF_binding_metalwise,MetPDB_biosig))

ele_wise_clusters_wanno <- merge(ele_wise_clusters_wanno, met_specific_metal_binders_meas,by = "ORF",all.x = "T")
colnames(ele_wise_clusters_wanno)[[ncol(ele_wise_clusters_wanno)]] <- "metal_specific_metalbinding"

## Plot summary of Metal related annotations ###

ele_wise_clusters_wanno_metalbinding_summary <- ele_wise_clusters_wanno%>%
  group_by(Element_Series, Cluster, metal_specific_metalbinding)%>%
  dplyr::summarise(num_protein = length(Gene))

pdf(paste0(ensclust_dir,"/Elements/metal_wise_clustering_metbinders_summary.pdf"),width = 16, height = 3)

ggplot(ele_wise_clusters_wanno_metalbinding_summary,
       aes(x = Cluster, 
           y = num_protein,
           fill = metal_specific_metalbinding))+
  geom_bar(stat = "identity", position = "stack", size = 0, width = 0.6)+
  facet_wrap("Element_Series", ncol = 4)+
  scale_fill_manual(values = colkey_Ele)+
  theme_SKA()+
  labs(y = "metal binding annotations in each cluster")

ele_wise_clusters_anymetbinding_summary <- ele_wise_clusters_wanno%>%
  group_by(Element_Series, Cluster, GO_MetPDB_metal_binding)%>%
  dplyr::summarise(num_protein = length(Gene))

ggplot(ele_wise_clusters_anymetbinding_summary,
       aes(x = Cluster, 
           y = num_protein,
           fill = GO_MetPDB_metal_binding))+
  geom_bar(stat = "identity", position = "stack", size = 0, width = 0.6, alpha = 0.7)+
  facet_wrap("Element_Series", ncol = 4)+
  scale_fill_manual(values = c("darkred","darkgreen"))+
  theme_SKA()+
  labs(y = "metal binding annotations in each cluster")
dev.off()

######################################
### Merge with cluster enrichments ###
######################################

ele_wise_clusters_wanno<- merge(ele_wise_clusters_wanno ,unique(ele_wise_cluster_enrichment[,c("source","term_name","Cluster","Element")]), 
                                all.x= T,
                                by.x = c("Element_Series","Cluster"), by.y = c("Element","Cluster"))

write.csv(ele_wise_clusters_wanno,paste0(ensclust_dir,"/Elements/all_element_wise_clusters_with_anno.csv"), row.names = F)

# only unchar proteions

write.csv(filter(ele_wise_clusters_wanno,Uniprot.Annotation.Score < 3),paste0(ensclust_dir,"/Elements/all_element_wise_clusters_with_anno_poorlycharprots.csv"))


#############################################
### Combine Full and ele-wise clustering ####
#############################################

colnames(ele_wise_clusters_wanno) <- c("Element","Elewise_Cluster","ORF","Uniprot_ID","Gene","Uniprot.Annotation.Score","GO_MetPDB_metal_binding",
                                       "metal_specific_metalbinding","Enrichment_database","Elewise_Clustering_Enrichment")

full_EC <- merge(full_EC,go_enrich_res_fullclustering_df_toplot[,c("source","term_name","Cluster")], by = "Cluster", all.x=T)

full_EC <- full_EC[,c("Cluster","ORF","source","term_name")]

colnames(full_EC) <- c("Full_Cluster","ORF","Enrichment_database","Ful_Clustering_Enrichment")

ele_wise_clusters_wanno_wFC <- merge(ele_wise_clusters_wanno, full_EC, by = "ORF", all = T)

write.csv(ele_wise_clusters_wanno_wFC,paste0(ensclust_dir,"/elewise_full_combined_clusters_wanno.csv"))
write.csv(filter(ele_wise_clusters_wanno_wFC,Uniprot.Annotation.Score < 3),paste0(ensclust_dir,"/elewise_full_combined_clusters_wanno_porly_charprots.csv"))

