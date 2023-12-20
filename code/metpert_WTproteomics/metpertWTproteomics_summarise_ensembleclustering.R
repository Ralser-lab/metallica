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

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/ensemble_clustering")
dir.create(paste0(plot_dir,"/allmetal/clusterwise/"), recursive = T)
dir.create(paste0(plot_dir,"/metalwise/clusterwise/"), recursive = T)

input_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering/clustering_results")

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering/clustering_results/summary")
dir.create(output_tables_dir,recursive = T)

####################################################################
### Read in Ensemble Clustering res using full proteome profiles ###
####################################################################

full_EC <- read.csv(paste0(input_tables_dir,"/allmetals/Full/Output_data/Clustered_full_solo.csv"),stringsAsFactors = F)
full_EC_SGD <- merge(full_EC, na.omit(unique(GenProt_SGD[,c("ORF","Uniprot.ID","Uniprot.Annotation.Score")])), 
                 by.x = "ORF",
                 by.y = "Uniprot.ID")

colnames(full_EC_SGD) <- c("Uniprot.ID","Cluster","ORF","Uniprot.Annotation.Score")

clus_annostat <- full_EC_SGD %>%
                 group_by(Cluster, Uniprot.Annotation.Score)%>%
                 summarize(Num_ORFs = length(ORF))%>%
                 ungroup()

## num poorly characterised proteins in full EV

poorly_char_prots_fullEC <- unique(filter(full_EC_SGD, Uniprot.Annotation.Score < 3)$ORF)
length(poorly_char_prots_fullEC)

pdf(paste0(plot_dir,"/allmetal/poorly_characterised_proteins_in_clusters_allmetalclustering.pdf"),width=5,height=5)
ggplot(clus_annostat,
       aes(x=Cluster,
           y=Num_ORFs,
           fill=factor(Uniprot.Annotation.Score)))+
  geom_bar(stat="identity", position = "stack", colour=NA, width = 0.6, size = 0, linewidth = 0.1)+
  scale_fill_viridis_d(begin =0.2,end=0.9, option = "B")+
  theme_metallica()+
  theme(legend.position = "bottom")+
  labs(fill = "UniProt Annotation Score", y = "Number of proteins")

ggplot(filter(clus_annostat,Uniprot.Annotation.Score <=2),
       aes(x=Cluster,
           y=Num_ORFs,
           fill=factor(Uniprot.Annotation.Score)))+
  geom_bar(stat="identity", position = "stack", colour=NA, width = 0.6, size = 0, linewidth = 0.1)+
  scale_fill_viridis_d(begin =0.2,end=0.47, option = "B")+
  theme_metallica()+
  theme(legend.position = "bottom")+
  labs(fill = "UniProt Annotation Score", y = "Number of proteins")

dev.off()

#############################################
### orphan metal binders and transporters ###
#############################################

full_EC_ORFnames <- full_EC%>%
                    mutate(ORF = as.character(lapply(ORF,convert_Uniprot2singleORF)))

## metal binders 

full_EC_orphan_binders <- merge(full_EC_ORFnames,rbind(GO_gset_MF_metalbinding_unspecific,GO_gset_MF_binding_metalwise), by = "ORF", all.x = T)%>%
                          mutate(term = ifelse(!is.na(term),term, "none"),
                                 color = ifelse(term == "none","beige",
                                                ifelse(term == "unspecific","black",
                                                       colkey_Ele[term])))


orphan_metb <- ggplot(full_EC_orphan_binders,
               aes(x=Cluster,
                   fill=color))+
                geom_bar(stat="count",
                         position = "stack", colour="black", 
                         width = 0.6, size = 0, linewidth = 0.1)+
                theme_metallica()+
                theme(legend.position = "bottom")+
                scale_fill_identity()+
                labs(fill = "metal binding annotation",
                     y = "Number of metal binders")


pdf(paste0(plot_dir,"/allmetal/orphan_metal_binders_transporter_per_cluster.pdf"),width = 7,height =5)
orphan_metb
dev.off()

#######################################################
### Read in protein quantities to plot cluster-wise ###
#######################################################

FCvsAE_allproteins_env<- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
                       dplyr::select(Genes, ORF,Element, BioSpecID, Element.Concentration,Mean_of_reps_Log2FC_vs_AE)%>%
                       unique()


FCvsAE_allproteins_cell<- read.csv(paste0(metpert_WTproteomics_dir,
                                          "/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                   stringsAsFactors = F)%>%
  dplyr::select(Genes, ORF,Element, median_relative_intracellular_concentration,median_log2_foldchangevsAE)%>%
  unique()


## count how many poorly characterised proteins we measure

poorly_char_ORFs <- unique(filter(GenProt_SGD, Uniprot.Annotation.Score < 3)$ORF)
all_meas_prots = unique(FCvsAE_allproteins_cell$ORF)
length(intersect(all_meas_prots,poorly_char_ORFs))

full_EC <- full_EC_SGD %>%
           mutate(GeneName = as.character(lapply(ORF,convert_ORF2SingleGeneName)))

allmetal_cluster_enrichments <- vector()

for(i in 1:length(unique(full_EC$Cluster))){
  
  Clus_df <- filter(FCvsAE_allproteins_env, Genes %in% filter(full_EC,Cluster == i)$GeneName)
  
  Clus_df_smry <- Clus_df%>%
             group_by(BioSpecID,Element.Concentration,Element)%>%
             summarise(mean_log2FCvsAE_piclust = mean(Mean_of_reps_Log2FC_vs_AE, na.rm = T))%>%
             ungroup()
  
  pdf(paste0(plot_dir,"/allmetal/clusterwise/mean_log2FCvsAE_all_proteins_cluster",i,".pdf"),width=10,height=6)
  print( ggplot(Clus_df_smry,
                aes(x = log2(Element.Concentration),
                    y = mean_log2FCvsAE_piclust,
                    colour = Element,
                    group=Element))+
           geom_point(size = 3)+
           geom_line(size = 1,alpha = 0.8)+
           scale_colour_manual(values = colkey_Ele)+
           theme_metallica()+
           labs(title = paste0("Cluster",i))
         
  )
  
  print(  ggplot(Clus_df_smry,
                 aes(x = log2(Element.Concentration),
                     y = mean_log2FCvsAE_piclust,
                     colour = Element,
                     group=Element))+
           geom_smooth(alpha = 0.1,method = "loess")+
           scale_colour_manual(values = colkey_Ele)+
           theme_metallica()+
           labs(title = paste0("Cluster ",i)))
  
  dev.off()
  
  nump = length(unique(Clus_df$Genes))
  
  pdf(paste0(plot_dir,"/allmetal/clusterwise/Mean log2FCvsAE per protein cluster",i,".pdf"),width = round(sqrt(nump),0)*5.6,height =round(sqrt(nump),0)*5)
  print(
    ggplot(Clus_df,
           aes(x=log2(Element.Concentration),
               y=Mean_of_reps_Log2FC_vs_AE,
               colour = Element,
               group=Element))+
      geom_point(size=3)+
      geom_line(size=1,alpha=0.7)+
      facet_wrap("Genes",ncol=round(sqrt(nump),0),scales="free")+
      scale_colour_manual(values = colkey_Ele)+
      theme_metallica()
  )
  
  # geom smooth
  print(
    ggplot(Clus_df,
                  aes(x=log2(Element.Concentration),
                      y=Mean_of_reps_Log2FC_vs_AE,
                      colour = Element,
                      group=Element))+
      geom_smooth(alpha = 0.1,method = "loess")+
      facet_wrap("Genes",ncol=round(sqrt(nump),0),scales="free")+
      scale_colour_manual(values = colkey_Ele)+
      theme_metallica()
  )
  dev.off()
  
  ## GO enrichment
  
  
  for_hyperGSA <- data.frame(ORF =  unique(full_EC$ORF))%>%
                  mutate(in_cluster = ifelse(ORF %in% Clus_df$ORF,1,0))
  
  hyperGSA_results_cluster <- run_HyperGSA(for_hyperGSA)
  
  if(!is.null(dim(hyperGSA_results_cluster))){
    hyperGSA_results_cluster$Cluster = i
    allmetal_cluster_enrichments = rbind(allmetal_cluster_enrichments,hyperGSA_results_cluster)
  }
}

write.csv(allmetal_cluster_enrichments,
          paste0(output_tables_dir,"/hyperGSA_enrichments_per_cluster_allmetalclustering.csv"), row.names = F)

#################################
### Summarise unchar proteins ###
#################################

DA_prots_env <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
              dplyr::select(Genes, Element, ORF, Significant)%>%
              filter(Significant == 1)%>%
              unique()

DA_prots_cell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),stringsAsFactors = F)%>%
              dplyr::select(Genes, Element, ORF, Significant)%>%
              filter(Significant == 1)%>%
              unique()


all_DA_proteins <- unique(rbind(DA_prots_cell, DA_prots_env))

# merge with differential abundance results
full_EC_SGD <- merge(full_EC, all_DA_proteins, by = "ORF")

full_EC_poorly_char <- full_EC_SGD%>%
                        filter(Uniprot.Annotation.Score  <3)

write.csv(full_EC_poorly_char,
          paste0(output_tables_dir,"/poorly_characterised_proteins_in_allmetal_clusters.csv"))


###############################
### ELEMENT WISE CLUSTERING ###
###############################

sig_metals <- unique(all_DA_proteins$Element)

metal_wise_clusters <- vector()

for(m in 1:length(sig_metals)){
  
  
  metal_clust_dir <- paste0(input_tables_dir,"/permetal/Output/",sig_metals[m])
  
  num_metal_clust <- length(dir(metal_clust_dir)[grep(".txt", dir(metal_clust_dir))])
  
  
  for(mc in 1:num_metal_clust){
    
    cl_uniprots = read.table(paste0(metal_clust_dir,"/C",mc,".txt"), stringsAsFactors = F)
    cl_uniprots$Cluster = mc
    cl_uniprots$metal_series = sig_metals[m]
    colnames(cl_uniprots)[[1]] = "Uniprot_ID"
    
    metal_wise_clusters <- rbind(metal_wise_clusters,cl_uniprots)
    
  }
}

## Get enrichments for each cluster in each metal series

metal_wise_cluster_enrichments <- vector()

for(m in 1:length(sig_metals)){
  
  mcs_w_enrich <- unique(filter(metal_wise_clusters,metal_series == sig_metals[m])$Cluster)
  
  for(mc in 1:length(mcs_w_enrich)){
  
    for_hyperGSA <- data.frame(Uniprot_ID = unique(metal_wise_clusters$Uniprot_ID))%>%
                    mutate(in_cluster = ifelse(Uniprot_ID %in% filter(metal_wise_clusters,
                                                                     metal_series == sig_metals[m] &
                                                                     Cluster == mcs_w_enrich[mc])$Uniprot_ID,1,0),
                          ORF = as.character(lapply(Uniprot_ID, convert_Uniprot2singleORF)))%>%
                    dplyr::select(ORF, in_cluster)
  
    hyperGSA_results_metalcluster <- run_HyperGSA(for_hyperGSA)
  
  if(!is.null(dim(hyperGSA_results_metalcluster))){
    hyperGSA_results_metalcluster$Cluster = mc
    hyperGSA_results_metalcluster$metal = sig_metals[m]
    metal_wise_cluster_enrichments = rbind(metal_wise_cluster_enrichments,hyperGSA_results_metalcluster)
    }
  }
}


metal_wise_clusters <- metal_wise_clusters %>%
  mutate(ORF = as.character(lapply(Uniprot_ID, convert_Uniprot2singleORF)),
         Gene = as.character(lapply(Uniprot_ID,convert_Uniprot2SingleGeneName)))

metal_wise_clusters <- merge(metal_wise_clusters,GenProt_SGD[,c("Uniprot.Annotation.Score","ORF")], by = "ORF")

Conc_Conversion <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
                   dplyr::select(BioSpecID,Element,Element.Concentration)%>%
                   unique()

####################################################################
### orphan metal binders and transporters per metal wise cluster ###
####################################################################

metal_wise_clusters_orphanmetb <- merge(metal_wise_clusters[,c("ORF", "Cluster","metal_series")],
                                        rbind(GO_gset_MF_metalbinding_unspecific, GO_gset_MF_binding_metalwise),
                                        by = "ORF", all.x = T)%>%
                                 group_by(Cluster, metal_series)%>%
                                 mutate(total = length(ORF),
                                        num_clusters = max(Cluster))%>%
                                 ungroup()%>%
                                 group_by(Cluster, metal_series,term,num_clusters)%>%
                                 mutate(percentage = 100*length(ORF)/total)%>%
                                 ungroup()%>%
                                 dplyr::select(Cluster,metal_series,term, percentage)%>%
                                 unique()%>%
                                 reshape2::dcast(Cluster + metal_series~term, value.var = "percentage")%>%
                                 dplyr::select(Cluster,metal_series,unspecific)%>%
                                 unique()

metal_wise_clusters_orphanmetb[is.na(metal_wise_clusters_orphanmetb)] <- 0

orphan_metb_metwise <- ggplot(metal_wise_clusters_orphanmetb,
                      aes(x=factor(Cluster),
                          y = metal_series,
                          fill = unspecific))+
  geom_tile(colour = "black", linewidth = 0.1)+
  theme_metallica()+
  theme(legend.position = "bottom")+
  scale_fill_viridis()+
  labs(fill = "% orphan metal binders",
       y = " ")

orphan_metb_metwise

pdf(paste0(plot_dir, "/metalwise/orphan_metal_binders_transportes_metalcluster_wise.pdf"),
width = 8,height = 4.5)
orphan_metb_metwise

dev.off()




### Plot protein quantities env and cellular of all DA proteins in clusters 

# env conc 

sig_metals_env <- unique(DA_prots_env$Element)
  
PQ_df_metalwise_clusters_env <- FCvsAE_allproteins_env%>% # Using fold change values from unimputed data to visualise results fo clustering
                          filter(Element %in% sig_metals)%>%
                          dplyr::select(Genes, ORF, Element, BioSpecID, Mean_of_reps_Log2FC_vs_AE)%>%
                          unique()

PQ_df_metalwise_clusters_env <- merge(PQ_df_metalwise_clusters_env, Conc_Conversion, by = c("BioSpecID","Element"))


PQ_df_metalwise_clusters_env <- merge(PQ_df_metalwise_clusters_env, 
                                  metal_wise_clusters[,c("ORF","Cluster","metal_series","Uniprot.Annotation.Score")], by = "ORF")%>%
                            filter(Element == metal_series)


PQ_df_metalwise_clusters_env <- PQ_df_metalwise_clusters_env%>%
                            group_by(Element,Genes)%>%
                            mutate(label = ifelse(grepl(pattern = "0.2",BioSpecID),
                                                  ifelse(Uniprot.Annotation.Score < 3,Genes,NA),NA))%>%
                            ungroup()



## Plot all element wise clusters 

pdf(paste0(plot_dir,"/metalwise/all_metalwise_clusters_noscalelim_envconc.pdf"),width = 30,height = 50)
ggplot(PQ_df_metalwise_clusters_env,
       aes(x = log2(as.numeric(as.character(Element.Concentration))),
           y = Mean_of_reps_Log2FC_vs_AE   ))+
  facet_wrap(c("metal_series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.25,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+
  geom_hline(yintercept = 0, size = 0.1)+
  geom_vline(xintercept = 0, size = 0.1)+
  geom_text(aes(label = label, colour = Uniprot.Annotation.Score),size = 3)+
  theme_metallica()
dev.off()



pdf(paste0(plot_dir,"/metalwise/all_metalwise_clusters_scalelim_envconc.pdf"),width = 30,height = 50)
ggplot(PQ_df_metalwise_clusters_env,
       aes(x = log2(as.numeric(as.character(Element.Concentration))),
           y = Mean_of_reps_Log2FC_vs_AE   ))+
  facet_wrap(c("metal_series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.5,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+ 
  geom_hline(yintercept = 0, size = 0.2)+
  geom_vline(xintercept = 0, size = 0.2)+
  ylim(-2,2)+
  geom_text(aes(label = label, colour = Uniprot.Annotation.Score),size = 2.5)+
  theme_metallica()+
  labs(x = "log2(metal concentration)", y = "log2(fold difference vs allele)")
dev.off()

# free scales 

PQ_df_metalwise_clusters_env <- filter(PQ_df_metalwise_clusters_env,
                                       !(Element == "Fe" & Mean_of_reps_Log2FC_vs_AE < -2))%>%
                                filter(! (Element == "Zn" & Mean_of_reps_Log2FC_vs_AE > 2.5))

pdf(paste0(plot_dir,"/metalwise/all_metalwise_clusters_freescales_envconc.pdf"),width = 30,height = 50)
ggplot(PQ_df_metalwise_clusters_env,
       aes(x = log2(as.numeric(as.character(Element.Concentration))),
           y = Mean_of_reps_Log2FC_vs_AE   ))+
  facet_wrap(c("metal_series","Cluster"),ncol = 6,scales = "free")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.2,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+ 
  geom_hline(yintercept = 0, size = 0.2)+
  geom_vline(xintercept = 0, size = 0.2)+
  geom_text(aes(label = label, colour = Uniprot.Annotation.Score),size = 2.5)+
  theme_metallica()+
  labs(x = "log2(metal concentration)", y = "log2(fold difference vs allele)")
dev.off()

## cellular conc

sig_metals_cellular <- unique(DA_prots_cell$Element)

PQ_df_metalwise_clusters_cell <- FCvsAE_allproteins_cell%>% # Using fold change values from unimputed data to visualise results fo clustering
  filter(Element %in% sig_metals_cellular)%>%
  dplyr::select(Genes, ORF, Element,median_relative_intracellular_concentration, median_log2_foldchangevsAE)%>%
  unique()


PQ_df_metalwise_clusters_cell <- merge(PQ_df_metalwise_clusters_cell, 
                                      metal_wise_clusters[,c("ORF","Cluster","metal_series","Uniprot.Annotation.Score")], by = "ORF")%>%
  filter(Element == metal_series)


PQ_df_metalwise_clusters_cell <- PQ_df_metalwise_clusters_cell%>%
  group_by(Element,Genes)%>%
  mutate(label = ifelse(median_log2_foldchangevsAE == min(median_log2_foldchangevsAE),
                        ifelse(Uniprot.Annotation.Score < 3,Genes,NA),NA))%>%
  ungroup()



## Plot all element wise clusters - cellular conc

pdf(paste0(plot_dir,"/metalwise/all_metalwise_clusters_noscalelim_cellconc.pdf"),width = 30,height = 50)
ggplot(PQ_df_metalwise_clusters_cell,
       aes(x = log2(as.numeric(as.character(median_relative_intracellular_concentration))),
           y = median_log2_foldchangevsAE   ))+
  facet_wrap(c("metal_series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.25,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+
  geom_hline(yintercept = 0, size = 0.1)+
  geom_vline(xintercept = 0, size = 0.1)+
  geom_text_repel(aes(label = label, colour = Uniprot.Annotation.Score),size = 3)+
  theme_metallica()+
  labs(x = "log2(metal concentration cellular)", y = "log2(fold difference vs allele)")
dev.off()



pdf(paste0(plot_dir,"/metalwise/all_metalwise_clusters_scalelim_cellconc.pdf"),width = 30,height = 50)
ggplot(PQ_df_metalwise_clusters_cell,
       aes(x = log2(as.numeric(as.character(median_relative_intracellular_concentration))),
           y = median_log2_foldchangevsAE   ))+
  facet_wrap(c("metal_series","Cluster"),ncol = 6,scales = "free_x")+
  geom_line(aes(colour = Uniprot.Annotation.Score, group = ORF),linewidth = 0.5,stat = "smooth", method = "loess", alpha = 0.7)+
  scale_colour_viridis_d(begin =0.2,end=0.9, option = "B")+ 
  geom_hline(yintercept = 0, size = 0.2)+
  geom_vline(xintercept = 0, size = 0.2)+
  ylim(-2,2)+
  geom_text_repel(aes(label = label, colour = Uniprot.Annotation.Score),size = 2.5)+
  theme_metallica()+
  labs(x = "log2(metal concentration cellular)", y = "log2(fold difference vs allele)")
dev.off()

#####################################
### Summarise metal wise clusters ###
#####################################

metalwise_cluster_numprots_summary <- metal_wise_clusters %>%
                                    dplyr::select(metal_series, Cluster, ORF,Uniprot.Annotation.Score)%>%
                                    unique()%>%
                                    group_by(metal_series, Cluster, Uniprot.Annotation.Score)%>%
                                    dplyr::summarise(num_protein = length(ORF))%>%
                                    ungroup()

pdf(paste0(plot_dir,"/metalwise/metal_wise_clustering_summary_env.pdf"),width = 16, height = 3)
ggplot(metalwise_cluster_numprots_summary,
       aes(x = Cluster, 
           y = num_protein,
           fill = Uniprot.Annotation.Score))+
  geom_bar(stat = "identity", position = "stack", size = 0, width = 0.6)+
  facet_wrap("metal_series", ncol = 4,scales = "free_x")+
  scale_fill_viridis_d(begin =0.2,end=0.9, option = "B")+
  theme_metallica()
dev.off()

PQ_df_metalwise_clusters_summary_unchar <- metal_wise_clusters %>%
                                           filter(Uniprot.Annotation.Score < 3)%>%
                                           dplyr::select(metal_series,ORF, Cluster, Uniprot.Annotation.Score)%>%
                                           unique()

pdf(paste0(plot_dir,"/metalwise/poorly_characterised_proteins_metal_wise_cluster_membership.pdf"),width = 15, height = 20)
ggplot(PQ_df_metalwise_clusters_summary_unchar,
       aes(x = ORF,
           y = factor(Cluster),
           fill = metal_series,
           colour = metal_series))+
  geom_tile()+
  scale_fill_manual(values = colkey_Ele)+
  scale_colour_manual(values = colkey_Ele)+
  coord_flip()+
  theme_metallica()+
  labs("cluster number")
dev.off()

#########################################################
### Merge dataset with GO MF metal and MetalPDB annos ###
#########################################################

## Plot summary of Metal related annotations ###

metal_wise_clusters_wanno_metalbinding_summary <- merge(metal_wise_clusters, GO_gset_MF_binding_metalwise,
                                                        by = "ORF")%>%
                                                  group_by(metal_series, Cluster, term)%>%
                                                  dplyr::summarise(num_protein = length(ORF))%>%
                                                  ungroup()

pdf(paste0(plot_dir,"/metalwise/metal_wise_clustering_metbinders_summary.pdf"),width = 16, height = 3)

ggplot(metal_wise_clusters_wanno_metalbinding_summary,
       aes(x = Cluster, 
           y = num_protein,
           fill = term))+
  geom_bar(stat = "identity", position = "stack", size = 0, width = 0.6)+
  facet_wrap("metal_series", ncol = 4)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(y = "metal binding annotations in each cluster")  

dev.off()

######################################################################################
### merge metal wise clustering with enrichment results for metal wise clustering ####
######################################################################################

metalwise_clusters_annotated <- merge(metal_wise_cluster_enrichments[,c("Cluster","metal","Gset.Term.Enriched", "Gset.Type")],
                                      metal_wise_clusters,by.x = c("Cluster","metal"), by.y = c("Cluster", "metal_series"),
                                      all.y = T)
#############################################################################################
### merge full allmetal clustering with enrichment results for full all metal clustering ####
#############################################################################################

allmetal_clusters_annotated <- merge( allmetal_cluster_enrichments[,c("Cluster","Gset.Term.Enriched","Gset.Type")], 
                                      full_EC[,c("Cluster","ORF","Uniprot.ID","GeneName")],
                                      by = "Cluster", all.y = T)

#############################################
### Combine Full and ele-wise clustering ####
#############################################
  
colnames(metalwise_clusters_annotated) <- c("cluster_metalwise","metal","Gset_term_enriched_metalwise","Gset_type_metalwise",
                                            "ORF", "Uniprot_ID","Gene","Uniprot.Annotation.Score")


colnames(allmetal_clusters_annotated) <- c("cluster_allmetal","Gset_term_enriched_allmetal","Gset_type_allmetal",
                                            "ORF", "Uniprot_ID","Gene")

all_ensemble_clustering_annotated <- merge(metalwise_clusters_annotated,allmetal_clusters_annotated,
                                           by = c("ORF","Uniprot_ID","Gene"), all = T)


write.csv(all_ensemble_clustering_annotated,
          paste0(output_tables_dir,"/ensemble_clustering_allresults_annotated.csv"),row.names = F)

poorly_characterise_all_ensemble_clustering_annotated <-
  unique(filter(all_ensemble_clustering_annotated,Uniprot.Annotation.Score < 3))

write.csv(poorly_characterise_all_ensemble_clustering_annotated,
          paste0(output_tables_dir,"/ensemble_clustering_allresults_annotated_poorly_charprots.csv"), row.names = F)


#######################################################################################################
### Count how many poorly characterised proteins are placed in clusters with functional enrichments ###
#######################################################################################################

# all metal clustering 

poorlychar_in_enriched_allm_clusters <- poorly_characterise_all_ensemble_clustering_annotated%>%
                                        dplyr::select(ORF, cluster_allmetal,Gset_term_enriched_allmetal)%>%
                                        unique()%>%
  filter(!is.na(cluster_allmetal))


enriched_df <- poorlychar_in_enriched_allm_clusters %>%
  group_by(ORF) %>%
  summarise(enriched = any(!is.na(Gset_term_enriched_allmetal)))

count_df <- enriched_df %>%
  summarise(num_enriched = sum(enriched))

# metalwise clustering

poorlychar_in_enriched_metalwise_clusters <- poorly_characterise_all_ensemble_clustering_annotated%>%
  dplyr::select(ORF, metal, cluster_metalwise,Gset_term_enriched_metalwise)%>%
  unique()%>%
  filter(!is.na(cluster_metalwise))


# Determine whether each ORF for each metal has a functional enrichment
enriched_metalwise_df <- poorlychar_in_enriched_metalwise_clusters %>%
  group_by(ORF, metal) %>%
  summarise(enriched = any(!is.na(Gset_term_enriched_metalwise)))

# Count for how many ORFs the cluster they were placed in is enriched in 
# at least one geneset term for each metal
count_metalwise_df <- enriched_metalwise_df %>%
  group_by(metal) %>%
  summarise(num_enriched = sum(enriched))


# Filter rows where cluster_metalwise is NA and where there is an enrichment
total_ORFs_with_enrichments_metalwise <- poorlychar_in_enriched_metalwise_clusters %>%
  filter(!is.na(cluster_metalwise) & !is.na(Gset_term_enriched_metalwise)) %>%
  distinct(ORF)

nrow(total_ORFs_with_enrichments_metalwise)


# number of ORFs in metalwise or allmetal clusters

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,
                     !(is.na(cluster_metalwise) & is.na(cluster_allmetal)))$ORF))

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,
                     !(is.na(cluster_metalwise) & is.na(cluster_allmetal)) &
                     !(is.na(Gset_term_enriched_allmetal) & is.na(Gset_term_enriched_metalwise)))$ORF))


# in metalwise clusters only 

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,!is.na(cluster_metalwise))$ORF))

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,
                     !(is.na(cluster_metalwise) | is.na(Gset_term_enriched_metalwise)))$ORF))

# allmetal

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,!is.na(cluster_allmetal))$ORF))

length(unique(filter(poorly_characterise_all_ensemble_clustering_annotated,
                     !(is.na(cluster_allmetal) | is.na(Gset_term_enriched_allmetal)))$ORF))

######################################################################################################################
### for the chosen example clusters - Fe Cluster 10 and all metal cluster 14 -- make aPEAR enrichment network plot ###
######################################################################################################################

library(AnnotationDbi)
library(S4Vectors)
library(org.Sc.sgd.db)

library(aPEAR)
library(clusterProfiler)

## ca cluster 6 


all_Cu_meas_prot <- filter(FCvsAE_allproteins_env,Element == "Cu")$ORF

Cu_cluster_2 <- read.table(paste0(input_tables_dir,"/permetal/Output/Cu//C2.txt"))%>%
  mutate(ORF = as.character(lapply(V1, convert_Uniprot2singleORF)))%>%
  pull(ORF)


enrich <- enrichGO(Cu_cluster_2, 
                   OrgDb = org.Sc.sgd.db, 
                   ont = 'BP',
                   keyType = "ORF",
                   universe = all_Ca_meas_prot)

pdf(paste0(plot_dir,"/metalwise/Cu_cluster_2_aPEARnetwork_BP.pdf"),width = 10,height =7)
enrichmentNetwork(filter(enrich@result,p.adjust < 0.05), 
                  drawEllipses = TRUE, 
                  colorBy = "p.adjust",
                  colorType = "pval",
                  fontSize = 3,
                  pCutoff = log(0.001))
dev.off()

all_Fe_meas_prot <- filter(FCvsAE_allproteins_env,Element == "Fe")$ORF

Fe_cluster_10 <- read.table(paste0(input_tables_dir,"/permetal/Output/Fe//C10.txt"))%>%
                mutate(ORF = as.character(lapply(V1, convert_Uniprot2singleORF)))%>%
                pull(ORF)


enrich <- enrichGO(Fe_cluster_10, 
                   OrgDb = org.Sc.sgd.db, 
                   ont = 'BP',
                   keyType = "ORF",
                   universe = all_Fe_meas_prot)

pdf(paste0(plot_dir,"/metalwise/Fe_cluster_10_aPEARnetwork_BP.pdf"),width = 10,height =7)
enrichmentNetwork(filter(enrich@result,p.adjust < 0.05), 
                       drawEllipses = TRUE, 
                       colorBy = "p.adjust",
                       colorType = "pval",
                       fontSize = 5,
                       pCutoff = log(0.001))
dev.off()



all_Mnmeas_prot <- filter(FCvsAE_allproteins_env,Element == "Mn")$ORF

Mn_cluster_3 <- read.table(paste0(input_tables_dir,"/permetal/Output/Mn//C3.txt"))%>%
  mutate(ORF = as.character(lapply(V1, convert_Uniprot2singleORF)))%>%
  pull(ORF)


enrich <- enrichGO(Mn_cluster_3, 
                   OrgDb = org.Sc.sgd.db, 
                   ont = 'CC',
                   keyType = "ORF",
                   universe = all_Mnmeas_prot)

pdf(paste0(plot_dir,"/metalwise/Mn_cluster_3_aPEARnetwork_CC.pdf"),width = 10,height =7)
enrichmentNetwork(filter(enrich@result,pvalue < 0.05), 
                  drawEllipses = TRUE, 
                  colorBy = "p.adjust",
                  colorType = "pval",
                  fontSize = 3,
                  pCutoff = log(0.001))
dev.off()



all_Znmeas_prot <- filter(FCvsAE_allproteins_env,Element == "Zn")$ORF

Zn_cluster_11 <- read.table(paste0(input_tables_dir,"/permetal/Output/Zn//C11.txt"))%>%
  mutate(ORF = as.character(lapply(V1, convert_Uniprot2singleORF)))%>%
  pull(ORF)


enrich <- enrichGO(Zn_cluster_11, 
                   OrgDb = org.Sc.sgd.db, 
                   ont = 'BP',
                   keyType = "ORF",
                   universe = all_Znmeas_prot)

pdf(paste0(plot_dir,"/Zn_cluster_11_aPEARnetwork_BP.pdf"),width = 10,height =7)
enrichmentNetwork(filter(enrich@result,p.adjust < 0.05), 
                  drawEllipses = TRUE, 
                  colorBy = "p.adjust",
                  colorType = "pval",
                  fontSize = 3,
                  pCutoff = log(0.001))
dev.off()



#############
### All metal clusters 


all_meas_prot <- full_EC$ORF

perform_enrichment <- function(cluster_number, ontology) {
  cluster <- full_EC %>%
    filter(Cluster == cluster_number) %>%
    mutate(ORF = as.character(lapply(Uniprot.ID, convert_Uniprot2singleORF))) %>%
    pull(ORF)
  
  enrich_result <- enrichGO(cluster, 
                            OrgDb = org.Sc.sgd.db, 
                            ont = ontology,
                            keyType = "ORF",
                            universe = all_meas_prot)
  return(enrich_result)
 
}


# Define a list of clusters and ontologies to iterate through
clusters <- unique(full_EC$Cluster)
ontologies <- c("CC", "MF", "BP")

pdf(paste0(plot_dir, "/allmetal/allmetal_clusters_aPEARnetwork.pdf"), width = 10, height = 7)

# Perform enrichment analysis for all clusters and ontologies
for (cluster_number in clusters) {
  for (ontology in ontologies) {
    
    enrich_result <- perform_enrichment(cluster_number, ontology)
    
    filtered_result <- filter(enrich_result@result, p.adjust < 0.05)
    
    # Check if the filtered result has more than 4 rows
    if (nrow(filtered_result) > 4) {
      tryCatch({
        
        plot_title = paste("cluster",cluster_number,"GO:",ontology)
        plot <- enrichmentNetwork(filtered_result, 
                                  drawEllipses = TRUE, 
                                  colorBy = "p.adjust",
                                  colorType = "pval",
                                  fontSize = 5,
                                  pCutoff = log(0.01))
        
        # Use ggtitle to add a title to the plot
        plot <- plot + ggtitle(plot_title)
        print(plot)
        
      }, error = function(e) {
        cat("Error in enrichmentNetwork for Cluster", cluster_number, "and GO-", ontology, ":", conditionMessage(e), "\n")
      })
    } else {
      cat("No significant enrichment for Cluster", cluster_number, "and GO-", ontology, "\n")
    }
  }
}

dev.off()
