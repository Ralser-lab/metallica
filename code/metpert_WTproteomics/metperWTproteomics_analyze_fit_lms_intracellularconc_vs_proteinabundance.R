#`---
#`  Title: "metpertWTproteomics_analyze_fit_lms_intracellularconc_vs_proteinabundance.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 30 May 2023 
#`  Description: Script to fit linear models to protein abundance data along intracellular metal concentraion

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# experiment specific


get_sigres_lms <- function(dfeleprot,gen,ele){
  
  dfeleprot <- na.omit(dfeleprot)
  colnames(dfeleprot) <- c("Genes","element_concentration","log2_protein_abundance")
  
  if(length(unique(dfeleprot$element_concentration)) >= 4 & 
     length(na.omit(unique(dfeleprot$log2_protein_abundance))) >= 8 ) {
    fit_null <-  lm(log2_protein_abundance ~1,data=dfeleprot)
    
    fit_degree_1 <- lm(log2_protein_abundance ~ log2(element_concentration),data=dfeleprot)
    fit_degree_2 <-  lm(log2_protein_abundance ~ poly(log2(element_concentration),2),data=dfeleprot)
    fit_degree_3 <-  lm(log2_protein_abundance ~ poly(log2(element_concentration),3),data=dfeleprot)
    
    aov_3vs0 <- broom::tidy(anova(fit_degree_3,fit_null))[2,c("statistic","p.value")]
    aov_3vs0$ModelTest <- "3_vs_0"
    aov_3vs0$Adj.R2.of.Higher.DF <- summary(fit_degree_3)$adj.r.squared
    
    aov_3vs1 <- broom::tidy(anova(fit_degree_3,fit_degree_1))[2,c("statistic","p.value")]
    aov_3vs1$ModelTest <- "3_vs_1"
    aov_3vs1$Adj.R2.of.Higher.DF <- summary(fit_degree_3)$adj.r.squared
    
    aov_3vs2 <- broom::tidy(anova(fit_degree_3,fit_degree_2))[2,c("statistic","p.value")]
    aov_3vs2$ModelTest <- "3_vs_2"
    aov_3vs2$Adj.R2.of.Higher.DF <- summary(fit_degree_3)$adj.r.squared
    
    
    aov_2vs0 <- broom::tidy(anova(fit_degree_2,fit_null))[2,c("statistic","p.value")]
    aov_2vs0$ModelTest <- "2_vs_0"
    aov_2vs0$Adj.R2.of.Higher.DF <- summary(fit_degree_2)$adj.r.squared
    
    aov_2vs1 <- broom::tidy(anova(fit_degree_2,fit_degree_1))[2,c("statistic","p.value")]
    aov_2vs1$ModelTest <- "2_vs_1"
    aov_2vs1$Adj.R2.of.Higher.DF <- summary(fit_degree_2)$adj.r.squared
    
    
    aov_1vs0 <- broom::tidy(anova(fit_degree_1,fit_null))[2,c("statistic","p.value")]
    aov_1vs0$ModelTest <- "1_vs_0"
    aov_1vs0$Adj.R2.of.Higher.DF <- summary(fit_degree_1)$adj.r.squared
    
    Ftest_res<-rbind(aov_3vs0,aov_3vs1,aov_3vs2,
                     aov_2vs0,aov_2vs1,
                     aov_1vs0)
    Ftest_res$Genes=gen
    Ftest_res$Element=ele
  }else{
    Ftest_res <- c(statistic=NA,p.value=NA,ModelTest=NA,Adj.R2.of.Higher.DF=NA,Genes=gen,Element=ele)
  }
  
  return(Ftest_res)
}

#################
### set paths ###
#################

# outputs

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/metallomics_proteomics_coanalysis")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/metallomics_proteomics_coanalysis")
dir.create(outputtables_dir, recursive = T)
  
  
##################
### Input data ###
##################

# protein abundances and relative ( fold difference values)

proteomics_data <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
                    dplyr::select(BioSpecID,Element, Genes,ORF, Log2.Protein.Quantity.,Log2FC_vs_AE)%>%
                    unique()

# measured intracellular relative metal concentrations

met_quants_df <- read.csv(paste0(metpert_WTmetallomics_dir,"/results/output_tables/metpertWTmetallomics_Pnorm_AEnorm.csv"),stringsAsFactors = F)%>%
                 filter(!element_measured %in% c("P","S") & BioSpecID != "AllEle" & !(element_measured == "Cu" & Ratio_to_AEngperwell > 2.10))%>%
                 dplyr::select(BioSpecID, element_measured, Ratio_to_AEngperwell)%>%
                 group_by(BioSpecID, element_measured)%>%
                 summarise(median_relative_intracellular_concentration = median(Ratio_to_AEngperwell,na.rm = T))%>%
                 ungroup()

# Mmrge proteomics and metallomics datasets 

prot_met_df <- merge(proteomics_data, met_quants_df, by = "BioSpecID")%>%
                ## convert intracellular concentrations into a smaller range of bins 
                mutate(median_relative_intracellular_concentration = round(median_relative_intracellular_concentration,2))%>%
                group_by(Genes,ORF, element_measured,median_relative_intracellular_concentration)%>%
                dplyr::summarise(median_log2_protein_abundance = median(Log2.Protein.Quantity.,na.rm = T),
                                 median_log2_foldchangevsAE = median(Log2FC_vs_AE,na.rm = T))%>%
                ungroup()


############################################################################################################################
### fit linear models to binned protein abundance data vs measured and binned intracellular concentration of metals data ###
############################################################################################################################

metals <- unique(prot_met_df$element_measured)

lm_res_df<-vector()


for(m in 1:length(metals)){
  
  metal_df <- filter(prot_met_df,element_measured == metals[m]) [,c("Genes","median_relative_intracellular_concentration","median_log2_protein_abundance")]
  
  gens = unique(metal_df$Genes)
  
  for( g in 1:length(gens)){
    
    df_ep <- filter(metal_df,Genes==gens[g]) 
    lm_res_df<-rbind(lm_res_df,get_sigres_lms(dfeleprot = df_ep,gen = gens[g],ele = metals[m]))
  }
}


lm_res_df <- lm_res_df%>%
             group_by(ModelTest,Element)%>%
             mutate(PValAdj_BH = p.adjust(p.value,method="BH"))

################################################################
### Note down which is the simplest model that explains data ###
################################################################

adjpvthresh <- 0.05
fcthresh <- 1.5

lm_res_df_ms <- na.omit(lm_res_df)%>%
                reshape2::dcast(Genes+Element ~ ModelTest, value.var= "PValAdj_BH")%>%
                mutate(LeastComplexModel = ifelse( `1_vs_0` < adjpvthresh &
                                                     `2_vs_1` >= adjpvthresh,"linear",
                                                   ifelse(
                                                     `2_vs_0` < adjpvthresh &
                                                       `3_vs_2` >= adjpvthresh,"quadratic",
                                                     ifelse(
                                                       `3_vs_0` < adjpvthresh,"cubic",
                                                       "null"
                                                     )
                                                   )
                ))



lm_res_df_ms$LeastComplexModel = factor(lm_res_df_ms$LeastComplexModel,
                                        levels=c("null","linear","quadratic","cubic"))


###########################################################
### Visualise how many of each type of model was chosen ###
###########################################################

pdf(plot_dir,"/intracellular_metalconc_vs_proteinabundance_number_proteins_by_least_complex_model.pdf",width=6,height=6)
ggplot(lm_res_df_ms,
       aes(x=Element,
           group=LeastComplexModel,
           fill= LeastComplexModel))+
  geom_bar(stat="count",position="stack",alpha=0.9, width = 0.5)+
  scale_fill_viridis_d(begin=0.1,end=0.9,direction = -1)+
  labs(title=paste("Pvalue of Ftest:",adjpvthresh), x= "")+
  theme_metallica()+
  theme(legend.position="bottom")
dev.off()

############################################################################
### Combine least complex model annotation with statistical test results ###
############################################################################

lm_res_df <- merge(lm_res_df,unique(lm_res_df_ms[,c("Genes","Element","LeastComplexModel")]),
                  by = c("Genes","Element"))

############################################################################
### Combine results with Protein quantities and assign significance hits ###
############################################################################

lmfit_results_w_FCvsAE <- merge(lm_res_df,prot_met_df,
                                by.x = c("Genes","Element"),
                                by.y = c("Genes","element_measured"))%>%
                          group_by( Genes,Element)%>%
                          mutate( MinLog2FC_vs_AE = min(median_log2_foldchangevsAE,na.rm=T),
                                  MaxLog2FC_vs_AE = max(median_log2_foldchangevsAE,na.rm=T),
                                  MaxMag_LogFC = ifelse( abs(max(MaxLog2FC_vs_AE)) > abs(min(MinLog2FC_vs_AE)),
                                                         max(MaxLog2FC_vs_AE),
                                                         min(MinLog2FC_vs_AE)))%>%
                          ungroup()%>%
                          mutate(Significant = ifelse(PValAdj_BH < adjpvthresh & Adj.R2.of.Higher.DF > 0.3 &
                                                          ( abs(MaxLog2FC_vs_AE - MinLog2FC_vs_AE) > log2(fcthresh) ),
                                                        1 ,0))%>%
                          mutate(ModelTest2Keep=ifelse(LeastComplexModel=="cubic","3_vs_0",
                                                         ifelse(LeastComplexModel =="quadratic","2_vs_0",
                                                                ifelse(LeastComplexModel =="linear","1_vs_0",
                                                                       "3_vs_0"))))%>%
                          filter(ModelTest == ModelTest2Keep)
  

write.csv(lmfit_results_w_FCvsAE,paste0(outputtables_dir,"/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_",
                              adjpvthresh,"_FCthresh_",fcthresh,".csv"), row.names = F)

######################################################
### Plot top 3 trajectories for each type of model ###
######################################################

# For the linear model

top3_linear <- lmfit_results_w_FCvsAE %>%
                filter(LeastComplexModel == "linear" & Significant == 1) %>%
                dplyr::select(Element,PValAdj_BH,Genes)%>%
                unique()%>%
                top_n(-3,PValAdj_BH)

pdf(paste0(plot_dir, "/top_3_linear_modelfits_intracellvsprot.pdf"), width = 10, height = 4)
print(ggplot(filter(prot_met_df, Genes %in% unique(top3_linear$Genes) & element_measured %in% unique(top3_linear$Element)),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = element_measured)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        facet_wrap(element_measured~Genes, scales = "free", ncol = 5) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories - Linear Model", 
             x = "Intracellular Metal Concentration", y = "Log2 Fold Change vs. AE"))
dev.off()

# For the quadratic model
top3_quadratic <- lmfit_results_w_FCvsAE %>%
                  filter(LeastComplexModel == "quadratic" & Significant == 1) %>%
                  dplyr::select(Element,PValAdj_BH,Genes)%>%
                  unique()%>%
                  top_n(-3,PValAdj_BH)


pdf(paste0(plot_dir, "/top_3_quadratic_modelfits_intracellvsprot.pdf"), width = 10, height = 4)
print(ggplot(filter(prot_met_df, Genes %in% unique(top3_quadratic$Genes) & 
                      element_measured %in% unique(top3_quadratic$Element)),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = element_measured)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula =y ~ poly(x, 2), color = "black") +
        facet_wrap(element_measured~Genes, scales = "free", ncol = 5) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories - Quadratic Model", x = "Intracellular Metal Concentration", y = "Log2 Fold Change vs. AE"))
dev.off()

# For the cubic model
top3_cubic <- lmfit_results_w_FCvsAE %>%
              filter(LeastComplexModel == "cubic" & Significant == 1) %>%
              dplyr::select(Element,PValAdj_BH,Genes)%>%
              unique()%>%
              top_n(-3,PValAdj_BH)

if(length(unique(top3_cubic$Genes)) >3){
  
  G2k <-unique(top3_cubic$Genes)[1:3]
  top3_cubic<- filter(top3_cubic,Genes %in% G2k)
}

pdf(paste0(plot_dir, "/top_3_cubic_modelfits_intracellvsprot.pdf"), width = 10, height = 4)
print(ggplot(filter(prot_met_df, Genes %in% unique(top3_cubic$Genes) & 
                      element_measured %in% unique(top3_cubic$Element)),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = element_measured)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula =y ~ poly(x, 3), color = "black") +
        facet_wrap(element_measured~Genes, scales = "free", ncol = 5) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories - Cubic Model", x = "Intracellular Metal Concentration", y = "Log2 Fold Change vs. AE"))
dev.off()

                              

