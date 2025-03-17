#################################################################
### EDP Script - 2A - Statistical analysis with linear models ###
#################################################################

library(gprofiler2)
#### Script to make fit linear models along concentration gradient for each protein-element combination

#################
### Set Paths ###
#################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))


###################################
### Get Functions and variables ###
###################################

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/extracellularconc_vs_proteinabundance")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/extracellularconc_vs_proteinabundance") 
dir.create(output_tables_dir,recursive = T)

###################################################################
### Set Fold Change and Adj PValue threshold for entire script ###
###################################################################

fcthresh <- 1.5 # Magnitude of change that will be considered significant for DE
adjpvthresh <- 0.05

#########################################################
### Function to run edgeR along a pseudotime gradient ###
#########################################################


get_sigres_lms <- function(dfeleprot,gen,ele){
  
  dfeleprot <- na.omit(dfeleprot)
  
  if(length(unique(dfeleprot$Element.Concentration)) >= 4 & 
     length(na.omit(unique(dfeleprot$Log2.Protein.Quantity.))) >= 8 ) {
    fit_null <-  lm(Log2.Protein.Quantity. ~1,data=dfeleprot)
    
    fit_degree_1 <- lm(Log2.Protein.Quantity. ~ log2(Element.Concentration),data=dfeleprot)
    fit_degree_2 <-  lm(Log2.Protein.Quantity. ~ poly(log2(Element.Concentration),2),data=dfeleprot)
    fit_degree_3 <-  lm(Log2.Protein.Quantity. ~ poly(log2(Element.Concentration),3),data=dfeleprot)
    
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


###################################################################
### Convert BioSpecID labels to measured concentration gradient ###
###################################################################

MeasConcAxis <- read.csv(paste0(metpert_WTmetallomics_dir,"/output/tables/measured_element_concentration_axis.csv"),stringsAsFactors = F)%>%
  mutate(BioSpecID = gsub("_"," ",BioSpecID))

### Input data from all samples  with NAs ###

AllProtData <- read.csv(paste0(metpert_WTproteomics_dir,"/output/forStat/metperWTproteomics_protein_quantities_w_NAs_alldata_unfiltered_for_stat.csv"),
                        stringsAsFactors = F)%>%
  filter(!BioSpecID %in% c( "Mo 5","Mo 10") &
           Element != "K")%>%
  mutate(Genes = as.character(lapply(Protein.Ids,convert_Uniprot2SingleGeneName )))

AllProtData_AE <- filter(AllProtData,BioSpecID == "AllEle 1")

AllProtData <- merge(AllProtData,MeasConcAxis[,c("BioSpecID","actual_AEnorm_conc")],by="BioSpecID")  %>%
  mutate(Element.Concentration = round(actual_AEnorm_conc,3))%>%
  filter(Element.Concentration >=0.0002 & Element.Concentration < 100)

## Since The measured media concentration data does not have any values for  the AllElement condition, at this point when we merge the two datasets we lose all the AllElement control samples


eles <- unique(AllProtData$Element)

lm_res_df<-vector()


for(e in 1:length(eles)){
  
  ele_df <- filter(AllProtData,Element == eles[e]) [,c("Genes","Element.Concentration","Log2.Protein.Quantity.")]
  
  gens = unique(ele_df$Genes)
  
  for( g in 1:length(gens)){
    
    df_ep <- filter(ele_df,Genes==gens[g]) 
    lm_res_df<-rbind(lm_res_df,get_sigres_lms(dfeleprot = df_ep,gen = gens[g],ele = eles[e]))
  }
}


lm_res_df<-lm_res_df%>%
  group_by(ModelTest,Element)%>%
  mutate(PValAdj_BH = p.adjust(p.value,method="BH"))

################################################################
### Note down which is the simplest model that explains data ###
################################################################

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

# Summarize number of proteins-ele combinations per each last complex model type
pdf(paste0(plot_dir,"/number_of_proteins_by_least_complexmodel.pdf"),width=6,height=6)
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

#######################################################################################
### Do an enrichment of genes that are more likely to be linear cubic and quadratic ###
#######################################################################################
dir.create(paste0(plot_dir,"/gene_set_enrichments/gprofiler/"),recursive = T)
background_proteome <- as.character(unique( lm_res_df_ms$Genes))


DE_genes_linear_models <-  as.character(unique(filter(lm_res_df_ms, LeastComplexModel != "null")$Genes))
  
## overrepresentation
go_enrich_res <- gost(query = DE_genes_linear_models,organism = "scerevisiae", 
                        multi_query = T, custom_bg = background_proteome,
                        domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))

go_enrichres_plot <- gostplot(go_enrich_res, capped = TRUE, interactive = F)

go_enrich_res_df <- go_enrich_res$result

publish_gostplot(go_enrichres_plot, highlight_terms = go_enrich_res_df$term_id, 
                 width = 10, height = 10, 
                 filename = paste0(plot_dir,"/gene_set_enrichments/gprofiler/gprofiler_enrich_result_DE in any lm genes.pdf") )


## CUBIC shapes only 

## set background to all DE genes -- 
background_proteome <- DE_genes_linear_models

cubic_genes <- as.character(unique(filter(lm_res_df_ms, LeastComplexModel == "cubic")$Genes))

## overrepresentation
go_enrich_res <- gost(query = cubic_genes,organism = "scerevisiae", 
                      multi_query = T, custom_bg = background_proteome,
                      domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))


##### Quadratic

quadratic_genes <- as.character(unique(filter(lm_res_df_ms, LeastComplexModel == "quadratic")$Genes))

## overrepresentation
go_enrich_res <- gost(query = quadratic_genes,organism = "scerevisiae", 
                      multi_query = TRUE, custom_bg = background_proteome,measure_underrepresentation = F,
                      domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))

go_enrichres_plot <- gostplot(go_enrich_res, capped = TRUE, interactive = F)

go_enrich_res_df <- go_enrich_res$result

publish_gostplot(go_enrichres_plot, highlight_terms = go_enrich_res_df$term_id, 
                 width = 10, height = 30, 
                 filename = paste0(plot_dir,"/gene_set_enrichments/gprofiler/gprofiler_enrich_result_quadratic genes.pdf") )

### linear

linear_genes <- as.character(unique(filter(lm_res_df_ms, LeastComplexModel == "linear")$Genes))

## overrepresentation
go_enrich_res <- gost(query = linear_genes,organism = "scerevisiae", 
                      multi_query = TRUE, custom_bg = background_proteome,measure_underrepresentation = T,
                      domain_scope = "custom",sources = c("GO:MF","GO:BP","KEGG","TF"))



lm_res_df<- merge(lm_res_df,unique(lm_res_df_ms[,c("Genes","Element","LeastComplexModel")]),
                  by = c("Genes","Element"))
            
write.csv(lm_res_df,paste0(output_tables_dir,"/lm_fit_res_pvthresh_",adjpvthresh,".csv"),row.names=F)


####################################################
### Combine Ttest Results with log2PQ data frame ###
####################################################

AEdat <- AllProtData_AE %>%
  group_by(Genes)%>%
  summarize(Median.AE.PQ = median(Log2.Protein.Quantity.,na.rm=T))%>%
  as.data.frame()
rownames(AEdat)<-AEdat$Genes

PQ_lm_res_df <- merge(AllProtData[,c("Genes","Element","BioSpecID",
                                     "Element.Concentration","Log2.Protein.Quantity.")],
                      lm_res_df,by=c("Genes","Element"))%>%
  mutate(Median.AE.value = AEdat[Genes,"Median.AE.PQ"])%>%
  mutate(Log2FC_vs_AE = Log2.Protein.Quantity.-Median.AE.value)%>%
  group_by(Genes,Element,Element.Concentration)%>%
  mutate(Mean_of_reps_Log2FC_vs_AE = mean(Log2FC_vs_AE,na.rm=T))%>%
  ungroup()%>%
  group_by(Genes,Element)%>%
  mutate(MinLog2FC_vs_AE = min(Mean_of_reps_Log2FC_vs_AE,na.rm=T),
         MaxLog2FC_vs_AE = max(Mean_of_reps_Log2FC_vs_AE,na.rm=T),
         MaxMag_LogFC =ifelse( abs(max(MaxLog2FC_vs_AE)) > abs(min(MinLog2FC_vs_AE)),
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
  filter(ModelTest == ModelTest2Keep)%>%
  mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF)))

write.csv(PQ_lm_res_df,paste0(output_tables_dir,"/lmfit_DE_res_with_PQ_and_SigNotSig_AdjPVthresh_",
                              adjpvthresh," FCthresh ",fcthresh,".csv"))

######################################################
### Plot top 3 trajectories for each type of model ###
######################################################

top3_lm <- PQ_lm_res_df %>%
  filter( Significant ==1 & LeastComplexModel =="linear" )%>%
  dplyr::select(Element,PValAdj_BH,Genes)%>%
  unique()%>%
  top_n(-3,PValAdj_BH)

top3_lm <- merge(top3_lm,PQ_lm_res_df, by = c("Element","Genes"))

if(length(unique(top3_lm$Genes)) >3){
  
  G2k <-unique(top3_lm$Genes)[1:3]
  top3_lm<- filter(top3_lm,Genes %in% G2k)
}


pdf(paste0(plot_dir,"/top_3_most_significant_linear.pdf"),width=10,height=4)
print(ggplot(top3_lm,
             aes(x=log2(Element.Concentration),
                 y=Log2FC_vs_AE,
                 colour=BioSpecID))+
        facet_wrap("Genes",scales="free",ncol = 5)+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm",formula=y~x,colour="black")+
        scale_colour_manual(values=colkey_BioSpecID)+
        theme(legend.position = "none")+
        labs(title="Top 3 trajectories linear")+
        theme_metallica()+
        theme(legend.position = "none"))
dev.off()



top3_qd <- PQ_lm_res_df %>%
  filter( Significant ==1 & LeastComplexModel =="quadratic" )%>%
  dplyr::select(Element,PValAdj_BH,Genes)%>%
  unique()%>%
  top_n(-3,PValAdj_BH)

top3_qd <- merge(top3_qd,PQ_lm_res_df, by = c("Element","Genes"))

if(length(unique(top3_qd$Genes)) >3){
  
  G2k <-unique(top3_qd$Genes)[1:3]
  top3_qd<- filter(top3_qd,Genes %in% G2k)
}

pdf(paste0(plot_dir,"/top_3_most_significant_quadratic.pdf"),width=10,height=4)
print(ggplot(top3_qd,
             aes(x=log2(Element.Concentration),
                 y=Log2FC_vs_AE,
                 colour=BioSpecID))+
        facet_wrap("Genes",scales="free",ncol = 5)+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm",formula=y~poly(x,2),colour="black")+
        scale_colour_manual(values=colkey_BioSpecID)+
        theme(legend.position = "none")+
        labs(title="Top 3 trajectories quadratic")+
        theme_metallica()+
        theme(legend.position = "none"))
dev.off()

top3_cub <- PQ_lm_res_df %>%
  filter( Significant ==1 & LeastComplexModel =="cubic" )%>%
  dplyr::select(Element,PValAdj_BH,Genes)%>%
  unique()%>%
  top_n(-3,PValAdj_BH)

top3_cub <- merge(top3_cub,PQ_lm_res_df, by = c("Element","Genes"))

if(length(unique(top3_cub$Genes)) >3){
  
  G2k <-unique(top3_cub$Genes)[1:3]
  top3_cub<- filter(top3_cub,Genes %in% G2k)
}

pdf(paste0(plot_dir,"/top_3_most_significant_cubic.pdf"),width=10,height=4)
print(ggplot(top3_cub,
             aes(x=log2(Element.Concentration),
                 y=Log2FC_vs_AE,
                 colour=BioSpecID))+
        facet_wrap("Genes",scales="free",ncol = 5)+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm",formula=y~poly(x,3),colour="black")+
        scale_colour_manual(values=colkey_BioSpecID)+
        theme(legend.position = "none")+
        labs(title="Top 3 trajectories cubic")+
        theme_metallica()+
        theme(legend.position = "none"))
dev.off()

##############################################################################
### Plot all significant and top 20 trajectories per element and per model ###
##############################################################################
dir.create(paste0(plot_dir,"/significant_fits/"))

plot_allsig_permodel_perele <- function(df,modeltype,ele){
  
  df <- df%>%
    filter(Element == ele, LeastComplexModel == modeltype & Significant == 1)%>%
    arrange(PValAdj_BH)
  
  df$Genes = factor(df$Genes,levels = unique(df$Genes))
  
  nsiggen <- length(unique(df$Genes))
  
  if(nrow(df)>1){
    
    if(modeltype=="linear"){
      
      pdf(paste0(plot_dir,"/significant_fits/AllSig_",modeltype,"_fits_along",ele,"_concgrad.pdf"),width = round(sqrt(nsiggen),0)*2.5,
          height = round(sqrt(nsiggen),0)*2.6)
      print(ggplot(df,
                   aes(x=log2(Element.Concentration),
                       y=Log2FC_vs_AE,
                       colour=BioSpecID))+
              facet_wrap("Genes",scales="free",ncol = round(sqrt(nsiggen),0))+
              geom_point(size=3,alpha=0.8)+
              geom_smooth(method="lm",formula=y~x,colour="black")+
              scale_colour_manual(values=colkey_BioSpecID)+
              theme(legend.position = "none")+
              labs(title=paste(ele,modeltype))+
              theme_metallica()+
              theme(legend.position = "none"))
      dev.off()
      
      pdf(paste0(plot_dir,"/significant_fits/Top20_",modeltype,"_fits_along",ele,"concgrad.pdf"),width = 20,
          height = 20)
      print(ggplot(filter(df,Genes %in% unique(df$Genes)[1:20]),
                   aes(x=log2(Element.Concentration),
                       y=Log2FC_vs_AE,
                       colour=BioSpecID))+
              facet_wrap("Genes",scales="free",ncol = 5)+
              geom_point(size=3,alpha=0.8)+
              geom_smooth(method="lm",formula=y~x,colour="black")+
              scale_colour_manual(values=colkey_BioSpecID)+
              theme(legend.position = "none")+
              labs(title=paste(ele,modeltype))+
              theme_metallica()+
              theme(legend.position = "none"))
      dev.off()
      
      
      
    }else if(modeltype=="quadratic"){
      
      pdf(paste0(plot_dir,"/significant_fits/AllSig_",modeltype,"_fits along",ele,"_concgrad.pdf"),width = round(sqrt(nsiggen),0)*2.5,
          height = round(sqrt(nsiggen),0)*2.5)
      print( ggplot(df,
                    aes(x=log2(Element.Concentration),
                        y=Log2FC_vs_AE,
                        colour=BioSpecID))+
               facet_wrap("Genes",scales="free",ncol = round(sqrt(nsiggen),0))+
               geom_point(size=3,alpha=0.8)+
               geom_smooth(method="lm",formula=y~poly(x,2),colour="black")+
               scale_colour_manual(values=colkey_BioSpecID)+
               labs(title=paste(ele,modeltype))+
               theme_metallica()+
               theme(legend.position = "none")
      )
      dev.off()
      
      pdf(paste0(plot_dir,"/significant_fits/Top20_",modeltype,"_fits_along",ele,"_concgrad.pdf"),width = 20,
          height = 20)
      print(ggplot(filter(df,Genes %in% unique(df$Genes)[1:20]),
                   aes(x=log2(Element.Concentration),
                       y=Log2FC_vs_AE,
                       colour=BioSpecID))+
              facet_wrap("Genes",scales="free",ncol = round(sqrt(nsiggen),0))+
              geom_point(size=3,alpha=0.8)+
              geom_smooth(method="lm",formula=y~poly(x,2),colour="black")+
              scale_colour_manual(values=colkey_BioSpecID)+
              labs(title=paste(ele,modeltype))+
              theme_metallica()+
              theme(legend.position = "none")
      )
      dev.off()
      
    }else if(modeltype=="cubic"){
      
      pdf(paste0(plot_dir,"/significant_fits/AllSig_",modeltype,"_fits_along",ele,"_concgrad.pdf"),width = round(sqrt(nsiggen),0)*2.5,
          height = round(sqrt(nsiggen),0)*2.5)
      print(
        ggplot(df,
               aes(x=log2(Element.Concentration),
                   y=Log2FC_vs_AE,
                   colour=BioSpecID))+
          facet_wrap("Genes",scales="free",ncol = round(sqrt(nsiggen),0))+
          geom_point(size=3,alpha=0.8)+
          geom_smooth(method="lm",formula=y~poly(x,3),colour="black")+
          scale_colour_manual(values=colkey_BioSpecID)+
          labs(title=paste(ele,modeltype))+
          theme_metallica()+
          theme(legend.position = "none")
      )
      dev.off()
      
      pdf(paste0(plot_dir,"/significant_fits/Top20_",modeltype,"_fits_along",ele,"_concgrad.pdf"),width = 20,
          height = 20)
      print(
        ggplot(filter(df,Genes %in% unique(df$Genes)[1:20]),
               aes(x=log2(Element.Concentration), 
                   y=Log2FC_vs_AE,
                   colour=BioSpecID))+
          facet_wrap("Genes",scales="free",ncol = 5)+
          geom_point(size=3,alpha=0.8)+
          geom_smooth(method="lm",formula=y~poly(x,3),colour="black")+
          scale_colour_manual(values=colkey_BioSpecID)+
          labs(title=paste(ele,modeltype))+
          theme_metallica()+
          theme(legend.position = "none")
      )
      dev.off()
    }
  }
}

plot_df <- PQ_lm_res_df%>%
  filter(LeastComplexModel!="null")%>%
  mutate(ModelTest2Keep=ifelse(LeastComplexModel=="cubic","3_vs_0",
                               ifelse(LeastComplexModel =="quadratic",
                                      "2_vs_0","1_vs_0")))%>%
  filter(ModelTest2Keep == ModelTest)


for(e in 1:length(eles)){
  
  mts = c("linear","quadratic","cubic")
  
  for(m in 1:length(mts)){
    
    
    plot_allsig_permodel_perele(df = plot_df ,
                                ele = eles[e],
                                modeltype = mts[m])
  }
  
}

##########################################################
## Plot Top 25 hit across all metals and all fit types ###
##########################################################

top30_alleles<- filter(PQ_lm_res_df,LeastComplexModel != "null")%>%
  dplyr::select(Element,Genes,PValAdj_BH,Significant,LeastComplexModel)%>%
  filter(Significant ==1 )%>%
  unique()%>%
  top_n(-30,PValAdj_BH)
write.csv(top30_alleles,paste0(output_tables_dir,"/Top_30_across_all_metals.csv"),row.names = F)

####################################################
### Summarize and Compare hits across Conditions ###
####################################################

DElmfit_res_sum <- PQ_lm_res_df%>%
  dplyr::select(Genes,Element,Significant,MaxMag_LogFC)%>%
  filter(Element != "AllEle")%>%
  unique()%>%
  group_by(Element)%>%
  mutate(Total_Significant = sum(Significant),
         Total_Tested = sum(Significant) + sum(!Significant),
         Fraction_Significant = sum(Significant) / (sum(Significant) + sum(!Significant) ))%>%
  ## Add info on how many conditions something is significant in 
  ungroup()%>%
  group_by(Genes)%>%
  mutate(Total_Conditions_Measured_in = length(Element),
         Total_Conditions_Significant_in = sum(Significant))

## Plot number of hits

pdf(paste0(plot_dir,"/number_of_hits_per_condition_lmDEres.pdf"),width=10,height=8)
ggplot(unique(DElmfit_res_sum[,c("Element","Total_Significant", "Total_Tested","Fraction_Significant")]),
       aes(x=Element))+
  geom_bar(aes(y=Total_Tested),
           fill=NA,
           colour="navyblue",
           stat="identity")+
  geom_bar(aes(y=Total_Significant),
           colour=NA,
           fill="navyblue",
           stat="identity")+
  geom_text(aes(y=Total_Significant+200,
                label=paste0("Tested:\n",Total_Tested,"\n","Significant:\n",Total_Significant)))+
  labs(title=paste0("DE genes padj < ",adjpvthresh, "\n abs(FC max minus min PQ) > ",fcthresh),
       y="Number of Genes")+
  theme_metallica()

dev.off()


DElmfit_res_sum_numsigcat <- PQ_lm_res_df%>%
  dplyr::select(Genes,Element,Significant,MaxMag_LogFC)%>%
  filter(Element != "AllEle")%>%
  unique()%>%
  group_by(Genes)%>%
  mutate(Fraction_Conditions_Significant = round(sum(Significant)/length(Element),1)) %>%
  ungroup()%>%
  filter(Significant == 1)%>%
  group_by(Element,Fraction_Conditions_Significant)%>%
  summarize(Total_Sig_per_EleFracSig = sum(Significant),
  )%>%
  ungroup()



####################################################################################
### Function to plot volcano plots and top 25 trajectories per metal for each lm ###
####################################################################################

plot_lmfitDEs<-function(df,ele){
  
  
  dfp <- filter(df, Element == ele)[,c("Genes","MaxMag_LogFC",
                                       "PValAdj_BH","Log2.Protein.Quantity.",
                                       "Element.Concentration","BioSpecID",
                                       "Significant")]
  
  ## Plot Summary Volcano plot 
  
  vp <- ggplot(unique(dfp[,c("Genes","MaxMag_LogFC","PValAdj_BH","Significant")]),
               aes(x=MaxMag_LogFC,
                   y=-log10(PValAdj_BH),
                   colour=factor(Significant)))+
    geom_point(size=3,alpha=0.6)+
    geom_text_repel(aes(label=Genes),colour="black")+
    scale_colour_manual(values=c("gray","darkred"))+
    theme_metallica()+
    theme(legend.position = "none")+
    labs(title=paste0("Differentially Expressed Genes ",ele))
  
  
  
  pdf(paste0(plot_dir,"/all_significant/DE_lm_",elename, i,"_degree_lm_summary.pdf"),width=8,height=8)
  print(vp)
  dev.off()
  
}

