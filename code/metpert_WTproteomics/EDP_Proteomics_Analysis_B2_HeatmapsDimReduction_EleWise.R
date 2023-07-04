##############################################################################
### EDP Script - 2B - Heatmaps & Dimensionality Reduction ### Element Wise ###
##############################################################################

# Need to run EDP_Proteomics_Analysis_ 0 & 1B before running this 
# This script will do the following :
# 1. 
# 2. 
# 3. 

require(ComplexHeatmap)
require(RColorBrewer)

# PCA plotting
require(PCAtools)

# umap and plotting
require(umap)
require(IDPmisc)


###########################################################################
### Function for Plotting Heatmap PCA and UMAP for each type of scaling ###
###########################################################################


getplot_heatmap_PCA_UMAP_perele <-function(scaled_df,scaletype="Log2FC",dftype = "whole"){
  
  setwd(perele_hTMP_pca_umap_dir)
  
  ###########################################
  ### Look at distribution of scaled data ###
  ###########################################
  
  rangetest<-scaled_df%>%
    mutate(Genes = rownames(scaled_df))%>%
    reshape2::melt(id.vars="Genes",value.name=scaletype,variable.name="BioSpecID")
  
  pdf(paste0(dftype," dataset ",scaletype,"scaled histogram all data.pdf"),width=8,height=5)
  
  ggplot(rangetest,
         aes(x=!!sym(scaletype)))+
    geom_histogram(bins=200,fill="maroon",colour="white")+
    #  facet_wrap("BioSpecID")+
    theme_SKA()
  dev.off()
  
  
  ##################################################
  ### Make Row and Column Annotation for Heatmap ###
  ##################################################
  
  
  row_anno.df <- merge(cbind(Genes=unique(rownames(scaled_df))),x,
                       by="Genes",all.x=T)%>%
    mutate(`ORF Type`=gsub("ORF[|]","",`ORF Type`))
  rownames(row_anno.df)<-row_anno.df$Genes
  row_anno.df<-row_anno.df[,-which(colnames(row_anno.df) %in% c("ORF","Genes"))]
  
  orftypecol=row_anno.df$`ORF Type`
  row_anno.df[!is.na(row_anno.df)] <- "1"
  row_anno.df[is.na(row_anno.df)] <- "0"
  row_anno.df$`ORF Type`=orftypecol
  
  row_anno <- rowAnnotation(df = row_anno.df,
                            col = list(
                              Mitochondria = c("0" ="#FFFFFF",
                                               "1" = RColorBrewer::brewer.pal(8,name="Dark2")[[1]]),
                              `Endoplasmic Reticulum`= c( "0"="#FFFFFF",
                                                          "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[2]] ),
                              Golgi = c( "0"="#FFFFFF",
                                         "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[3]] ),
                              Nucleus = c( "0"="#FFFFFF",
                                           "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[4]] ),
                              Peroxisome = c( "0"="#FFFFFF",
                                              "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[5]] ),
                              Vacuole = c( "0"="#FFFFFF",
                                           "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[6]] ),
                              `Metal related GO MF` = c( "0"="#FFFFFF",
                                                         "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[7]] ),
                              
                              GO_BP_metabolic_process_anno=c("0"="#FFFFFF",
                                                             "1"= RColorBrewer::brewer.pal(8,name="Dark2")[[8]] ),
                              Ca = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[1]]),
                              Cu = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[2]]),
                              Fe = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[3]]),
                              K = c("0"="#FFFFFF",
                                    "1"=element2colour()[,"Colour"][[4]]),
                              Mg = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[5]]),
                              Mn = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[6]]),
                              Mo = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[7]]),
                              Na = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[8]]),
                              Zn = c("0"="#FFFFFF",
                                     "1"=element2colour()[,"Colour"][[9]]),
                              `ORF Type`= c("Verified" = "#D9D9D9",
                                            "Uncharacterized"="#030303",
                                            "Dubious"="#030303")
                            ))
  
  ## Colour for heatmap cells
  
  ## Make annotation data frame for columns
  
  col_anno.df <- as.data.frame(cbind(BioSpecID = colnames(scaled_df)))%>%
    separate(BioSpecID, into = c("Element","Element.Concentration"),
             sep=" ",remove=F)%>%
    mutate(Element.Concentration = ifelse(Element %in% c("AllEleNew","AllEle"),1,
                                          Element.Concentration))%>%
    mutate(Direction = ifelse(Element.Concentration < 1 ,"Depletion",
                              ifelse(Element.Concentration > 1, "Excess",
                                     "Control")),
           Ele_Dir = paste(Element,Direction))
  
  rownames(col_anno.df)<-col_anno.df$BioSpecID
  col_anno.df <- col_anno.df[,c("BioSpecID","Ele_Dir")]
  
  col_anno = HeatmapAnnotation(df=col_anno.df,
                               col = list ( BioSpecID = colkey_BioSpecID,
                                            Ele_Dir = colkey_EleDir)
                               )
  
  # row clustering 
  dendr <- hclust(dist(scaled_df))
  dendr <- color_branches(dendr, k = 30)
  
  
  if(scaletype =="Log2FC"){ col_fun = colorRamp2(c(-4,0,4), c("blue","white","red"))}
  if(scaletype =="Zscore_Vs_AE"){ col_fun = colorRamp2(c(-4,0,4), c("blue","white","red"))}
  if(scaletype =="Zscore_Vs_AE_AEstats"){ col_fun = colorRamp2(c(-4,0,4), c("blue","white","red"))}
  if(scaletype =="Zscore_Vs_AS"){ col_fun = colorRamp2(c(-4,0,4), c("blue","white","red"))}
  
  
  pdf(paste0("Heatmap All Proteomics Data Imp Values Hclust ",dftype," dataset ",scaletype,".pdf"),width = 10,height = 110)
  
  draw(ComplexHeatmap::Heatmap(as.matrix(scaled_df),
                               col=col_fun,
                             
                               cluster_rows = dendr,
                               use_raster = T,
                               raster_device = "CairoPNG",
     
                               row_split = 30,
                               row_gap = unit(3, "mm"),
                               row_names_gp = gpar(fontsize=2),   
                               
                               column_names_gp = gpar(fontsize = 11),
                               column_gap = unit(2, "mm"),
                               column_names_side = "top",
                               
                               heatmap_legend_param = list(title = scaletype),
                               top_annotation = col_anno,
                               right_annotation = row_anno
  ),
  merge_legend = TRUE
  
  )
  
  dev.off()
  
  #####################################################################
  #############  PCA and UMAP with proteins as dimensions ############# 
  #####################################################################
  
  metadata_pca <- as.data.frame(cbind(BioSpecID = colnames(scaled_df)))%>%
    separate(BioSpecID, into = c("Element","Element.Concentration"),
             sep=" ",remove=F)%>%
    mutate(Direction = ifelse(Element.Concentration < 1, "Depletion",
                              ifelse (Element.Concentration > 1, "Excess","Control")),
           EleDir = paste(Element,Direction))
  
  rownames(metadata_pca) <- metadata_pca$BioSpecID
  
  # matrix for PCA = scaled_df
  
  all(colnames(scaled_df) == rownames(metadata_pca))     
  
  ## do PCA
  
  pca_res <- pca(scaled_df,metadata_pca)
  
  horn <- parallelPCA(scaled_df)
  horn$n
  
  elbow <- findElbowPoint(pca_res$variance)
  
  ## scree plot 
  
  scree<-PCAtools::screeplot(pca_res,
                   components = getComponents(pca_res, 1:20),
                   axisLabSize = 18, titleLabSize = 22,  
                   vline = c(horn$n,elbow))+
    
    geom_label(aes(x = horn$n + 1, y = 50,
                   label = 'Horn\'s', vjust = -0.5, size = 8))+
    geom_label(aes(x = elbow + 1, y = 50,
                   label = 'Elbow method', vjust = -3, size = 8))+
    geom_label(aes(x=16,y=80,label=paste("Num PCs for 80% var: ",
                                         which(cumsum(pca_res$variance) > 80)[1])))
  
  
  pp<-PCAtools::pairsplot(pca_res,
                components = getComponents(pca_res, c(1:4)),
                triangle = TRUE, trianglelabSize = 12,
                gridlines.major = FALSE, gridlines.minor = FALSE,
                colby = 'BioSpecID',
                title = 'Pairs plot', plotaxes = FALSE,
                margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
  
  ## biplot
  
  pc1_vs_pc2<-PCAtools::biplot(pca_res,
                               lab = pca_res$metadata$BioSpecID,
                               colby = 'BioSpecID',
                               colkey = colkey_BioSpecID,
                               showLoadings = T,
                               pointSize = 4,
                               colLoadingsNames = "red",
                               legendPosition = 'none')
  
  
  pc2_vs_pc3<-PCAtools::biplot(pca_res,
                               x="PC2",
                               y="PC3",
                               lab = pca_res$metadata$BioSpecID,
                               colby = 'BioSpecID',
                               colkey = colkey_BioSpecID,
                               showLoadings = T,
                               colLoadingsNames = "red",
                               pointSize = 4,
                               legendPosition = 'none')
  
  
  pc1_vs_pc3<-PCAtools::biplot(pca_res,
                     x="PC1",
                     y="PC3",
                     lab = pca_res$metadata$BioSpecID,
                     colby = 'BioSpecID',
                     colkey = colkey_BioSpecID,
                     showLoadings = T,
                     colLoadingsNames = "red",
                     pointSize = 4,
                     legendPosition = 'none')
  
  
  
  loadingPlot<-plotloadings(pca_res,components = c("PC1","PC2","PC3","PC4","PC5","PC6"),
  labSize=6,
  rangeRetain=0.01,
  shapeSizeRange=c(10,10),
  caption= "Top5% variables",
  drawConnectors = F,
  labhjust=-.5,
  gridlines.major = F,
  gridlines.minor=F,
  col = c(viridis(3)[[1]],
          viridis(3)[[2]], viridis(3)[[3]]))
  
  
  pca_loading_df<- pca_res$loadings%>% 
    mutate(Genes = rownames(pca_res$loadings))%>%
    mutate(ORF=as.character(lapply(Genes,convert_GeneName2ORF )))%>%
    reshape2::melt(id.vars=c("Genes","ORF"),
                   variable.name = 'PC',value.name = "PC loading")%>%
    mutate(PC = as.numeric(gsub("PC","",PC)))%>%
    filter(PC < 20)
  
  ### Run HyperGSA on top 1% genes in top 5 loadings ###
  
  setwd(HyperGSA_onPCA_elewise_dir)
  
  PCA_loadings_HyperGSA.df <- data.frame()
  Genes_top5_loadings_df <- data.frame()
  
  for(p in 1:5){
    
    df2bind<-data.frame()
    gene_df<-data.frame()
    
    top1pcnt = filter(pca_loading_df,PC == p)%>%
      top_frac(0.01,`PC loading`)%>%
      pull(ORF)
    
    top20 = filter(pca_loading_df,PC == p)%>%
      top_n(20,`PC loading`)%>%
      pull(ORF)
    
    gene_df<-rbind(gene_df,cbind(ORF=top20,PC=p))
    
    top_for_HyperGSA = pca_loading_df%>%
      mutate(UseforHyperGSA = ifelse(ORF %in% top1pcnt,1,0))%>%
      dplyr::select(ORF,UseforHyperGSA)
    
    HypGsa_Topres<-get_plot_HyperGSA(top_for_HyperGSA,paste0("Top1pcnt_fullPCA_PC",p," ",dftype," dataset ",scaletype),
                                     EnrichPV_thresh = 0.05)
    
    if(is.data.frame(HypGsa_Topres)){
      
      HypGsa_Topres$PC =p
      HypGsa_Topres$Loading_Direction ="Top1pcnt"
      df2bind<-rbind(df2bind,HypGsa_Topres)
      
    }
    
    bottom1pcnt = filter(pca_loading_df,PC == p)%>%
      top_frac(-0.1,`PC loading`)%>%
      pull(ORF)
    
    bottom20 = filter(pca_loading_df,PC == p)%>%
      top_n(-20,`PC loading`)%>%
      pull(ORF)
    gene_df<-rbind(gene_df,cbind(ORF=bottom20,PC=p))
    
    bottom_for_HyperGSA = pca_loading_df%>%
      mutate(UseforHyperGSA = ifelse(ORF %in% bottom1pcnt,1,0))%>%
      dplyr::select(ORF,UseforHyperGSA)
    
    HypGsa_Botres<-get_plot_HyperGSA(bottom_for_HyperGSA,paste0("Bottom1pcnt_fullPCA_PC",p," ",dftype," dataset ",scaletype)
                                     ,EnrichPV_thresh = 0.05)
    
    if(is.data.frame(HypGsa_Botres)){
      HypGsa_Botres$PC =p
      HypGsa_Botres$Loading_Direction ="Bottom1pcnt"
      df2bind<-rbind(df2bind,HypGsa_Botres)
      
    }
    
    PCA_loadings_HyperGSA.df<-rbind(PCA_loadings_HyperGSA.df,df2bind)
    Genes_top10_loadings_df<-rbind(Genes_top5_loadings_df,gene_df)
  }
  
  setwd(FinalOutCSV_perele_dir)
  write.csv(PCA_loadings_HyperGSA.df,paste0("Hyper GSA on PC loadings AllData Proteins as dimensions",dftype," dataset ",scaletype,".csv"),row.names=F)
  
  setwd(HyperGSA_onPCA_elewise_dir)
  
  if(nrow(PCA_loadings_HyperGSA.df) !=0){

  gsetypes <- unique(PCA_loadings_HyperGSA.df$Gset.Type)
  
  for(gs in 1:length(gsetypes)){
    
    df= filter(PCA_loadings_HyperGSA.df,
               Gset.Type == gsetypes[[gs]])%>%
      mutate(Loading_Direction = ifelse(Loading_Direction=="Bottom1pcnt","negative","positive"))
    
    if(is.data.frame(df)){
      pdf(paste0("Hyper GSA on PC loadings ",gsetypes[[gs]],dftype," dataset ",scaletype,".pdf"),width=20,height=30)
      
      print(
        ggplot(df,
               aes(x=paste("PC",PC),
                   y=Gset.Term.Enriched,
                   fill=Loading_Direction))+
          geom_tile()+
          #scale_fill_manual(values=c("darkblue","maroon"))+
          scale_fill_viridis_d(option="A",end=0.6,begin=0.2)+
          theme_SKA()+
          theme(axis.text.x = element_text(angle=90))+
          labs(x="Principal Component",title = gsetypes[[gs]])
      )
      dev.off()
    }
  }
  }
  
  ######################################
  ### Plot PCA results per elements  ###
  ######################################
  
  setwd(perele_hTMP_pca_umap_dir)
  
  pdf(paste0("PCA Alldata Proteins as Dimensions",dftype, " dataset ", scaletype,".pdf"),width=20,height=25)
  print(
    (scree+pp)/(pc1_vs_pc2+plot_spacer()+pc2_vs_pc3)/(pc1_vs_pc3+plot_spacer()+plot_spacer())/plot_spacer()
  )
  print( 
    loadingPlot/plot_spacer()
  )
  dev.off()
  
  ###############################
  ### Plot loadings heatmap  ###
  ###############################
  
  pca_loading_forhtmp <- pca_loading_df%>%
    reshape2::dcast(Genes ~ PC,value.var = "PC loading")
  
  rownames(pca_loading_forhtmp)<-pca_loading_forhtmp$Genes
  pca_loading_forhtmp<-pca_loading_forhtmp[,-1]
  
  # row clustering 
  dendr <- hclust(dist(pca_loading_forhtmp))
  dendr <- color_branches(dendr, k = 20)
  
  ## Colour for heatmap cells
  col_fun = colorRamp2(c(range(pca_loading_forhtmp)[1],0,range(pca_loading_forhtmp)[2]), c("blue","white","red"))
  
  
  pdf(paste0("PCA loadings heatmap",dftype," dataset ",scaletype,".pdf"),width=20,height=40)
  draw(ComplexHeatmap::Heatmap(as.matrix(pca_loading_forhtmp),
                               col=col_fun,
                               
                               cluster_rows = dendr,
                               use_raster = T,
                               raster_device = "CairoPNG",
                              
                               
                               row_split = 20,
                               row_gap = unit(3, "mm"),
                               row_names_gp = gpar(fontsize=2),   
                               
                               column_names_gp = gpar(fontsize = 11),
                               column_gap = unit(2, "mm"),
                               column_names_side = "top",
                               
                               heatmap_legend_param = list(title = "PCA loadings"),
                               right_annotation = row_anno
  ),
  merge_legend = TRUE
  
  )
  
  dev.off()
  
  
}

#################
### Set Paths ###
#################

projdir <- "/camp/lab/ralserm/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020"

proteomics_dir = "/camp/home/aulakhs/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020/Proteomics"
protRscript_dir <- paste0(proteomics_dir,"/RScripts")

FinalOutCSV_perele_dir <-paste0(proteomics_dir,"/Final Output Stat Imp/Per Element")

htmp_pca_umap_plotdir = paste0(proteomics_dir,"/Plots/Heatmaps_PCA_UMAP")

perele_hTMP_pca_umap_dir =paste0(htmp_pca_umap_plotdir,"/PerElement")
dir.create(perele_hTMP_pca_umap_dir,recursive = T)

HyperGSA_onPCA_dir=paste0(htmp_pca_umap_plotdir,"/HyperGSAonLoadings")
HyperGSA_onPCA_elewise_dir= paste0(HyperGSA_onPCA_dir,"/PerElement")
dir.create(HyperGSA_onPCA_elewise_dir,recursive = T)

############################################################
### Get Functions and variables from EDP Analysis script ###
############################################################

source(paste0(protRscript_dir,'/EDP_Analysis_Functions.R'))
source(paste0(projdir,"/Common Functions/Graphic_Parameters.R"))

##################################
### Make annotation dataframes ###
##################################

## Filter GOslim for Mitochondria, Golgi, Nucleus and ER

GO_CC_Anno<- filter(GO_gset_CC,  grepl("mitochon",term)| 
                      grepl("ER",term) | grepl("endoplasmic",term) |
                      grepl("nucleo",term) | grepl("nuclea",term)|
                      grepl("perox",term) | grepl("Golgi",term) |
                      grepl("vacuol",term)) %>%
  mutate(term=ifelse( grepl("mitochon",term),"Mitochondria",
                      ifelse(grepl("ER",term) | grepl("endoplasmic",term), "Endoplasmic Reticulum",
                             ifelse(grepl("nucleo",term) | grepl("nuclea",term),"Nucleus",
                                    ifelse(grepl("perox",term),"Peroxisome",
                                           ifelse(grepl("Golgi",term),"Golgi",
                                                  ifelse( grepl("vacuol",term),"Vacuole","None")))))))%>% 
  unique()%>%
  reshape2::dcast(ORF~term)

# Metal related terms from GO

GO_MFmetal_Anno <- filter(GO_gset_MF,  grepl("metal",term))%>%
  mutate(term="Metal related GO MF")%>%
  unique()%>%
  reshape2::dcast(ORF~term)

# Metal PDB

MetPDBanno<- filter(MetPDB,term %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  reshape2::dcast(ORF~term)

## Metabolic Process

GO_BP_metabolic_process_anno <- filter(GO_gset_BP,term =="metabolic process")%>%
  mutate(term="Metabolism")%>%
  unique()%>%
  reshape2::dcast(ORF~term)


## Make annotation data frame for rows

x = merge(GO_CC_Anno,GO_MFmetal_Anno,by="ORF",all=T)
x = merge(x,GO_BP_metabolic_process_anno,by="ORF",all=T)
x = merge(x,MetPDBanno,by="ORF",all=T)

x= merge(x,annostatdb,by="ORF",all=T)%>%
  mutate(Genes = as.character(lapply(ORF,convert_ORF2SingleGeneName)))



#################################################################
### Read in ESR & Growth Rate datasets for filtering heatmaps ###
#################################################################

setwd(dbdir)

StressResp_Brauer2007 <- read.csv("Brauer2007_SupplementaryTable1.csv",stringsAsFactors = F)

ESRgenes <- filter(StressResp_Brauer2007, ESR %in% c("up","down"))%>%
  pull(ORF)%>%
  unique()

ESRgenes_up <- filter(StressResp_Brauer2007, ESR == "up")%>%
  pull(ORF)%>%
  unique()

ESRgenes_down <- filter(StressResp_Brauer2007, ESR == "down")%>%
  pull(ORF)%>%
  unique()

GrowthRategenes <- filter(StressResp_Brauer2007, Growth.Rate.Response..1.5SD %in% c("up","down"))%>%
  pull(ORF)%>%
  unique()

GrowthRategenes_up <- filter(StressResp_Brauer2007, Growth.Rate.Response..1.5SD == "up")%>%
  pull(ORF)%>%
  unique()

GrowthRategenes_down <- filter(StressResp_Brauer2007, Growth.Rate.Response..1.5SD == "down")%>%
  pull(ORF)%>%
  unique()

## filter out all growth rate , ESR and TMR genes

genes2filtout <- unique(c(GrowthRategenes, ESRgenes))

genes2filtout<- as.character(lapply(genes2filtout, convert_ORF2SingleGeneName))

######################################
### Input NAfree data element wise ###
######################################


eles =c("Ca","Cu","Fe","Mg","Mn","Mo","Na","Zn")

for(e in 1:length(eles)){
  setwd(FinalOutCSV_perele_dir)
  
  EleProtData<-read.csv(paste0(eles[e]," Protein Quantities wImpValues FullMatrix.csv"),stringsAsFactors = F)
  
  EleProtData_sum <- EleProtData%>%
    group_by(BioSpecID,Protein.Ids,Element,Element.Concentration)%>%
    summarize(Mean.Protein.Quantity = log2(mean(2^(Log2.Protein.Quantity.),na.rm=T)))%>%
    mutate(Genes = as.character(lapply(Protein.Ids,convert_Uniprot2SingleGeneName)))
  
  EleProtData_mat <- EleProtData_sum%>%
    reshape2::dcast(Genes ~ BioSpecID, value.var =  "Mean.Protein.Quantity")
  rownames(EleProtData_mat)<-EleProtData_mat$Genes
  
  EleProtData_mat <- EleProtData_mat[,-1]
    
  #########################################
  ### Scale data with different methods ###
  #########################################

  setwd(FinalOutCSV_perele_dir)
  
  ### Scale with FC vs AllEle ###
  
  # Gene Names
    EleProtData_mat_FCAE <- EleProtData_mat - EleProtData_mat$`AllEle 1`
    EleProtData_mat_FCAE <- EleProtData_mat_FCAE[,-which(colnames(EleProtData_mat_FCAE)=="AllEle 1")]
  
  ### Scale with standard scalar centered around AllEle ###
    
    ## Read in stats from allele replicates only calculated with NA free dataframe
   
    AE_stats<- read.csv(paste0(eles[e]," AllEle Replicates Protein Quantity Stats FullMatrix.csv"),stringsAsFactors = F)[,c("Protein.Ids","SD.Protein.Quant",
                                                                                                                                        "Mean.Protein.Quant")]%>%
      mutate(Genes = as.character(lapply(Protein.Ids,convert_Uniprot2SingleGeneName)))
    
  ## Read in stats from all samples calculated with NA free dataframe
    AS_stats <- read.csv(paste0(eles[e]," AllSample Protein Quantity Stats FullMatrix.csv"),stringsAsFactors = F)[,c("Protein.Ids","SD.Protein.Quant",
                                                                                       "Mean.Protein.Quant")]%>%
    mutate(Genes = as.character(lapply(Protein.Ids,convert_Uniprot2SingleGeneName)))
  
    # Genes names as ids
    EleProtData_mat_ZscVsAE <- scale_rows_vsMeanAEsdAS(EleProtData_mat,"AllEle 1",AS_stats,idtype = "Genes")
    EleProtData_mat_ZscVsAE <- EleProtData_mat_ZscVsAE[,-which(colnames(EleProtData_mat_ZscVsAE)=="AllEle 1")]
  
    ## Zscore vs AE and AE stats  ( sd )used to calculte Z score
    
    EleProtData_mat_ZscVsAE_AEstats <- scale_rows_vsMeanAEsdAS(EleProtData_mat,"AllEle 1",AE_stats,idtype = "Genes")
    EleProtData_mat_ZscVsAE_AEstats <- EleProtData_mat_ZscVsAE[,-which(colnames(EleProtData_mat_ZscVsAE_AEstats)=="AllEle 1")]
    
    
  ### Scale with standard scalar centered mean of all samples ### stats calculated on whole with NA df
  
    # Genes names as ids
    EleProtData_mat_ZscVsAS <- scale_rows_vsMeanASsdAS(EleProtData_mat,AS_stats,idtype = "Genes")
  
  
  ##########################################################################
  ### Heatmaps PCA and UMAP for each type of scaling ### Elewise Dataset ###
  ##########################################################################
  
  
  scdfs<- list(EleProtData_mat_FCAE,
               EleProtData_mat_ZscVsAE,
               EleProtData_mat_ZscVsAE_AEstats,
               EleProtData_mat_ZscVsAS
               )
  
  scaletypes <- list("Log2FC",
                     "Zscore_Vs_AE",
                     "Zscore_Vs_AE_AEstats",
                     "Zscore_Vs_AS"
                     )
  
  setwd(perele_hTMP_pca_umap_dir)
  
  for(scdf in 1:length(scdfs)){
    
    getplot_heatmap_PCA_UMAP_perele(scaled_df = scdfs[[scdf]],
                                    scaletype = scaletypes[[scdf]], 
                                    dftype = eles[e])
    
  }
  
  
  scdfs_filt<-list()
  
  for(scdf in 1:length(scdfs)){
    
    scdfs_filt[[scdf]] = scdfs[[scdf]][-which(rownames(scdfs[[scdf]]) %in% genes2filtout ),]
    
  }
  
  for(scdf in 1:length(scdfs_filt)){
    
    getplot_heatmap_PCA_UMAP_perele(scdfs_filt[[scdf]],
                                    scaletype =scaletypes[[scdf]], 
                                    dftype = paste0(eles[e]," Growth_ESR free"))
    
  }

} # Close element loop
