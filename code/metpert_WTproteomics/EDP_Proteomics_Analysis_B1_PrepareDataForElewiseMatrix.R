#######################################################################
### EDP Script - 1B - Data preparation for complete matrix analysis ###
#######################################################################

#################
### Set Paths ###
#################

projdir <- "/camp/lab/ralserm/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020"
proteomics_dir = "/camp/home/aulakhs/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020/Proteomics"

protRscript_dir <- paste0(proteomics_dir,"/RScripts")

source(paste0(projdir,"/Common Functions/Initialize_Paths.R"))
source(paste0(projdir,"/Common Functions/Graphic_Parameters.R"))
source(paste0(protRscript_dir,'/EDP_Analysis_Functions.R'))

datdir = paste0(proteomics_dir,"/DIANNout")
plotdir=paste0(proteomics_dir,"/Plots")
QCplotdir=paste0(plotdir,"/QC_plots")
Int_datdir =paste0(paste0(proteomics_dir,"/IntermediateDFs"))

FinalOutCSV_perele_dir <-paste0(proteomics_dir,"/Final Output Stat Imp/Per Element")
dir.create(FinalOutCSV_perele_dir,recursive = T)

#################################################################################################################
#################################################################################################################
#################################################################################################################
########################################### All samples - per Element ########################################### 
#################################################################################################################
#################################################################################################################
#################################################################################################################

PEle_QCplot_forFullMatrix_dir = paste0(QCplotdir,"/forEleWiseFullMatrix/All Samples")

dir.create(PEle_QCplot_forFullMatrix_dir,recursive = T)

###################
### Input Data ###
###################

## Read in uniprot - sgd - gene name conversion file to add gene names ###

GenProt_SGD <- read.csv(paste0(dbdir,"SGD_Uniprot.csv"),stringsAsFactors = F)

#######################################
### Read in Experiment Design files ###
#######################################

setwd(Int_datdir)

Exp_Design<-read.csv("EDP_ExperimentDesign.csv",stringsAsFactors = F)

Exp_Design_AEonly<-read.csv("EDP_ExperimentDesign_AEonly.csv",stringsAsFactors = F)

###############################
### Read in All Sample data ###
###############################

Sample_data<-read.csv("AllData_prepared_noPrecFilter.csv",stringsAsFactors = F)

## Loop through each element and carryout all analysis per element 

eles = unique(Sample_data$Element)[-grep("AllEle",unique(Sample_data$Element))]

numprot_allmetser <- vector()

for( e in 1:length(eles)){
  
  Sample_data_ele <- filter(Sample_data,Element %in% c(eles[e],"AllEle","AllEleprepQC") )
  
  ele_los <- unique(filter(Sample_data_ele,Element ==eles[e])$Layout)
  
  Sample_data_ele <- filter(Sample_data_ele, Layout %in% ele_los)
  
  ###########################################################################
  ### Filter out Precursors that have InterBatch CVs > 30 % in AE samples ###
  ###########################################################################
  
  ## Compute CVs Intra Batch and Inter Batch CVs
  
  AE_CVs_preBC <- Sample_data_ele%>%
    filter(Sample.Type == "Control")%>%
    mutate(Lo_Prot_PrecID = paste(Layout, Protein.Ids,Precursor.Id,sep="_"),
           Prot_PrecID=paste(Protein.Ids,Precursor.Id,sep="_"))%>%
    group_by(Lo_Prot_PrecID)%>%
    mutate(CV.AE.Intra.Batch = cv(Precursor.Normalised,na.rm = T))%>%
    ungroup()%>%
    group_by(Prot_PrecID)%>%
    mutate(CV.AE.Inter.Batch = cv(Precursor.Normalised, na.rm = T))%>%
    mutate(QCtype="Pre Batch Correction")
  
  PrecIDList<-unique(filter(AE_CVs_preBC,CV.AE.Intra.Batch <= 30)$Prot_PrecID)
  
  Sample_data_ele <- Sample_data_ele%>% 
    mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
    filter(Prot_PrecID %in% PrecIDList)
  
  ### Compute Protein Group quantities ###
  
  PG_quants<-get_PGQ_df(Sample_data_ele)
  
  ## Merge PG_qunats with annotation
  
  PG_quants<-merge(PG_quants,
                   unique(Sample_data_ele[,c("File.Name","Layout","Element",
                                             "Element.Concentration",
                                             "Digestion_Batch",  "TimePoint",
                                             "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                   by="File.Name")
  
  
  ######################################################################
  ### Filter out proteins that are detected in < 50 % of the samples ###
  ######################################################################
  
  protcompletness_filter <-PG_quants %>%
    filter(!grepl("AllEle",BioSpecID ))%>%
    dplyr::select(Protein.Ids,File.Name,Protein.Quantity)%>%
    reshape2::dcast(Protein.Ids~File.Name)%>%
    reshape2::melt()%>%
    group_by(Protein.Ids)%>%
    summarize(Fraction_not_NA = 1- (sum(is.na(value))/length(value)) )
  
  Prot2keep <- unique(filter(protcompletness_filter, Fraction_not_NA > 0.60 )$Protein.Ids)
  
  ## Remove samples with < 1300 proteins quantified 
  
  numprot_filter <- PG_quants %>%
    dplyr::select(File.Name,Protein.Ids,Protein.Quantity)%>%
    na.omit()%>%
    group_by(File.Name)%>%
    summarize(numprot_detected = length(Protein.Ids))%>%
    ungroup()%>%
    filter(numprot_detected >=1300)%>%
    pull(File.Name)
  
  PG_quants <- filter(PG_quants, File.Name %in% numprot_filter)      
  PG_quants <- filter(PG_quants, Protein.Ids %in% Prot2keep)
  
  Exp_Design <- read.csv(paste0(Int_datdir,"/EDP_ExperimentDesign.csv"),stringsAsFactors = F)
  
  Exp_Design <- filter( Exp_Design, label %in% numprot_filter)
  
  ## Setup annotations for proBatch
  
  setup_proBatch_annos(PG_quants,onlyAEdat = F, perele = T)
  essential_columns <- c(PGcol, PGintcol, fname_col)
  
  
  
  ### Prepare data for use with DEP ###
  
  DEP_PGquants <- reshape2::dcast ( PG_quants,
                                    Protein.Ids~File.Name,
                                    value.var = "Protein.Quantity")%>%
    mutate ( name = Protein.Ids,
             ID = Protein.Ids)
  
  ### Create Summarized Experiment Object ##
  
  ## Note : this will log2 transform all values
  
  SumExp <- make_se(DEP_PGquants,seq(2,ncol(DEP_PGquants)-2),Exp_Design )
  
  ############################################
  ### Batch Correct at the precursor level ###
  ############################################
  
  Samdata_fMedCorr <- AEmed_correction_Precursor(Sample_data_ele)
  
  ### Convert batch corrected data into protein quantities ###
  
  PG_quants_BC<-get_PGQ_df(Samdata_fMedCorr)
  
  ## Merge PQ_quants with annotation
  
  PG_quants_BC<-merge(PG_quants_BC,
                      unique(Sample_data_ele[,c("File.Name","Layout","Element","Element.Concentration",
                                                "Digestion_Batch","TimePoint",
                                                "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                      by="File.Name")
  
  
  PG_quants_BC <- filter(PG_quants_BC, File.Name %in% numprot_filter)      
  PG_quants_BC <- filter(PG_quants_BC, Protein.Ids %in% Prot2keep)
  
  DEP_PGquants_BC <- reshape2::dcast ( PG_quants_BC,
                                       Protein.Ids~File.Name,
                                       value.var = "Protein.Quantity")%>%
    mutate ( name = Protein.Ids,
             ID = Protein.Ids)
  
  SumExp_BatchCorr <- make_se(DEP_PGquants_BC,seq(2,ncol(DEP_PGquants_BC)-2),Exp_Design)  
  
  ######################
  ### Median Scaling ###
  ######################
  
  PG_quants_BCMS <- PG_quants_BC%>%
    group_by(File.Name)%>%
    mutate( MPQ_perFile= median(Protein.Quantity,na.rm=T))%>%
    ungroup()%>%
    # Calculate mean of median protein intensities across all files
    mutate( Median_MPQ = median(MPQ_perFile,na.rm=T),
            # Scale each protein quantity such that all files will have same Median 
            Median_Scaled_Protein_Quantity = (Protein.Quantity/MPQ_perFile)*Median_MPQ)
  
  
  DEP_PGquants_BCMS <- reshape2::dcast ( PG_quants_BCMS,
                                         Protein.Ids~File.Name,
                                         value.var = "Median_Scaled_Protein_Quantity")%>%
    mutate ( name = Protein.Ids,
             ID = Protein.Ids)
  
  SumExp_BatchCorrMedSc <- make_se(DEP_PGquants_BCMS,seq(2,ncol(DEP_PGquants_BCMS)-2),Exp_Design)  
  
  ##############################
  ### Fill in missing values ###
  ##############################
  
  for_NAfill <- PG_quants_BCMS %>%
    reshape2::dcast(Protein.Ids~File.Name, value.var = "Median_Scaled_Protein_Quantity")
  
  # If missing in all replicates of a condition --> fill in with minimum value measured in dataset
  # Else - replace with median of all measured replicates
  
  NAfree.df <- fill_in_NAs(for_NAfill)
  
  NAfree.df.fSE <- data.frame(NAfree.df)%>%
    mutate(ID = rownames(NAfree.df),
           name = rownames(NAfree.df))
  
  SumExp_noNA <- DEP::make_se(NAfree.df.fSE,seq(1:(ncol(NAfree.df.fSE)-2)),Exp_Design)
  
  ##################################
  ### Collect all SumExp objects ### 
  ##################################
  
  SEos<-list(
    SumExp,
    SumExp_BatchCorr,
    SumExp_BatchCorrMedSc,
    SumExp_noNA
    
  )
  
  names(SEos) <- c("Log2 Raw",
                   "BatchCorrected",
                   "BatchCorrected_medianScaled",
                   "BatchCorrected_medianScaled_noNA"
  )
  
  # Get DFs out of se object ( these will have log2 scale values )
  
  SEdfs<-lapply(SEos, SumExpObj_to_df)
  
  #################################################################
  ## Check intensity profiles to find samples with bad profiles ###
  #################################################################
  
  setwd(PEle_QCplot_forFullMatrix_dir)
  

  pdf(paste0("Protein Intensity profiles allsamp elewise",eles[e],".pdf"),width=20,height=15)
  for(i in 1:length(SEdfs)){
    plot_protint_dist(SEdfs[[i]], names(SEdfs)[[i]])
  }
  
  ##############################
  ### Make proBatch QC plots ###
  ##############################
  
  ### Set up Pro batch annotations
  
  
  for( i in 1:length(SEdfs)){
    
    if(length(unique(sample_anno$Digestion_Batch)) >1){
      
      sel_annos <- c("LO_TP","Sample.Type","Digestion_Batch") ## annotations to plot by
      
    }else if(length(unique(sample_anno$LO_TP)) > 1){
      
      sel_annos <- c("LO_TP","Sample.Type") ## annotations to plot by
      
    }else{
      
      sel_annos <- c("Sample.Type") ## annotations to plot by
      
    }
    
    make_proBatchQCplots(SEdfs[[i]],
                         sel_anno_hclust=sel_annos,
                         onlyAEdat = F,
                         dsname=names(SEdfs)[[i]],
                         AE_or_samp = "Samples",
                         sanno = sample_anno )
    
  }
  
  
  
  ######################################################
  ######################################################
  ### AE median Batch Corrected DFs to write to disk ###
  ######################################################
  ######################################################
  
  setwd(FinalOutCSV_perele_dir)
  
  Fin_wNA_df <- merge(PG_quants_BCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_BCMS)))],by="File.Name")%>%
    mutate(BioSpecID=trimws(BioSpecID),
           `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
    filter(BioSpecID != "AllEleprepQC 1")
  
  write.csv(Fin_wNA_df,paste0(eles[e]," Protein Quantities w NAs FullMatrix.csv"),row.names = F)
  
  
  # Calculate per Gene mean, sd, median, mad for AllEle replicates across dataset and write to CSV
  
  AERep_Stats <- Fin_wNA_df%>%
    filter(BioSpecID =="AllEle 1")%>%
    group_by(Protein.Ids)%>%
    mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
    ungroup()
  
  ## Plot distribution of fraction of controls that are NAs
  
  ggplot(unique(AERep_Stats[,c("Protein.Ids","Num.Non.NAs")]),
         aes(x=Num.Non.NAs))+
    geom_histogram(bins=25,colour="darkred",
                   fill="darkred")+
    theme_SKA()
  
  AERep_Stats<-AERep_Stats%>%
#    filter(Num.Non.NAs > 10)%>%
    group_by(Protein.Ids)%>%
    summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
              Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
    )%>%
    ungroup()%>%
    mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
    as.data.frame()
  
  rownames(AERep_Stats)<-AERep_Stats$Protein.Ids
  
  write.csv(AERep_Stats,paste0(eles[e]," AllEle Replicates Protein Quantity Stats FullMatrix.csv"),row.names=F)
  
  
  # Calculate per Gene mean, sd, median, mad for all samples AllEle replicates across dataset and write to CSV
  
  AllSample_Stats <- Fin_wNA_df%>%
    filter(!grepl("AllEle",BioSpecID))%>%
    group_by(Protein.Ids)%>%
    mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
    ungroup()
  
  ## Plot distribution of fraction of controls that are NAs
  
  ggplot(unique(AllSample_Stats[,c("Protein.Ids","Num.Non.NAs")]),
         aes(x=Num.Non.NAs))+
    geom_histogram(bins=25,colour="darkred",
                   fill="darkred")+
    theme_SKA()
  
  AllSample_Stats<-AllSample_Stats%>%
  #  filter(Num.Non.NAs >10)%>%
    group_by(Protein.Ids)%>%
    summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
              Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
    )%>%
    ungroup()%>%
    mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
    as.data.frame()
  
  rownames(AllSample_Stats)<-AllSample_Stats$Protein.Ids
  
  write.csv(AllSample_Stats,paste0(eles[e]," AllSample Protein Quantity Stats FullMatrix.csv"),row.names=F)
  
  #### Replicate stats
  
  Ele_RepStats <- Fin_wNA_df%>%
    filter(!grepl("AllEle",BioSpecID))%>%
    group_by(Protein.Ids,Element,BioSpecID)%>%
    summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
              Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
    )%>%
    mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
    ungroup()%>%
    group_by(Protein.Ids, Element)%>%
    summarize(Mean.COVreps = mean(Cov,na.rm=T))%>%
    ungroup()%>%
    as.data.frame()
  
  rownames(Ele_RepStats)<-paste(Ele_RepStats$Protein.Ids,Ele_RepStats$Element)
  
  write.csv(Ele_RepStats,paste0(eles[e]," Protein Quantity Replicate Stats.csv"),row.names=F)
  
  ###################################################
  ### Prepare data without NAs for writing to CSV ### AE median Batch Corrected
  ###################################################
  
  FinnoNA_df <- NAfree.df %>%
    data.frame()%>%
    mutate(Protein.Ids=rownames(NAfree.df))%>%
    reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
                   value.name="Protein.Quantity")%>%
    mutate(`Log2(Protein.Quantity)`=log2(Protein.Quantity))
  
  FinnoNA_df <- merge(FinnoNA_df,sample_anno,by="File.Name")%>%
    mutate(BioSpecID = trimws(BioSpecID))%>%
    filter(BioSpecID != "AllEleprepQC 1")
  
  write.csv(FinnoNA_df,paste0(eles[e]," Protein Quantities wImpValues FullMatrix.csv"),row.names = F)
  
  
  ### Write unimputed data for clustering ###
  
  FinwNA_df_forclus <- Fin_wNA_df%>%
    filter(!BioSpecID %in% c("Cu 50","Cu 100",
                             "Mo 10","Mo 20",
                             "Mo 25","Mo 50","Mo 100","Fe 100") 
    )
  
  
  write.csv(FinwNA_df_forclus,paste0( eles[e]," Protein Quantities noImp for clustering FullMatrix.csv"),row.names=F)
  
  numprot_BioSpecID <- FinwNA_df_forclus %>%
    na.omit()%>%
    dplyr::select(Protein.Ids,BioSpecID)%>%
    unique()%>%
    group_by(BioSpecID)%>%
    summarize(Num_Prot = length(Protein.Ids))

  numprot_BioSpecID$Element = eles[e]
  
  numprot_allmetser <- rbind(numprot_allmetser, numprot_BioSpecID)
}


setwd(QCplotdir)
pdf("Number of proteins per condition Branch B ylim.pdf",width = 12,height=5)

ggplot(numprot_allmetser,
       aes(x=BioSpecID,y=Num_Prot,
           colour=BioSpecID,
       ))+
  geom_point(size=3,alpha=1.0)+
  scale_color_manual(values=colkey_BioSpecID)+
  theme_SKA()+
  theme(
    panel.grid.major = element_line(colour="lightgray",size = 0.1),
    panel.border = element_rect(size=0.15),
    legend.position = "none",
    axis.text.x = element_text(angle=90,size=9.5,hjust = 0.9),
    axis.ticks = element_line(size=0.75))+
  labs(x="", y="Number of proteins quantified")+
  ylim(1300,2500)
dev.off()
