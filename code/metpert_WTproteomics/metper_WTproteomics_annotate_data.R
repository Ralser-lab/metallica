##`---
#`  Title: "metpertWTproteomics_annotate_data.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 June 2020
#`  Description: Script to annotate proteomics data

#################
### Set Paths ###
#################

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# proteomics specific functions 

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))

DIANN_dir = paste0(metpert_WTproteomics_dir,"/DIANN_output")
proteomics_ODatharvest_dir = paste0(metpert_WTproteomics_dir,"/OD600_at_harvest")

plot_dir = paste0(metpert_WTproteomics_dir,"/output/plots")
qc_plot_dir = paste0(plot_dir,"/qc_plots")

dir.create(qc_plot_dir,recursive = T)
dir.create(paste0(qc_plot_dir,"/alldata_stats"))

intm_data_dir = paste0(paste0(metpert_WTproteomics_dir,"/output/intermediate_datafiles"))
dir.create(intm_data_dir,recursive = T)


######################################
### Read in stats file from DIA-NN ###
######################################

DIANNout_stats <- read.delim(paste0(DIANN_dir,"/metpertWTproteomics_rawdata_DIANNout.stats.tsv"),stringsAsFactors = F)%>%
                  filter(!File.Name %in% c(
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_5_C4.wiff", # based on intensity profile
                    # Based on media measurements
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_B6.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_B9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_C6.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_C9.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_D6.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_D8.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_D9.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_E6.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_E8.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_E9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_E12.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_F9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_F12.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_G5.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_G9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_6_G12.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_20210210_EDP_LO_7_24hr_D9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_20210210_EDP_LO_7_24hr_E9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_20210210_EDP_LO_7_24hr_F9.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_20210210_EDP_LO_7_24hr_F11.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_20210210_EDP_LO_7_24hr_G11.wiff",
                    
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_3_F6.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_3_G7.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_3_G2.wiff",
                    "Z:\\SimranQTOF2\\DatafromQTOF2\\EDP\\Samples\\SKA_202008_EDP_LO_5_C5.wiff"))

##############################################################################################################################
### Filter out files from LO8 ### We decided to not use chelator data and the controls in this plate have too few proteins ###
##############################################################################################################################

setwd(paste0(qc_plot_dir,"/alldata_stats"))
# Plot QC from raw DIANN data
plot_DIANNstat_distrb(DIANNout_stats,"DIANN output")

####################################################
### Filter out files based on DIANN stats report ###
####################################################

DIANNout_stats_filt <- filter(DIANNout_stats,
                              ## Filter out files which have total number of precursors < 60% of max number of precursors ###
                              Precursors.Identified >  0.4*max(DIANNout_stats$Precursors.Identified)&
                                Proteins.Identified > 1000  &
                                Total.Quantity > 1000000 &
                                MS1.Signal > 1000000 &
                                MS2.Signal > 10000000 &
                                Normalisation.Instability < 0.5 
)
# Plot QC after filtering

plot_DIANNstat_distrb(DIANNout_stats_filt,"DIANN output stat filter")

files_filtered_out <- unique(DIANNout_stats$File.Name)[which(!unique(DIANNout_stats$File.Name) %in% unique(DIANNout_stats_filt$File.Name))]
# Store names of raw files to keep based on DIA-NN stat filer

files2keep_DIANNstatfilter <- DIANNout_stats_filt$File.Name


####################################
### Read in layout of experiments ##
####################################

layout <- read.csv(paste0(metadata_dir,"/proteomics/layouts_proteomics_allcombined.csv"),stringsAsFactors = F)

#######################
### Read in OD data ###
#######################

ODfiles = dir(proteomics_ODatharvest_dir)
ODfiles = ODfiles[grep("24hr",ODfiles)]
ODfiles = ODfiles[!grepl("L8",ODfiles)]

OD_df <- vector()

for(i in 1:length(ODfiles)){
  
  
  odf = read_excel(paste0(proteomics_ODatharvest_dir,"/",ODfiles[i]))[24:31,]
  colnames(odf) = c("Row",1:12)  
  odf = reshape2::melt(odf,id.vars="Row",variable.name="Column",value.name="OD600")
  lo = gsub("L","",unlist(strsplit(ODfiles[i],"_"))[[4]])
  tp =  gsub("L","",unlist(strsplit(ODfiles[i],"_"))[[5]])
  odf$Plate_Position_TimePoint=paste(lo,paste0(odf$Row,odf$Column),tp,sep="_")
  
  OD_df = rbind(OD_df,odf)
}

OD_df$OD600 <- as.numeric(OD_df$OD600)

## Read in TECAN to spectrophotometer calibration file

Spec2Tecan = read.csv(paste0(acc_file_dir,"/od600_converter_spectrometer_to_NE3tecan.csv"),
                      stringsAsFactors = F) # Read in Spectrophotoeter to Tecan conversion file

Spec2Tecan = filter(Spec2Tecan,TECAN_NE3_OD < 0.8)

S2T_lm = lm(Spec2Tecan$SpecOD~Spec2Tecan$TECAN_NE3_OD)

ggplot(Spec2Tecan,aes(x=TECAN_NE3_OD,y=SpecOD))+
  geom_point(size=2,alpha=0.8)+
  geom_abline(slope = S2T_lm$coefficients[[2]],intercept = S2T_lm$coefficients[[1]],colour="red")+
  theme_minimal(base_size = 25)

OD_df <- mutate(OD_df,SpecOD=4*OD600*S2T_lm$coefficients[[2]]+S2T_lm$coefficients[[1]])

##################################
### Read-in DIA-NN output file ###
##################################

mpWTp_data <- as.data.frame(diann_load(paste0(DIANN_dir,"/metpertWTproteomics_rawdata_DIANNout.tsv")))

head(mpWTp_data)

mpWTp_data <- mpWTp_data %>%
              filter(File.Name %in% files2keep_DIANNstatfilter)%>%
              ## Filter to keep only proteotypic peptides
              filter(Proteotypic==1)%>%
              ## Filter on the basis of Q Value GG Q Value PG Q Value 
              filter(Q.Value <=0.01 & GG.Q.Value <=0.01 & PG.Q.Value <= 0.01)%>%
              dplyr::select(File.Name,Run,Protein.Ids,Genes,Protein.Group,
                            Precursor.Id,Precursor.Normalised,Q.Value,PG.Q.Value,Protein.Q.Value,GG.Q.Value)%>%
              unique()%>%
              mutate(Run=gsub("7_24hr","7-24hr",Run),
                     File.Name=gsub("[\\]","", File.Name),
                     File.Name=gsub("Z:SimranQTOF2DatafromQTOF2EDPSamples","" ,File.Name))%>%
              separate(Run, into = c(NA,NA,NA,NA,"Layout","Position"),sep="_" )%>%
              separate(Layout, into = c("Layout","TimePoint"),sep="-")%>%
              mutate(TimePoint = ifelse(is.na(TimePoint),"24hr",TimePoint))%>%
              mutate(Plate_Position = paste(Layout,Position,sep="_"))%>%
              mutate(Plate_Position_TimePoint =paste(Plate_Position,TimePoint,sep="_"))


## Merge with OD data

mpWTp_data <- merge(mpWTp_data,OD_df[,c(4,5)],by="Plate_Position_TimePoint")

## Merge with layout

mpWTp_data$Sample.Name = layout[match(mpWTp_data$Plate_Position,layout$Plate_Position),"Condition.or.Strain"]

## Fix sample names

mpWTp_data <- mpWTp_data%>%
              separate(Sample.Name,into = c("Element","Element.Concentration"),sep="_")%>%
              mutate(Element.Concentration = ifelse(grepl("AllEle",Element),1,Element.Concentration),
                     Element.Concentration = as.numeric(Element.Concentration))%>%
              mutate(ODunits_sampled = ifelse(Layout %in% c(7,8), SpecOD*2*1.6,SpecOD*1.6))%>%
              mutate(LO_TP=paste(Layout,TimePoint),
                     BioSpecID = paste(Element,Element.Concentration),
                     Sample.Type= ifelse(Element == "AllEleprepQC","Extraction Control",
                                         ifelse(grepl("AllEle",Element) ,
                                                "Control","Sample")))%>%
              ## Add digestion batch
              mutate(Digestion_Batch = ifelse(LO_TP %in% c("1 24hr","2 24hr","3 24hr","4 24hr"),"DB_1",
                                              ifelse(LO_TP %in% c("5 24hr","6 24hr","7 24hr","8 24hr"),"DB_2",
                                                     "DB_3" )))
######################
### Plot OD values ###
######################

OD_2p <- mpWTp_data%>%
  filter(!Element %in% c("AllEleprepQC","Empty"))%>%
  dplyr::select(Plate_Position,Element,TimePoint,Element.Concentration,
                SpecOD,ODunits_sampled)%>%
  unique()

ggplot(OD_2p, aes(x=SpecOD, colour=TimePoint))+
  geom_density(aes(y=..scaled..))+
  facet_wrap("Element")+
  theme_minimal()


ggplot(OD_2p, aes(x=ODunits_sampled, colour=TimePoint))+
  geom_density(aes(y=..scaled..))+
  facet_wrap("Element")+
  theme_minimal()

################################################################################################
### Filter samples based on OD units sampled & time of sampling - leaving out 36 hrs for now ###
################################################################################################

## NOTE : 

mpWTp_data <- mpWTp_data%>%
  mutate(ODunits_sampled = ifelse(Element =="AllEleprepQC",2,ODunits_sampled))%>%
  # -- Make sure to not filter out AllEleprepQC samples as they were added right before sample prep and don't have OD measurements
  filter(ODunits_sampled > 0.75 & Element != "Empty")%>%
  # Filter out AllEleNew samples that were on the edges, convert all other AllEleNew to AllEle
  filter(!(Element =="AllEleNew" & Position %in%c("C1","D1","E1")))%>%
  filter(Element != "B")%>%
  mutate(BioSpecID = gsub("AllEleNew","AllEle",BioSpecID),
         Element = gsub("AllEleNew","AllEle",Element) )%>%
  filter(TimePoint == "24hr" & !BioSpecID %in% c("Cu 50","Cu 100",
                                                             "Fe 100",
                                                             "K 10",
                                                             "Mg 20","Mg 50",
                                                             "Mo 20","Mo 50",
                                                             "Zn 2e-04","Zn 0.001","Zn 0.002") )
# Filter out some excess and depletion conditions 

#filter(Element != "B")


### Filter DIANN stats report for files remaining after filtering for OD units and replot stats

files_after_ODfilter <- unique(mpWTp_data$File.Name)

DIANNout_stats_filt<-DIANNout_stats_filt%>%
  mutate(File.Name=gsub("[\\]","", File.Name),
         File.Name=gsub("Z:SimranQTOF2DatafromQTOF2EDPSamples","" ,File.Name))

DIANNout_stats_OD_filt <- filter(DIANNout_stats_filt,File.Name %in% files_after_ODfilter)

setwd(paste0(qc_plot_dir,"/AllDataStats"))
plot_DIANNstat_distrb(DIANNout_stats_OD_filt,"stats and OD filter")


###################################################################
### Read-in Acquisition sequence file & merge with DIANN output ###
###################################################################


AcSeq <- read.delim(paste0(metadata_dir,"/proteomics/acquistion_sequence.txt"),sep="\t",stringsAsFactors = F)
colnames(AcSeq)<-"file_metdat"
AcSeq <- AcSeq%>%
         mutate(file_metdat = stringr::str_squish(file_metdat))%>%
         separate(file_metdat, into = c(NA,NA,NA,NA,NA,"Date","Time","File.Name"),remove = F,sep=" ")

# Add Acquisition Data and time to main Data Frame : 
mpWTp_data <- merge(mpWTp_data,AcSeq,by="File.Name")

##########################################################################
### Write unfiltered data for detecting proteins that switch on / off  ###
##########################################################################


write.csv(mpWTp_data, paste0(intm_data_dir,"/AllData_prepared_noPrecFilter.csv"), row.names = F)

####################################
### Make Experiment Design table ###
####################################

ExpDesign <- mpWTp_data%>%
              mutate(condition = trimws(gsub("NA","",paste(Element,Element.Concentration))))%>%
              dplyr::select(File.Name,condition)%>%
              unique()%>%
              mutate(label=File.Name)%>%
              mutate(number=1)%>%
              group_by(condition)%>%
              mutate(replicate = cumsum(number))%>%
              ungroup()


ExpDesign_AE<- mpWTp_data%>%
  filter(Sample.Type != "Sample")%>%
  mutate(condition = trimws(gsub("NA","",paste(Sample.Type,LO_TP))))%>%
  dplyr::select(File.Name,condition)%>%
  unique()%>%
  mutate(label=File.Name)%>%
  mutate(number=1)%>%
  group_by(condition)%>%
  mutate(replicate = cumsum(number))%>%
  ungroup()


write.csv(ExpDesign[,c("label","condition","replicate")],
          paste0(intm_data_dir,"/EDP_ExperimentDesign.csv"),row.names = F)

write.csv(ExpDesign_AE[,c("label","condition","replicate")],
          paste0(intm_data_dir,"/EDP_ExperimentDesign_AEonly.csv"),row.names = F)

