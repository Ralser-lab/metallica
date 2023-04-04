#`---
#`  Title: "metallomics_2_process_od600atsamplingdata.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 April 2021 
#`  Description: Script to process OD600 measurements of S.cerevisiae cultures recorded before sampling for intracellular metal quantification
#`---

#############################################
### source paths functions and libraries  ###
#############################################

source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/metpert_WTmetallomics/metallomics_0_libraries_functions.R"))

#################
### set paths ###
#################

# script specific paths 
od_dir <- paste0(metpert_WTmetallomics_dir,"/od600_at_sampling")

#######################
### Read in OD data ###
#######################

# Read in Spectrophotoeter to Tecan conversion file
Spec2Tecan <- read.csv(paste0(acc_file_dir,"/od600_converter_spectrometer_to_NE3tecan.csv"),stringsAsFactors = F) 
Spec2Tecan <- filter(Spec2Tecan,TECAN_NE3_OD < 1.2)

# make a linear regression model to convert between spectrophotometer and tecan od600 values

S2T_lm <- lm(Spec2Tecan$SpecOD~Spec2Tecan$TECAN_NE3_OD)

pdf(paste0(acc_file_dir,"/TECAN3NE_to_Spectrophotometer.pdf"),width=5,height=5)
ggplot(Spec2Tecan,aes(x=TECAN_NE3_OD,y=SpecOD))+
  geom_point(size=2,alpha=0.8)+
  geom_abline(slope = S2T_lm$coefficients[[2]],intercept = S2T_lm$coefficients[[1]],colour="red")+
  theme_metallica()
dev.off()

# Read in all csvs with OD600 data 

od_csvs <- dir(od_dir)[grep("metallomics_OD600as_",dir(od_dir))]

od600_allplates = data.frame()

for(i in 1:length(od_csvs)){
  
  f <- read.csv(paste0(od_dir,"/",od_csvs[i]),stringsAsFactors = F)[24:31,1:13]
  colnames(f) <- c("Row",1:12)
  rownames(f) <- f$Row
  
  f <- reshape2::melt(f, id.vars = "Row",variable.name = "Col",value.name="OD600")
  f$PlateID <- unlist(strsplit(unlist(strsplit(od_csvs[i],"[.]"))[[1]],"_OD600as_"))[[2]]
  od600_allplates <-  rbind (od600_allplates,f)
}

od600_allplates <- od600_allplates%>%
                   mutate(Position=paste0(Row,Col))%>%
                   separate(PlateID,into=c("Layout","Dilution" ),remove=F)%>%
                   mutate(lop=paste(Layout,Position,sep="_"))%>%
                   mutate(lop=gsub("LO","",lop))%>%
                   mutate(OD600dil=as.numeric(OD600)*as.numeric(Dilution))

### Check dilutions ###
pdf(paste0(od_dir,"/dilutioncheck_od600_before_metallomics.pdf"),width = 10, height = 6)
ggplot(filter(od600_allplates,Position %in% c("D2","D5","B4","B8","G6","F10")),
       aes(x=PlateID,y=as.numeric(OD600),
           colour=Dilution))+
  geom_point()+
  theme_metallica()+
  labs(title = "OD600")

ggplot(filter(od600_allplates,Position %in% c("D2","D5","B4","B8","G6","F10")),
       aes(x=PlateID,y=as.numeric(OD600dil),
           colour=Dilution))+
  geom_point()+
  theme_metallica()+
  labs(title = "OD600 normalised to dilution")
dev.off()

od600_allplates <- od600_allplates%>%
                    dplyr::select(-c(Row,Col,Layout,Position,Dilution,PlateID))%>%
                    group_by(lop)%>%
                    summarize(OD600corrected = mean(OD600dil,na.rm = T))%>%
                    ungroup()%>%
                    mutate(SpecOD=OD600corrected*S2T_lm$coefficients[[2]]+S2T_lm$coefficients[[1]])

write.csv(od600_allplates,paste0(od_dir,"/OD_at_sampling_Metallomics_ATP.csv"),row.names = F)

