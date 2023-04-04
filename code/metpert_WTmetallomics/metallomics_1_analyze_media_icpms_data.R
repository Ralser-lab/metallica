#`---
#`  Title: "metallomics_1_analyze_media_icpms_data"
#`  @author: Simran Kaur Aulakh
#`  Date created: 8 April 2021 
#`  Description: Script to analyze all ICP-MS data from cultivation media
#`---

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# experiment specific
source(paste0(code_dir,"/metpert_WTmetallomics/metallomics_0_libraries_functions.R"))

#################
### set paths ###
#################

# inputs
icpms_data_dir <- paste0(metpert_WTmetallomics_dir,"/data_from_agilent_masshunter/media_metallomics")
lo_dir <- paste0(metpert_WTmetallomics_dir,"/layouts_96wellplates/media_icpms_layouts")

# outputs
plot_dir <- paste0(metpert_WTmetallomics_dir,"/results/plots/media_metallomics")
outputtables_dir <- paste0(metpert_WTmetallomics_dir,"/results/output_tables")

dir.create(plot_dir, recursive = T)
dir.create(outputtables_dir, recursive = T)

#################
### set paths ###
#################

media_icpms_data <- data.frame()

for( i in 1:4){
  
  pl <- read.csv(paste0(icpms_data_dir,"/data_media_icpms_from_AgilentMassHunter_plate",i,".csv") ,stringsAsFactors = F)
  lo <- read.csv(paste0(lo_dir,"/layout_media_icpms_plate",i,"_wide.csv"),stringsAsFactors = F)
  
  colnames(lo) <- c("row",1:12)
  
  plate_layout <- na.omit(melt(lo,id.vars = c("row"), value.name = "condition", variable.name = "column"))%>%
                  mutate(position = paste0(row, column))%>%
                  dplyr::select(-c(row,column))%>%
                  filter(condition != "Empty")

  pl <- pl%>%
        tidyr::separate(Sample.Name,into = c("rack","plate","position"),sep="-")%>%
        dplyr::select(-c("rack","plate"))
      
  colnames(pl) <- gsub("[.]{1}"," ",gsub("[.]{2}","_",gsub("X","",colnames(pl))))
  
  pl$media <- plate_layout[match(pl$position,plate_layout$position),"condition"]
  pl$plate <- i
  media_icpms_data <- rbind(media_icpms_data,pl)
}


media_icpms_data_mlt <- melt(media_icpms_data[,-1],id.vars = c("media","plate"), value.name = "CPS", variable.name = "metal")%>%
                        tidyr::separate(metal,into=c("atomic_mass","metal","gas_mode"),sep="_")%>%
                        tidyr::separate(metal,into=c("metal","ISTD"),sep=" ")%>%
                        tidyr::separate(media,into = c("element","concentration","chelator_concentration"),sep="_",remove = F)%>%
                        tidyr::separate(concentration,into=c("element_concentration","chelator"),sep="[+]")%>%
                        mutate(element_concentration = as.numeric(lapply(element_concentration,Conc2NumConc)),
                               chelator_concentration = as.numeric(lapply(chelator_concentration, ChelCon2Num)),
                               CPS = as.numeric(CPS))%>%
                        ## filter out media that was excluded from whole analysis because cells died or not a metal
                        filter(!media %in% c("Empty" , "K_0x","Mg_x/1000","Zn_x/5000") &
                                                     !grepl("B",media))

## Filter out media that were remade and remeasured (plate 4) from other plates

pl4med <- unique(filter(media_icpms_data_mlt,plate ==4)$media)
pl4med <- pl4med[grepl("Fe",pl4med)|grepl("Zn",pl4med)|grepl("Dfp",pl4med)|grepl("DiP",pl4med)]

media_icpms_data_mlt <- filter(media_icpms_data_mlt,! (media %in% pl4med & plate != 4))

eles=unique(media_icpms_data_mlt$element)


pdf(paste0(plot_dir,"/media_measurements_metal_wise_plots_unnormalised.pdf"),width = 20,height = 20)

for( i in 1:length(eles)) {
  
  df <- filter(media_icpms_data_mlt,element %in% eles[[i]] & is.na(chelator))
  df <- rbind(df,filter(media_icpms_data_mlt,element =="AllEle" & is.na(chelator)))
  
  print(ggplot(df,aes(x=log10(element_concentration),y=log10(CPS), colour=gas_mode))+
        geom_point(size=3)+
        geom_line(aes(group=gas_mode))+
        facet_wrap(c("metal"),ncol=4)+
        theme_metallica(base_size = 21)+
        theme( panel.grid.major = element_line(size=0.5,colour="lightgray"), 
               panel.grid.minor = element_line(size=0.5,colour="lightgray"))+
        scale_color_viridis_d(end=0.7)+
        geom_vline(xintercept = log10(1))+
        labs(title = paste0("metal perturbed",eles[i]),
                 x = "log10(concentration perturbed metal)",
                 y = "log10 (CPS of measured metal)")
        )
}
dev.off()


## Plot chelator media

chels <- as.character(na.omit(unique(media_icpms_data_mlt$chelator)))

pdf(paste0(plot_dir,"/media_measurements_chelator_wise_plots_unnormalised.pdf"),width = 20,height = 20)

for(c in 1:length(chels)){
  
  df=filter(media_icpms_data_mlt,chelator==chels[c])
  df=rbind(df,filter(media_icpms_data_mlt,element =="AllEle" & is.na(chelator)))
  
  print(
    ggplot(df,aes(x=chelator_concentration,y=log10(CPS), colour=element))+
    geom_point()+
    geom_line(aes(group=element))+
    facet_wrap(c("metal","gas_mode"),scales="free",ncol=5)+
    scale_color_viridis_d(end=0.7)+
    theme_metallica(base_size = 21)+
    theme( panel.grid.major = element_line(size=0.5,colour="lightgray"), 
           panel.grid.minor = element_line(size=0.5,colour="lightgray"))+
    geom_vline(xintercept = 0)+
    labs(title=paste0("chelator ",chels[c]),x= "log10(Concentration of perturbed metal)",
         y ="log10 (CPS of measured metal)")
  )
}


########################
### Normalize to  Sc ###
########################

for( i in 1:nrow(media_icpms_data_mlt)){
  
  media_icpms_data_mlt[i,"Sc.CPS"]=median(filter(media_icpms_data_mlt,media == media_icpms_data_mlt[i,"media"] & plate ==media_icpms_data_mlt[i,"plate"]&gas_mode==media_icpms_data_mlt[i,"gas_mode"] & metal == "Sc" )$CPS)
  
}

media_icpms_data_mlt$Sc.Normalized.CPS = media_icpms_data_mlt$CPS / media_icpms_data_mlt$Sc.CPS

pdf(paste0(plot_dir,"/media_measurements_metal_wise_plots_Scnormalised.pdf"),width = 20,height = 20)

for( i in 1:length(eles)) {

      df <- filter(media_icpms_data_mlt,element %in% eles[[i]] & is.na(chelator))
      df <- rbind(df,filter(media_icpms_data_mlt,element =="AllEle" & is.na(chelator)))
    
  print(
      ggplot(df,aes(x=log10(element_concentration),y=log10(Sc.Normalized.CPS), colour=gas_mode))+
      geom_point(size=3)+
      geom_line(aes(group=gas_mode))+
      facet_wrap(c("metal"),ncol=4)+
      theme_metallica(base_size = 21)+
      theme( panel.grid.major = element_line(size=0.5,colour="lightgray"), 
             panel.grid.minor = element_line(size=0.5,colour="lightgray"))+
      scale_color_viridis_d(end=0.7)+
      geom_vline(xintercept = log10(1))+
      labs(title=paste0("perturbed metal ",eles[i]),
           x= "log10 (Concentration of perturbed metal)",
           y ="Sc Normalized CPS")
      )
  }

## Plot chelator media

chels <- as.character(na.omit(unique(media_icpms_data_mlt$chelator)))

pdf(paste0(plot_dir,"/media_measurements_chelator_wise_plots_Scnormalised.pdf"),width = 20,height = 20)

for(c in 1:length(chels)){
  if(chels[c] != "DMSO"){
    
    df <- filter(media_icpms_data_mlt,chelator==chels[c])
    df <- rbind(df,filter(media_icpms_data_mlt,element =="AllEle" & is.na(chelator)))
    
    print(
      ggplot(df,aes(x=chelator_concentration,y=Sc.Normalized.CPS, colour=element))+
      geom_point()+
      geom_line(aes(group=element))+
      facet_wrap(c("metal","gas_mode"),scales="free",ncol=5)+
      scale_color_viridis_d(end=0.7)+
      theme_metallica(base_size = 21)+
      theme( panel.grid.major = element_line(size=0.5,colour="lightgray"), 
             panel.grid.minor = element_line(size=0.5,colour="lightgray"))+
      geom_vline(xintercept = 0)+
      labs(title=paste0("chelator ",chels[c]),
           x= "log 10 (Concentration of perturbed metal)", y ="Sc Normalized CPS")
    )

  }
}


##################################
### Plot all media on one plot ###
##################################

mets <- unique(media_icpms_data_mlt$metal)

pdf(paste0(plot_dir,"/media_measurements_alldata.pdf"),width = 28,height = 8)

for( m in 1:length(mets)){
  
  print(
    ggplot(filter(media_icpms_data_mlt,metal==mets[m]),aes(x = media, y=log10(CPS), colour=element))+
    geom_point(size = 3)+
    facet_wrap(c("gas_mode","metal"),scales="free",ncol=1)+
    theme_metallica(base_size = 21)+
    theme( panel.grid.major = element_line(size=0.5,colour="lightgray"), 
           panel.grid.minor = element_line(size=0.5,colour="lightgray"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_color_manual(values = colkey_Ele)+
    labs(title=paste0(mets[m]," in all samples"),x= "media", y ="log10(CPS)")
  )
}
dev.off()

############################################
### Normalize all samples to AllEle data ###
############################################

ae_dat <- filter( media_icpms_data_mlt, media=="AllEle_1x_0uM" & metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn","P","S"))%>%
          group_by(metal,gas_mode)%>%
          summarise(median_ScNorm_AE_CPS = median(Sc.Normalized.CPS,na.rm=T))%>%
          ungroup()

ae_norm_dat <- merge(media_icpms_data_mlt, ae_dat, by = c("metal","gas_mode"))%>%
               mutate(AE_norm_CPS = Sc.Normalized.CPS/median_ScNorm_AE_CPS)

### Plot AE normalized data

###################################
### Create output for conc axis ###
###################################

conc_axis_df <- filter(ae_norm_dat,metal == element & is.na(chelator))%>%
                group_by(metal,media,element_concentration)%>%
                summarise(Actual_AEnorm_Conc = round(mean(AE_norm_CPS,na.rm=T),3))%>%
                ungroup()

patterns <- c("x/100", "x/50", "x/20", "x/10", "x/5", "x/2", "0x","2x", "5x", "10x", "20x", "50x", "100x")
replacements <- c("0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "0","2", "5", "10", "20", "50", "100")

# fix all patterns
for (i in seq_along(patterns)) {
  conc_axis_df$media <- gsub(pattern = patterns[i], replacement = replacements[i], x = conc_axis_df$media)
}

colnames(conc_axis_df) <- c("metal","BioSpecID","element_concentration","actual_AEnorm_conc")

write.csv(conc_axis_df,paste0(outputtables_dir,"/measured_element_concentration_axis.csv"),row.names=F)

conc_axis_df <- conc_axis_df%>%
  group_by(metal)%>%
  mutate(Lowest_Conc = min(element_concentration,na.rm=T))%>%
  ungroup()%>%
  mutate(element = metal,
         Lowest_Conc = ifelse(element_concentration == Lowest_Conc,
                              paste(element,round(actual_AEnorm_conc,3)),NA))%>%
  group_by(element)%>%
  mutate(Highest_Conc =  max(element_concentration))%>%
  ungroup()%>%
  mutate(Highest_Conc = ifelse(Highest_Conc == element_concentration,
                               paste(element,round(actual_AEnorm_conc,2)),NA))

pdf(paste0(plot_dir,"/AE_normalized_media_concentrations.pdf"),width = 8,height=6)
ggplot(filter(conc_axis_df,element != "AllEle"),
       aes(x = log2(element_concentration),
           y= log2(actual_AEnorm_conc),
           colour = element))+
  geom_smooth(se=F,show.legend = F)+
  geom_point(size=2.5,alpha=0.8)+
  xlim(-8,7)+
  geom_label_repel(size=3,xlim=c(-8.5,-4.5),
                   aes( x= log2(element_concentration),
                        label = Lowest_Conc),show.legend = F)+
  geom_label_repel(size=3,
                   aes(  label = Highest_Conc),show.legend = F)+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x="log2(Theoretical concentration)",
       y="log2(Measured element concentration)")+
  theme( panel.grid.major = element_line(colour="lightgray"), 
         panel.grid.minor = element_line(colour="lightgray"),
         legend.title = element_blank())
dev.off()





