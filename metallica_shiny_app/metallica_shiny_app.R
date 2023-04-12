######################################################################
### Shiny app for visualising all data in the metallica manuscript ###
######################################################################


######################
### Load libraries ###
######################

#install.packages("shiny")
#install.packages("ggplot2")
#install.packages("webshot")
#webshot::install_phantomjs()
library(shiny)
library(plotly)
library(ggplot2)
library(viridis)


##########################
### Graphic Parameters ###
##########################


theme_metallica <- function(base_size = 25) {
  # Starts with theme_grey and then modify some parts
  theme_minimal(base_size = base_size) %+replace%
    theme(
      strip.background = element_rect(colour="black",fill=NA,size = 0.15),
      strip.text.x = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.text.y = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.switch.pad.grid = unit(0.2,"cm"),
      strip.placement = "outside",
      
      axis.text.x = element_text(size=16,vjust = 0.4),
      axis.text.y = element_text(size=16,hjust = 0.5),
      axis.ticks =  element_line(colour = "black", size = 0.2), 
      axis.title.x= element_text(size=18, vjust = -2),
      axis.title.y= element_text(size=18,angle=90,vjust=2),
      
      panel.background = element_blank(), 
      panel.border = element_rect(colour="black",fill = NA,size=0.15), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      
      plot.background = element_blank(), 
      plot.margin = unit(c(0.75,  0.75, 0.75, 0.75), "lines"),
      
      axis.line.x = element_line(color="black", size = 0.1),
      axis.line.y = element_line(color="black", size = 0.1)
    )
}


colkey_Eleconc<-c("0" = viridis(15,begin = 0.1)[[1]],
                  "2e-04" = viridis(15,begin=0.1)[[2]],
                  "0.01" = viridis(15,begin=0.1)[[3]],
                  "0.02" = viridis(15,begin=0.1)[[4]],
                  "0.05" = viridis(15,begin=0.1)[[5]],
                  "0.1" = viridis(15,begin=0.1)[[6]],
                  "0.2" = viridis(15,begin=0.1)[[7]],
                  "0.5" = viridis(15,begin=0.1)[[8]],
                  "1" = viridis(15,begin=0.1)[[9]],
                  "2" = viridis(15,begin=0.1)[[10]],
                  "5" = viridis(15,begin=0.1)[[11]],
                  "10" = viridis(15,begin=0.1)[[12]],
                  "20" = viridis(15,begin=0.1)[[13]],
                  "50" = viridis(15,begin=0.1)[[14]],
                  "100" = viridis(15,begin=0.1)[[15]])

colkey_Ele <-c("AllEle Control" = "#EEE0E5",
               "Ca" = RColorBrewer::brewer.pal(12,"Paired")[[1]],
               "Cu" = RColorBrewer::brewer.pal(12,"Paired")[[4]],
               "Fe" = RColorBrewer::brewer.pal(12,"Paired")[[6]],
               "K" = "#838B83",
               "Mg" = RColorBrewer::brewer.pal(12,"Paired")[[8]],
               "Mn" = RColorBrewer::brewer.pal(12,"Paired")[[10]],
               "Mo" = "#00868B",
               "Na" = "#DAA520",
               "Zn" = "#000080")

colkey_EleDir <-c("AllEle Control" = "#EEE0E5",
                  "Ca Depletion" = RColorBrewer::brewer.pal(12,"Paired")[[1]],
                  "Ca Excess" = RColorBrewer::brewer.pal(12,"Paired")[[2]],
                  "Cu Depletion" = RColorBrewer::brewer.pal(12,"Paired")[[3]] ,
                  "Cu Excess" = RColorBrewer::brewer.pal(12,"Paired")[[4]],
                  "Fe Depletion" = RColorBrewer::brewer.pal(12,"Paired")[[5]],
                  "Fe Excess" = RColorBrewer::brewer.pal(12,"Paired")[[6]],
                  "K Depletion" = "#838B83",
                  "K Excess" = "#C1CDC1",
                  "Mg Depletion" = RColorBrewer::brewer.pal(12,"Paired")[[7]],
                  "Mg Excess" = RColorBrewer::brewer.pal(12,"Paired")[[8]],
                  "Mn Depletion" = RColorBrewer::brewer.pal(12,"Paired")[[9]],
                  "Mn Excess" = RColorBrewer::brewer.pal(12,"Paired")[[10]],
                  "Mo Depletion" = "#00C5CD",
                  "Mo Excess" = "#00868B",
                  "Na Depletion" = "#FFD700",
                  "Na Excess" = "#DAA520",
                  "Zn Depletion" = "#4876FF",
                  "Zn Excess" = "#000080")


colkey_BioSpecID <-c("AllEle 1" = "#EEE0E5",
                     
                     # Boron
                     "B 0"  =  grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[1]],
                     "B 0.01"  = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[2]],
                     "B 0.02" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[3]],    
                     "B 0.05" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[4]],
                     "B 0.1" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[5]],
                     "B 0.2" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[6]],
                     "B 0.5" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[7]],
                     "B 2" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[8]],
                     "B 5" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[9]],
                     "B 10" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[10]],
                     "B 20" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[11]],
                     "B 50" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[12]],
                     "B 100" = grDevices::colorRampPalette(c("#878787","#0D0D0D"))(13)[[13]],
                     
                     # Calcium
                     "Ca 0"  =  grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[1]],
                     "Ca 0.01"  = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[2]],
                     "Ca 0.02" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[3]],    
                     "Ca 0.05" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[4]],
                     "Ca 0.1" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[5]],
                     "Ca 0.2" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[6]],
                     "Ca 0.5" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[7]],
                     "Ca 2" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[8]],
                     "Ca 5" = grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[9]],
                     
                     # Copper 
                     
                     "Cu 0"  =  grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     "Cu 0.01"  = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[2]],
                     "Cu 0.02" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[3]],    
                     "Cu 0.05" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[4]],
                     "Cu 0.1" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[5]],
                     "Cu 0.2" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[6]],
                     "Cu 0.5" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[7]],
                     "Cu 2" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[8]],
                     "Cu 5" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[9]],
                     "Cu 10" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[10]],
                     "Cu 20" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[11]],
                     "Cu 50" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[12]],
                     "Cu 100" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[13]],
                     
                     # Copper Chelators
                     
                     "Cu 0 BCS 100" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     "Cu 0 BCS 200" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     "Cu 0 BCS 300" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     "Cu 0 BCS 400" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     "Cu 0 BCS 500" = grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                     
                     # Iron 
                     
                     "Fe 0"  =  grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                     "Fe 0.01"  = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[2]],
                     "Fe 0.02" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[3]],    
                     "Fe 0.05" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[4]],
                     "Fe 0.1" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[5]],
                     "Fe 0.2" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[6]],
                     "Fe 0.5" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[7]],
                     "Fe 2" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[8]],
                     "Fe 5" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[9]],
                     "Fe 10" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[10]],
                     "Fe 20" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[11]],
                     "Fe 50" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[12]],
                     "Fe 100" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[13]],
                     
                     # Iron Chelators
                     
                     "Fe 0 DiP 50" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                     "Fe 0 DiP 100" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                     "Fe 0 DiP 150" = grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                     
                     # Potassium 
                     
                     "K 0.1"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[1]],
                     "K 0.2"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[2]],
                     "K 0.5"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[3]],
                     "K 2"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[4]],
                     "K 5"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[5]],
                     "K 10"  =  grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[6]],
                     
                     # Magnesium 
                     "Mg 0.05"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[1]],
                     "Mg 0.1"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[2]],
                     "Mg 0.2"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[3]],
                     "Mg 0.5"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[4]],
                     "Mg 2"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[5]],
                     "Mg 5"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[6]],
                     "Mg 10"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[7]],
                     "Mg 20"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[8]],
                     "Mg 50"  =  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[9]],
                     
                     # Manganese 
                     
                     "Mn 0"  =  grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[1]],
                     "Mn 0.01"  = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[2]],
                     "Mn 0.02" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[3]],    
                     "Mn 0.05" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[4]],
                     "Mn 0.1" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[5]],
                     "Mn 0.2" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[6]],
                     "Mn 0.5" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[7]],
                     "Mn 2" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[8]],
                     "Mn 5" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[9]],
                     "Mn 10" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[10]],
                     "Mn 20" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[11]],
                     "Mn 50" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[12]],
                     "Mn 100" = grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[13]],
                     
                     # Molybdenum 
                     
                     "Mo 0"  =  grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[1]],
                     "Mo 0.01"  = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[2]],
                     "Mo 0.02" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[3]],    
                     "Mo 0.05" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[4]],
                     "Mo 0.1" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[5]],
                     "Mo 0.2" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[6]],
                     "Mo 0.5" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[7]],
                     "Mo 2" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[8]],
                     "Mo 5" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[9]],
                     "Mo 10" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[10]],
                     "Mo 20" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[11]],
                     "Mo 50" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[12]],
                     "Mo 100" = grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[13]],
                     
                     # Sodium 
                     
                     "Na 0"  =  grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[1]],
                     "Na 0.01"  = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[2]],
                     "Na 0.02" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[3]],    
                     "Na 0.05" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[4]],
                     "Na 0.1" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[5]],
                     "Na 0.2" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[6]],
                     "Na 0.5" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[7]],
                     "Na 2" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[8]],
                     "Na 5" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[9]],
                     "Na 10" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[10]],
                     "Na 20" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[11]],
                     "Na 50" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[12]],
                     "Na 100" = grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[13]],
                     
                     # Zinc
                     "Zn 0" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[1]],
                     "Zn 2e-04" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[2]],
                     "Zn 0.01" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[3]],    
                     "Zn 0.02" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[4]],
                     "Zn 0.05" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[5]],
                     "Zn 0.1" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[6]],
                     "Zn 0.2" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[7]],
                     "Zn 0.5" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[8]],
                     "Zn 2" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[9]],
                     "Zn 5" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[10]],
                     "Zn 10" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[11]],
                     "Zn 20" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[12]],
                     "Zn 50" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[13]],
                     "Zn 100" = grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[14]]
                     
                     
)

### Colour Key for Datasets ###

dataset_colkey <- c("Meas_in_any" = "#EEDC82",
                    "Hit_in_any"  = 	"#B0171F",
                    "SGD_anno_unchar" = "#CDCDC1",
                    "In_GO_MF" = "#9BCD9B",
                    "In_MetPDBessmet" = "#DDA0DD" ,
                    "In_KOgrowth" = "#CD8500",
                    "In_hdPCA" = "#A2CD5A",
                    "In_proteomics"= "#104E8B")

dataset_colkey_darker <- c("Meas_in_any" = "#EEC900",
                           "Hit_in_any"  = 	"#B0171F",
                           "SGD_anno_unchar" = "#CDCDC1",
                           "In_GO_MF" = "#9BCD9B",
                           "In_MetPDBessmet" = "#DDA0DD" ,
                           "In_KOgrowth" = "#CD8500",
                           "In_hdPCA" = "#A2CD5A",
                           "In_proteomics"= "#104E8B")

###################3
### input datasets 
#################

app_data_dir <-"/Users/aulakhs/Documents/Ralser Lab/metallica/metallica_shiny_app/data"

# metpert WT metallomics dataset
metpertWTmetallomics_data <-  read.csv(paste0(app_data_dir,"/metallica_app_metpertWTmetallomics.csv"),stringsAsFactors = F)%>%
                     mutate(`metal perturbed` = element_perturbed,
                            `relative environmental concentration` = rel_element_concentration_actual,
                            `metal measured`= element_measured,
                             `relative intracellular concentration` = Ratio_to_AEngperwell)%>%
                     dplyr::select(BioSpecID,`metal perturbed`,`relative environmental concentration`,`metal measured`,`relative intracellular concentration`)


# metpert WT growth dataset -- growth rates only

metpertWTgrowth_data <-  read.csv(paste0(app_data_dir,"/metallica_app_metpertWTgrowthrate.csv"),stringsAsFactors = F)
  


# y5k metrel met specific data 
metpertWTproteomics_data <- read.csv(paste0(app_data_dir,"/metallica_app_y5kmetspecificKOs.csv"),stringsAsFactors = F)

#y5k_metrelmetspecific_data <-  read.csv(paste0(app_data_dir,"/metallica_app_y5kmetrelmetspecificKOs.csv"),stringsAsFactors = F)

# y5k KO metallomics 

KOmetallomics_data <-  read.csv(paste0(app_data_dir,"/metallica_app_KOmetallomics.csv"),stringsAsFactors = F)


# OE data

OEmetallomics_data <-  read.csv(paste0(app_data_dir,"/metallica_app_OEmetallomics.csv"),stringsAsFactors = F)

######################
### User Interface ###
######################

ui <- fluidPage(
  tags$head(
    tags$link(href = "https://fonts.googleapis.com/css2?family=Roboto&display=swap", rel = "stylesheet"),
  tags$style(HTML("
            body {
              background-color: #152238;
              font-family: 'Roboto', sans-serif;
              color: white;
            }
            .titlePanel {
              color: white;
            }
            .well {
              background-color: white;
              color: black;
              border-radius: 15px;
            }
            .nav-tabs > li > a {
              color: white;
            }
            .plot-container {
              border-radius: 10px;
              margin: 15px 0px; 
            }
            .description-container {
              border-radius: 10px;
              overflow: visible;
               margin-top: 20px;
            }
            .video-container {
              border-radius: 10px;
              overflow: hidden;
              margin-top: 20px;
              position: relative;
              padding-bottom: 56.25%;
              padding-top: 25px;
              height: 0;
            }
            .video-container iframe {
              position: absolute;
              top: 0;
              left: 0;
              width: 100%;
              height: 100%;
            }
            
      "))
  ),
  titlePanel("The Proteomic Landscape of Cellular Metal Ion Homeostasis", windowTitle = NULL),
  mainPanel(
    width = 12,
    tabsetPanel(
      tabPanel("Overview",
               fluidRow(
                 column(12, 
                        div(class = "description-container",
                        wellPanel(
                          p(HTML("<i>Simran Kaur Aulakh <sup>1,2</sup></i>,
                                  <i>Stephan Kamrad <sup>1</sup></i>, 
                                  <i>Oliver Lemke <sup>3</sup></i>,
                                  <i>Lukasz Szyrwiel <sup>3</sup></i>,
                                  <i>Johannes Hartl <sup>3</sup></i>,
                                  <i>Michael Muelleder <sup>3</sup></i> 
                                  <i>and Markus Ralser <sup>1,2,3,#</sup></i>")),
                          p(HTML("<sup>1</sup> Molecular Biology of Metabolism Laboratory,
                                               The Francis Crick Institute, 
                                               1 Midland Road, 
                                               NW1 1AT, 
                                               London,
                                               United Kingdom <br>
                          <sup>2</sup> Wellcome Centre for Human Genetics,
                                       University of Oxford
                                       Roosevelt Drive, 
                                       OX3 7BN
                                       Oxford,
                                       United Kingdom <br>
                          <sup>3</sup> Department of Biochemistry,
                                               Charité – Universitätsmedizin Berlin,
                                               Charitéplatz 1,
                                               10117, 
                                                Berlin, 
                                                Germany <br>
                          <sup>#</sup> All correspondence should be addressed to markus.ralser@charite.de ")),
                     
                          tags$b("Abstract"), 
                          p(HTML("Metal ions have played a crucial role in the origin and evolution of
                            metabolism and are required for 80% of metabolic pathways.
                            Extracellular metal ion concentrations can vary greatly and 
                            are buffered by cells.  However, the underlying biochemical 
                            network that controls intracellular metal homeostasis is still 
                            poorly understood. In this study, we manipulated the concentration 
                            of each bioessential metal ion in <i>S.cerevisiae</i> and generated a 
                            systematic, quantitative interaction network monitoring growth, 
                            ionome and proteome. Our resource quantifies the cellular ability to
                            buffer ion concentrations and identifies new protein-metal interactions,
                            including 67 proteins whose function remained elusive.
                            Our study provides a comprehensive description of the proteome-wide 
                            landscape of the
                            cellular response to altered metal ion concentrations.")),
                          tags$ul(
                            tags$li(tags$a(href = "https://doi.org/your_doi_here", target = "_blank", "DOI Link")),
                            tags$li(tags$a(href = "https://github.com/your_username/your_repository", target = "_blank", "GitHub Repository"))
                          )
                        ),
                        h4("Web App Tutorial"),
                        div(class = "video-container",
                            tags$iframe(src = paste0("https://www.youtube.com/embed/T75IKSXVXlc"),
                                        frameborder = "0",
                                        allow = "autoplay; fullscreen",
                                        allowfullscreen = "allowfullscreen")
                        ),
                        column(12, style = "height: 20px;"), # add blank column with increased height for spacing
                        div(class = "select-container", style = "display: inline-block;",
                                selectInput("selected_main_metal", 
                                            label = "Choose a metal:",
                                            choices = c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"),
                                            selected = "Ca",
                                            width = "150px")
                            
                        )
                        )
                 )
               ),
               fluidRow(
                 column(6,
                        h4(HTML("Intracellular metallomics -  WT <i>S.cerevisiae</i>")),
                        div(class = "plot-container", plotOutput("pointPlot1"))),
                 column(6, 
                        h4(HTML("Growth rate - WT <i>S.cerevisiae</i>")),
                        div(class = "plot-container", plotOutput("pointPlot2")))
               ),
               fluidRow(
                 column(3,
                        h4(HTML("Proteomcis -  WT <i>S.cerevisiae</i> in metal perturbation media")),
                        div(class = "plot-container", plotOutput("barPlot1"))),
                 column(3,
                        h4(HTML("Proteomcis - deletions of metal related <i>S.cerevisiae</i> genes")),
                        div(class = "plot-container", plotOutput("barPlot2"))),
                 column(3,
                        h4(HTML("Intracellular metallomics - <i>S.cerevisiae</i> KO library")),
                        div(class = "plot-container", plotOutput("barPlot3"))),
                 column(3,
                        h4(HTML("Intracellular metallomics - <i>S.cerevisiae</i> OE library")),
                        div(class = "plot-container", plotOutput("barPlot4")))
               )),
      
      ## metal perturbation proteomics tab 

      tabPanel("WT - Metal Perturbation - Proteomics",
               fluidRow(
                 column(12,
                        div(class = "description-container",
                            wellPanel(
                              p("This is a description of the MetPert Proteomics dataset. For more information, please visit the following resources:"),
                              tags$ul(
                                tags$li(tags$a(href = "https://doi.org/your_doi_here", target = "_blank", "DOI Link")),
                                tags$li(tags$a(href = "https://github.com/your_username/your_repository", target = "_blank", "GitHub Repository"))
                              )
                            ),
                            selectInput("selected_metpert_protein", 
                                        label = "Choose a protein:",
                                        choices = sprintf("Protein %04d", 1:1800),
                                        selected = "Protein 0001")
                        )
                 )
               ),
               fluidRow(
                 lapply(c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"), function(metal) {
                   column(4, div(class = "plot-container", plotOutput(paste0("proteomicsPlot", metal))))
                 })
               )
      ),
      
      #########################################################################
      ### y5k proteomics data - filtered KOs with metal specific annotation ###
      #########################################################################
      
      tabPanel("y5k KO Proteomics",
               fluidRow(
                 column(12,
                        div(class = "description-container",
                            wellPanel(
                              # Add other elements or text here if needed
                            ),
                            div(style = "display: inline-block;text-align: center;",
                                selectInput("selected_metal_proteomics", 
                                            label = "Choose a metal:",
                                            choices = c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"),
                                            selected = "Ca")
                            ),
                            div(style = "display: inline-block;text-align: center; margin-left: 40px",
                                selectizeInput("selected_KO_gene_proteomics", 
                                               label = "Choose a metal associated knock-out:",
                                               choices = NULL,
                                               selected = NULL)
                            )
                        )
                 )
               ),
               fluidRow(
                 column(6,
                        plotOutput("volcano_plot")
                 ),
                 column(6,
                        div(style = "background-color: white; padding: 15px;",
                            fluidRow(
                              column(6,
                                     tableOutput("top10_table_part1")
                              ),
                              column(6,
                                     tableOutput("top10_table_part2")
                              )
                            )
                        )
                 )
               )
      ),
      
      tabPanel("y5k KO Growth in Metal Depletion",),
      tabPanel("hdPCA in Metal Depletion",),
      tabPanel("y5k KO + OE -  Metallomics",
               fluidRow(
                 column(6,
                        div(style = "text-align: center;",
                            div(class = "description-container",
                                tags$div(style = "display: inline-block;",
                                         selectInput("selected_KOmetallomics_metal", 
                                                     label = "Choose metal- KO mutant data:",
                                                     choices = c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"),
                                                     selected = "Mn")
                                ),
                                tags$div(style = "display: inline-block; margin-left: 20px;",
                                         selectizeInput("selected_KOmetallomics_KOgene", 
                                                        label = "Choose KO Gene:",
                                                        choices = NULL, # Initially empty, will be populated in the server function
                                                        selected = NULL)
                                )
                            )
                        )
                 ),
                 column(6,
                        div(style = "text-align: center;",
                            div(class = "description-container",
                                tags$div(style = "display: inline-block;",
                                         selectInput("selected_OEmetallomics_metal", 
                                                     label = "Choose metal - OE mutant data:",
                                                     choices = c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"),
                                                     selected = "Mn")
                                ),
                                tags$div(style = "display: inline-block; margin-left: 20px;",
                                         selectizeInput("selected_OEmetallomics_OEgene", 
                                                        label = "Choose OE Gene:",
                                                        choices = NULL, # Initially empty, will be populated in the server function
                                                        selected = NULL)
                                )
                            )
                        )
                 )
               )
               ,
               fluidRow(
                 column(6,
                        h4("metallomics Z-scores of KO mutants"),
                        plotOutput("KOzscore_density_plot"),
                        tags$div(
                          style = "text-align: center; margin-top: 20px;",
                          downloadButton("download_ko_png", "Download KO Plot as PNG"),
                          tags$span(" ", style = "padding: 0 65px;"),
                          downloadButton("download_ko_pdf", "Download KO Plot as PDF")
                        )
                 ),
                 column(6,
                        h4("metallomics Z-scores of OE mutants"),
                        plotOutput("OEzscore_density_plot"),
                        tags$div(
                          style = "text-align: center; margin-top: 20px;",
                          downloadButton("download_oe_png", "Download OE Plot as PNG"),
                          tags$span(" ", style = "padding: 0 65px;"),
                          downloadButton("download_oe_pdf", "Download OE Plot as PDF")
                        )
                 )
               )
      )
      
      ,
      fluidRow(
        column(12, style = "height: 18px;"), # add blank column with increased height for spacing
        column(12,
               div(style = "background-color: white; padding: 20px; border-radius: 10px; color: black; text-align: center;",
                   div(style = "display: flex; justify-content: center; align-items: center;",
                       img(src = "cricklogo.png", alt = "Francis Crick Institute", width = "12%", height = "12%",style = "margin-right: 50px;"),
                       img(src = "charitelogo.png", alt = "Charité-Universitätsmedizin Berlin", width = "18%", height = "18%",style = "margin-right: 40px;"),
                       img(src = "oxfordlogo.png", alt = "University of Oxford", width = "25%", height = "25%")
                   ),
                   br(),
                   div("This webpage was created by Simran Kaur Aulakh using the",
                       tags$a(href = "https://shiny.rstudio.com", target = "_blank", "Shiny"), "R library")
               )
        )
      )
    )
  )
)


# Define the server
server <- function(input, output, session) {
  
  #####################
  ### summary plots ###
  #####################
  
  # metallomics env vs intracellular
  
  metpertWTmetallomics_data_filtered_pointPlot1 <- reactive({
                                   metpertWTmetallomics_data%>%
                                   filter(`metal measured` == input$selected_main_metal,
                                          `metal perturbed` == input$selected_main_metal)
  })
  
  output$pointPlot1 <- renderPlot({
    ggplot(metpertWTmetallomics_data_filtered_pointPlot1(), aes(x = log2(`relative environmental concentration`), 
                                                       y = log2(`relative intracellular concentration`),
                                                       color = BioSpecID)) +
      geom_point(size = 4,alpha = 0.8) +
      scale_color_manual(values = colkey_BioSpecID) +
      theme_metallica()+
      theme(legend.position = "none")+
      labs(x = expression(log[2]~relative~environmental~metal~concentration),
           y = expression(log[2]~relative~intracellular~metal~concentration))
    
  })
  
  metpertWTgrowth_data_filtered_pointPlot2 <- reactive({
    metpertWTgrowth_data%>%
      filter(`metal` == input$selected_main_metal)
  })
  
  output$pointPlot2 <- renderPlot({
    ggplot(metpertWTgrowth_data_filtered_pointPlot2(), aes(x = log2(`env_metal_concentration`), 
                                                                y = growth_rate,
                                                                color = BioSpecID)) +
      geom_point(size = 4,alpha = 0.8) +
      scale_color_manual(values = colkey_BioSpecID) +
      theme_metallica()+
      theme(legend.position = "none")+
      ylim(range(metpertWTgrowth_data$growth_rate))+
      labs(x = expression(log[2]~relative~environmental~metal~concentration),
           y = expression(growth~reate))
    
  })
  
  output$barPlot1 <- renderPlot({
    ggplot() +
      geom_bar(aes(x = factor(1), y = rnorm(1, input$selected_main_metal == "Fe")), stat = "identity", width = 0.5) +
      labs(title = "", x = "", y = "number of hits") +
      theme_metallica()
  })
  
  output$barPlot2 <- renderPlot({
    ggplot() +
      geom_bar(aes(x = factor(1), y = rnorm(1, input$selected_main_metal == "Fe")), stat = "identity", width = 0.5) +
      labs(title = "", x = "", y = "number of hits") +
      theme_metallica()
  })
  
  
  output$barPlot3 <- renderPlot({
    ggplot() +
      geom_bar(aes(x = factor(1), y = rnorm(1, input$selected_main_metal == "Fe")), stat = "identity", width = 0.5) +
      labs(title = "", x = "", y = "number of hits") +
      theme_metallica()
  })
  
  output$barPlot4 <- renderPlot({
    ggplot() +
      geom_bar(aes(x = factor(1), y = rnorm(1, input$selected_main_metal == "Fe")), stat = "identity", width = 0.5) +
      labs(title = "", x = "", y = "number of hits") +
      theme_metallica()
  })
  
  
  ###################################################
  ### Data for metal perturbation  proteomics tab ###
  ###################################################
  
  proteomicsData <- reactive({
    # Generate random data for the example; replace with your actual data
    lapply(c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"), function(metal) {
      list(x = 1:10, y = rnorm(10))
    })
  })
  
  # Create proteomics plots
  for (i in seq_along(c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn"))) {
    metal <- c("Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Zn")[i]
    output[[paste0("proteomicsPlot", metal)]] <- renderPlot({
      ggplot() +
        geom_point(aes(x = proteomicsData()[[i]]$x, y = proteomicsData()[[i]]$y), color = colkey_Ele[metal]) +
        geom_smooth(aes(x = proteomicsData()[[i]]$x, y = proteomicsData()[[i]]$y), method = "loess", se = FALSE, color = "black") +
        labs(title = paste("Proteomics", metal, "Plot"),
             x = "X-axis",
             y = "Y-axis") +
        theme_metallica()
    })
  }
  
  
  #############################################################################################
  ### Server code y5k proteomics data filtered for KOs with known metal specific annotation ###
  #############################################################################################
  
  # Update the choices for the KO gene input based on the selected metal
  observe({
    filtered_data <- subset(metpertWTproteomics_data, metal_connec2KO == input$selected_metal_proteomics)
    unique_KO_genes <- unique(filtered_data$KO_gene)
    updateSelectizeInput(session, "selected_KO_gene_proteomics", choices = unique_KO_genes, selected = unique_KO_genes[1])
  })
  
  # Filter the data based on the selected metal and KO gene
  filtered_proteomics_data <- reactive({
    subset(metpertWTproteomics_data, metal_connec2KO == input$selected_metal_proteomics & KO_gene == input$selected_KO_gene_proteomics)
  })
  
  
  output$volcano_plot <- renderPlot({
    filtered_data <- filtered_proteomics_data()
    threshold <- subset(filtered_data, p_value < 0.05 & abs(log2FC) > log2(1.5))
    
    ggplot(filtered_data, aes(x = log2FC, y = -log10(p_value))) +
      geom_point(aes(color = (p_value < 0.05) & (abs(log2FC) > log2(1.5)), alpha = 0.8)) +
      geom_text(data = threshold, aes(label = measured_protein), hjust = 1.5, vjust = 0.5, size = 3) +
      scale_color_manual(values = c("grey", colkey_Ele[input$selected_metal_proteomics])) +
      labs(x = expression(paste("log"[2], "(fold difference vs control)")),
           y = expression(paste("-log"[10], "(p-value)"))) +
      theme_metallica() +
      theme(legend.position = "none",
            plot.title = element_blank())
  })
  
  
  output$top10_table <- renderTable({
    y5k_prot_filtered_data <- filtered_proteomics_data()
    significant_data <- subset(y5k_prot_filtered_data, p_value < 0.05)
    sorted_data <- significant_data[order(abs(significant_data$log2FC), decreasing = TRUE),]
    top10 <- head(sorted_data, 10)
    top10[, c("measured_protein", "log2FC", "p_value")]
  })
  
  
  ###################################
  ###################################
  ### Metallomics data to display ###
  ###################################
  ###################################
  
  ######################
  ### KO metallomics ###
  ######################
  
  # filter KO metallomics df based on metal selected
  KOmetallomics_filtered_data <- reactive({
    subset(KOmetallomics_data, metal == input$selected_KOmetallomics_metal)
  })
  
  ## list of genes available 
  observe({
    unique_genes <- unique(KOmetallomics_data$KOgenename)
    updateSelectizeInput(session, "selected_KOmetallomics_KOgene", choices = unique_genes, selected = "SOD2")
  })
  
  # note z score of selected gene
  selected_gene_KOzscore <- reactive({
    KOmetallomics_filtered_data()$zscore_KOmetallomics[KOmetallomics_filtered_data()$KOgenename == input$selected_KOmetallomics_KOgene]
  })
  
  create_ko_density_plot <- function() {
    filtered_data <- KOmetallomics_filtered_data()
    p <- ggplot(filtered_data,
                aes(x = zscore_KOmetallomics)) +
      geom_density(aes(fill = metal,
                       colour = metal),
                   alpha = 0.3) +
      geom_vline(xintercept = selected_gene_KOzscore(),
                 color = "black",
                 linetype = "dashed",
                 size = 0.7) +
      geom_label( x = selected_gene_KOzscore(),
                  y = 0.4,
                  size = 5,
                  label = paste(input$selected_KOmetallomics_KOgene, "\n", round(selected_gene_KOzscore(), 3)),
                  color = "black") +
      labs(title = " ",
           x = "Zscores of all KOs",
           y = "density") +
      xlim(-10, 10) +
      scale_fill_manual(values = colkey_Ele) +
      scale_colour_manual(values = colkey_Ele) +
      theme_metallica() +
      theme(legend.position = "none")
    
    return(p)
  }
  
  # Modify the output$KOzscore_density_plot to use the create_ko_density_plot() function
  output$KOzscore_density_plot <- renderPlot({
    create_ko_density_plot()
  })
  
  # Modify the downloadHandler content to use webshot for downloading the PNG
  output$download_ko_png <- downloadHandler(
    filename = function() {
      paste("KO_density_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Save the plot as a temporary HTML file
      temp_html_file <- tempfile(fileext = ".html")
      htmlwidgets::saveWidget(create_ko_density_plot(), temp_html_file, selfcontained = FALSE)
      
      # Use webshot to convert the HTML file to a PNG file
      webshot::webshot(url = paste0("file:///", temp_html_file),
                       file = file,
                       selector = ".plot-container",
                       delay = 0.5)
    }
  )
  
  # Modify the downloadHandler content to use the create_ko_density_plot() function
  output$download_ko_pdf <- downloadHandler(
    filename = function() {
      paste("KO_density_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(create_ko_density_plot())
      dev.off()
    }
  )
  
  ######################
  ### OE metallomics ###
  ######################
  
  # filter KO metallomics df based on metal selected
  OEmetallomics_filtered_data <- reactive({
    subset(OEmetallomics_data, metal == input$selected_OEmetallomics_metal)
  })
  
  ## list of genes available 
  observe({
    unique_genes <- unique(OEmetallomics_data$OEgenename)
    updateSelectizeInput(session, "selected_OEmetallomics_OEgene", choices = unique_genes, selected = "SOD2")
  })
  
  # note z score of selected gene
  selected_gene_OEzscore <- reactive({
    OEmetallomics_filtered_data()$zscore_OEmetallomics[OEmetallomics_filtered_data()$OEgenename == input$selected_OEmetallomics_OEgene]
  })
  
  # Create a separate function to generate the ggplot object
  create_oe_density_plot <- function() {
    filtered_data <- OEmetallomics_filtered_data()
    ggplot(filtered_data,
           aes(x = zscore_OEmetallomics)) +
      geom_density(aes(fill = metal,
                       colour = metal),
                   alpha = 0.3) +
      geom_vline(xintercept = selected_gene_OEzscore(),
                 color = "black",
                 linetype = "dashed",
                 size = 0.7) +
      geom_label(x = selected_gene_OEzscore(),
                 y = 0.4,
                 size = 5,
                 label = paste(input$selected_OEmetallomics_OEgene, "\n", round(selected_gene_OEzscore(), 3)),
                 color = "black") +
      labs(title = " ",
           x = "Zscores of all OEs",
           y = "density") +
      xlim(-10, 10) +
      scale_fill_manual(values = colkey_Ele) +
      scale_colour_manual(values = colkey_Ele) +
      theme_metallica() +
      theme(legend.position = "none")
  }
  
  # Modify the renderPlot to use the create_oe_density_plot() function
  output$OEzscore_density_plot <- renderPlot({
    create_oe_density_plot()
  })
  
  ## download as png button 
  output$download_oe_png <- downloadHandler(
    filename = function() {
      paste("OE_density_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Save the plot as a temporary HTML file
      temp_html_file <- tempfile(fileext = ".html")
      htmlwidgets::saveWidget(create_oe_density_plot(), temp_html_file, selfcontained = FALSE)
      
      # Use webshot to convert the HTML file to a PNG file
      webshot::webshot(url = paste0("file:///", temp_html_file),
                       file = file,
                       selector = ".plot-container",
                       delay = 0.5)
    }
  )
  
  output$download_oe_pdf <- downloadHandler(
    filename = function() {
      paste("OE_density_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(create_oe_density_plot())
      dev.off()
    }
  )
  
  
  
}

# Run the app
shinyApp(ui = ui, server = server)