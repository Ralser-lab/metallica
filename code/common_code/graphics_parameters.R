library(showtext)
library(viridis)

font_add_google("Tinos")
showtext_auto()


theme_metallica <- function(base_size = 25) {
  # Starts with theme_grey and then modify some parts
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "Tinos"),
      strip.background = element_rect(colour="black",fill=NA,size = 0.15),
      strip.text.x = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.text.y = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.switch.pad.grid = unit(0.2,"cm"),
      strip.placement = "outside",
      
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16,hjust=1),
      axis.ticks =  element_line(colour = "black", size = 0.2), 
      axis.title.x= element_text(size=18),
      axis.title.y= element_text(size=18,angle=90,vjust=3),
      
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



element2colour<- function(){
  
  coldf<-rbind( cbind("Ca,","#87CEFF"),
         cbind("Cu","#2E8B57"),
         cbind("Fe","#8B3626"),
         cbind("K","#7A378B"),
         cbind("Mg","#C76114"),
         cbind("Mn","#EE799F"),
         cbind("Mo","#BC8F8F"),
         cbind("Na","#EEC900"),
         cbind("Zn","#104E8B"),
         cbind("B,","#1C1C1C"))
  
  colnames(coldf) <- c("Element","Colour")
  return(coldf)
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


cbio <- data.frame(cbind(hexcode=colkey_BioSpecID))
cbio$BioSpecID <- rownames(cbio)

write.csv(cbio,paste0(acc_file_dir,"/EPP_Colourkey_BioSpecID.csv"),row.names = F)



###############################
### Colour Key for Datasets ###
###############################

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




      