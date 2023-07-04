library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)

library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(readxl)


plot_DIANNstat_distrb <- function(df,plotname){
  
  df <- reshape2::melt(df,id.vars = "File.Name")
  pdf(paste0("QC stats ",plotname,".pdf"),width=20,height=30)
  print(
    ggplot(df,
           aes(x=value))+
      geom_histogram(bins=75,fill ="#191970")+
      facet_wrap("variable",scales="free",ncol=3)+
      theme_metallica()+
      theme(axis.text.x = element_text(angle=90))
  )
  dev.off()
}
