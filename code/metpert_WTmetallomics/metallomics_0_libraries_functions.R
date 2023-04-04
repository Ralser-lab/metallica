require(dplyr)
require(tidyr)
require(ggplot2)
require(readxl)
require(robustbase)
require(reshape2)
require(RColorBrewer)
require(ggrepel)

Conc2NumConc=function(a){
  
  if(grepl("[0-9]x",a)){
    y=gsub("x","*1",a)
  }else if(grepl("x[/]",a)){
    y=gsub("x","1",a)
  }else{
    y=NA
  }
  y=eval(parse(text=y))
  return(y)
}


ChelCon2Num=function(a){
  
  if(grepl("mM",a)){
    y=1000*as.numeric(gsub("mM","",a))
  }else if(grepl("uM",a)){
    y=as.numeric(gsub("uM","",a))
  }else{
    y=NA
  }
}
