require(dplyr)
require(tidyr)
require(ggplot2)
require(readxl)
require(reshape2)
require(RColorBrewer)
require(ggrepel)
library(grid)
library(robustbase)
library(corrplot)

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


# Function to compute correlations and p-values for metal- metal correlations
compute_correlation_and_pvalue <- function(df) {
  pearson_test <- cor.test(df$rel_env_element_concentration_actual, df$Ratio_to_AEngperwell, method = "pearson")
  spearman_test <- cor.test(df$rel_env_element_concentration_actual, df$Ratio_to_AEngperwell, method = "spearman")
  
  return(data.frame(pearson = pearson_test$estimate, 
                    p_value_pearson = pearson_test$p.value,
                    spearman = spearman_test$estimate,
                    p_value_spearman = spearman_test$p.value))
}
