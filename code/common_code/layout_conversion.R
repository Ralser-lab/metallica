library(tidyverse)
library(dplyr)

# Convert Letters to number
LETTER2num <- function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}


convert_384row_to_96row <- function(x){
  
  Row_96 <- ifelse(LETTER2num(x) < 3,"A",
                   ifelse(LETTER2num(x) < 5, "B",
                          ifelse(LETTER2num(x) < 7 ,"C",
                                 ifelse(LETTER2num(x) < 9, "D",
                                        ifelse(LETTER2num(x) < 11,"E",
                                               ifelse(LETTER2num(x) < 13,"F",
                                                      ifelse(LETTER2num(x) < 15,"G",
                                                             "H")
                                               )
                                        )
                                 )
                          )
                   )
  )
  
  return(Row_96)
}


convert_384col_to_96col <- function(x){
  
  Col_96 <- ifelse(x<3,1,
                   ifelse(x<5,2,
                          ifelse(x<7,3,
                                 ifelse(x<9,4,
                                        ifelse(x<11,5,
                                               ifelse(x<13,6,
                                                      ifelse(x<15,7,
                                                             ifelse(x<17, 8, 
                                                                    ifelse(x<19,9,
                                                                           ifelse(x<21,10,
                                                                                  ifelse(x<23,11,
                                                                                         12)))))))))))
  return(Col_96)
  
}

convert_384pos_to_96post<-function(x){
  
  rw_384 <- unlist(strsplit(x,""))[[1]]
  cl_384 <- as.numeric(unlist(strsplit(x,"[A-Z]"))[[2]])
  
  pos_96 <-paste0(convert_384row_to_96row(rw_384),convert_384col_to_96col(cl_384)) 
  
  return(pos_96)
}
