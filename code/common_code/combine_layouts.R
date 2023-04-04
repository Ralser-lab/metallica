#`---
#`  Title: "combine_layouts.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 22 April 2021
#`  Description: Script to combine all 96 well layouts used for cell growth into one table
#`---

####################
### source paths ###
####################

source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")

### Make layout File ###

lo_fns <- dir(paste0(lo_dir,"/media4cellgrowth_layouts"))

layout_allplates <- vector()

for( i in 1:length(lo_fns)){
  
  layout <- read.csv(paste0(lo_dir,"/media4cellgrowth_layouts/",lo_fns[i]), stringsAsFactors = F)
  colnames(layout) <- c("Row",1:12)
  layout_melted <- na.omit(melt(layout,id.vars = c("Row")))
  layout_melted$position <- paste0(layout_melted$Row,layout_melted$variable)
  layout_melted <- layout_melted[,-c(1,2)]
  colnames(layout_melted) <- c("condition","position")
  layout_melted$plate <- i
  
  layout_allplates <- rbind(layout_allplates,layout_melted)
  
}

layout_allplates$plate_position <- paste(layout_allplates$plate,layout_allplates$position,sep="_")
layout_allplates$condition <- gsub("AllEle_prepQC","AllEleprepQC",layout_allplates$condition)


patterns <- c("x/100", "x/50", "x/20", "x/10", "x/5", "x/2", "0x","2x", "5x", "10x", "20x", "50x", "100x")
replacements <- c("0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "0","2", "5", "10", "20", "50", "100")

# fix all patterns
for (i in seq_along(patterns)) {
  layout_allplates$condition <- gsub(pattern = patterns[i], 
                                     replacement = replacements[i], 
                                     x =  layout_allplates$condition)
}

layout_allplates$condition <- gsub("_"," ",layout_allplates$condition)

write.csv(layout_allplates,paste0(proj_dir,"/accessory_files/all_layouts_combined.csv"),row.names = F)

