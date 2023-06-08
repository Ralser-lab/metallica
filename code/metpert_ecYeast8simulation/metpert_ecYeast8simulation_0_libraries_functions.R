library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(readxl)


library(purrr)
library(stringr)


# function to map colours to ORFs based on median_perc_abs_deviation_ORF

convert_percabsdeviation_to_colour <- function(pad) {
  #  pad - is a percentage absolute deviation of an ORF
  # check input
  if(!is.numeric(pad) || any(pad < 0)) {
    stop("Input should be numeric and greater than 0")
  }
  
  # Define color palette
  pad_palette <- viridis(11, option = "C",direction = -1, begin = 0.2)
  
  # For values between 0-100, choose a color linearly from first 10 colors
  # For values > 100, choose the 11th color
  if (pad <= 100) {
    index <- floor(pad / 10) + 1
    color <- pad_palette[index]
  } else {
    color <- pad_palette[11]
  }
  
  return(color)
}
