#load the library
library(dplyr)
library(magrittr)
### functions 
##functions of read the ascn segmentation file 
read.cn <- function(file) {
  cn <- read.delim(file=file, stringsAsFactors = F, header = T)
  cn %<>% set_colnames(c("SampleName", colnames(cn)[2:length(colnames(cn))]))
  return(cn)
  
}

##functions of read the ascn model file to get the ploidy value, cellurarity value
read.model <- function(file) {
  model <- read.delim(file=file, stringsAsFactors = F, header = F)
  SampleName = model[1, 2]
  ploidy.ASCAT = model[5, 2]
  ploidy.median = model[6, 2]
  cellularity=model[9, 2]
  df = as.data.frame(SampleName)
  df$Celluarity = cellularity
  df$ploidy.ASCAT = ploidy.ASCAT
  df$ploidy.median = ploidy.median
  return(df)
  
}

####functions of read the csb segmentation file 
read.cbs <- function(file) {
  cbs <- read.delim(file=file, stringsAsFactors = F, header = T)
  cbs %<>% set_colnames(c("SampleName", colnames(cbs)[2:length(colnames(cbs))]))
  return(cbs)
  
}