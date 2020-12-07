#functions for oncoscan process pinelines

library(dplyr)
library(magrittr)

##function to generate the pairfile data for batch process
os.pair <- function() {
  celfiles <- list.files(path = getwd(), pattern = "*CEL", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)
  if(length(celfiles) == 0 ) { 
    celfiles <- list.files(path = getwd(), pattern = "*cel", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)
  }
  at.gc.files <- split(celfiles, 1:2)
  at.files <- at.gc.files$`1`
  gc.files <- at.gc.files$`2`
  
  os.data <- as.data.frame(at.files, check.names = F) %>% set_colnames("ATChannelCel")
  os.data$ATChannelCel <- as.character(os.data$ATChannelCel)
  os.data$GCChannelCel <- gc.files
  sapply(os.data, class)
  os.data %<>% mutate(SampleName = sapply(.$ATChannelCel, function(x) {
    xx = strsplit(x, "/") %>% unlist() 
    xname = xx[length(xx)]
    sname = strsplit(xname, "_") %>% unlist()
    name = sname[1]
    return(name)
  }
  ))
  return(os.data)
}

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

####functions of read the cbs segmentation file to get the log2R and probes number of the segments 
read.cbs <- function(file) {
  cbs <- read.delim(file=file, stringsAsFactors = F, header = T)
  cbs %<>% set_colnames(c("SampleName", colnames(cbs)[2:length(colnames(cbs))]))
  return(cbs)
  
}