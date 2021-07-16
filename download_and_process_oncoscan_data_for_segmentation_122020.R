
###Read Me####
### This code provide the pipeline for the following processes
#1) download raw data and matrix data from GEO using GEOquery
#2) process to get the ASCAT segment data using EaCoN package
#3) process to get the CBS segment data using EaCoN package

##befeore start, please install the packages 
#require(GEOquery, EaCoN, dplyr, magrittr, etc)
############please install EaCoN using the following command lines
#require(devtools)
#devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
#devtools::install_github("mskcc/facets")
#devtools::install_github("gustaveroussy/EaCoN")
##refer to the EaCoN website (https://github.com/gustaveroussy/EaCoN) for install of other related R packages

########## install and load the library
library(GEOquery)
library(plyr)
library(dplyr)
library(magrittr)
library(EaCoN)

##copy the script to your working directory and load the functions
setwd("/home/rstudio/oncoscan_pipeline") #or your working directory path
source("oncoscan_process_functions.R")

##########set up the GEO series to be processed:
geo.series <- c("GSE125700", "GSE83916", "GSE80806", "GSE107394" )
#gse = "GSE125700"
##create a folder call "data" to store the oncoscan raw data
dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
mainDir = paste0(getwd(), "/data", sep="")

#start the loop for data process, ignore errors, but print the error messages
for(gse in geo.series) {
  print(gse)
  tryCatch({
    #create the series directory and set it as working directory for the processing 
    dir.create(file.path(mainDir, gse), showWarnings = FALSE)
    setwd(file.path(mainDir, gse))
    ##download the matrix with META data 
    gse.mat <- getGEO(paste0(gse), GSEMatrix = TRUE)
    meta.data <- pData(gse.mat[[1]])
    if(length(gse.mat) > 1) {
      for(i in 2:length(gse.mat)) {
        meta.data2 = pData(gse.mat[[i]])
        meta.data = rbind(meta.data, meta.data2)
      }
      
    }
    #head(pData(gse.mat[[1]])[, 1:3])
    ##get the sample and its platform
    meta.data %<>% select(c("geo_accession", "title", "platform_id")) 
    meta.data %<>% mutate(apt.build = ifelse(platform_id %in% "GPL18602", "na33.r4", ifelse(platform_id %in% "GPL21558", "na33.r2", "unknown")))
    #remove sample in the platform can't be processed in EaCoN
    meta.data %<>% filter(apt.build != "unknown")
    ##if all the samples of the series not in platforms (GPL18602 and GPL21558 ), quit the loop and continue with the next series 
    if(nrow(meta.data)==0) {
      print(paste0(gse, " can't be processed as the samples not in GPL18602 and GPL21558 platforms.", sep= " "))
      next
    }
    
    #download raw data and extract the data
    getGEOSuppFiles(gse)
    tarfiles = list.files(path = getwd(), pattern = "*.tar", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)
    ldply(.data = tarfiles, .fun = untar)
    gzfiles = list.files(path = getwd(), pattern = "*CEL.gz", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)
    if(length(gzfiles) == 0 ) { 
      gzfiles = list.files(path = getwd(), patern = "*cel.gz", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)
    }
    ldply(.data = gzfiles, .fun = gunzip)
    
    
    #get the os.pair data and start the batch process
    pair.data <- os.pair()
    colnames(meta.data)[1] <- "SampleName"
    meta.data %<>% left_join(pair.data, by= "SampleName")
    # process the all the samples
    for(i in 1:nrow(meta.data)) {
      OS.Process(ATChannelCel = meta.data[i, 5], GCChannelCel = meta.data[i, 6], samplename = meta.data[i, 1], apt.build = meta.data[i, 4])
    }
    
    Segment.ff.Batch(RDS.files = list.files(path = getwd(), pattern = ".*_processed.RDS$", full.names = TRUE, recursive = TRUE), segmenter = "ASCAT",  recenter = "l2r.median",
                     smooth.k = 5, SER.pen = 20, nrf = 1.0, nthread = 2)
    ASCN.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "SEG\\..*\\.RDS$", full.names = TRUE, recursive = TRUE), nthread = 2)
    ##generate the ASCAT segmentation data #use the data in gamma = 0.55
    cn.files <- list.files(path = getwd(), pattern = "*.gamma0.55.cn$", full.names = TRUE, recursive = TRUE)
    model.files <- list.files(path = getwd(), pattern = "*.gamma0.55_model.txt$", full.names = TRUE, recursive = TRUE)
    cn.list <- lapply(cn.files, read.cn)
    cn.merge <- do.call(rbind,cn.list)
    model.list <- lapply(model.files, read.model)
    model.merge <- do.call(rbind, model.list)
    print(nrow(model.merge)) ##sample number
    cn.seg = cn.merge %>% left_join(model.merge)
    #save the ASCAT segmentation data
    #setwd("/home/bchchen/projects/Jason_group/GEO_BC-GI/R_outputs/CN_segmentation_data")
    write.table(cn.seg, paste0("ASCN_segmenation_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
    ##generate CBS segmentation data
    cbs.files <- list.files(path = getwd(), pattern = "*.NoCut.cbs$", full.names = TRUE, recursive = TRUE)
    cbs.list <-  lapply(cbs.files, read.cbs)
    cbs.merge =  do.call(rbind, cbs.list)
    cbs.seg = cbs.merge 
    write.table(cbs.seg, paste0("CBS_segmenation_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##find out the series have been suceesfully processed or non-processed 
processed.series = c()
unprocessed.series = c()
for(gse in geo.series){ 
  print(gse)
  setwd(file.path(mainDir, gse))
  RDS.files = list.files(path = getwd(), pattern = ".*_processed.RDS$", full.names = TRUE, recursive = TRUE)
  if(length(RDS.files)==0) {
    unprocessed.series = c(unprocessed.series, gse)
  }else{
    processed.series = c(processed.series, gse)
  }
  
}
