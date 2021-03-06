## This file aim to calculate HRD scores, GI Scores and median segment size for each sample from the ASCAT segmentation

#load the library
library(dplyr)
library(magrittr)

###load the functions
setwd("/home/rstudio/oncoscan_pipeline") ##or the directory where the R files locate at
source("HRD_calculation_Functions_multiple_methods_csy.R")
source("GI_calculation_Functions_multiple_methods_csy.R")
source("read_segmentation_Functions_csy.R")

##function of calculate the median segment size in kb from the segmenation output
median.segSize.cal <- function(seg) {
  df = as.data.frame(unique(seg$SampleName))
  df %<>% set_colnames(c("SampleName")) %>% distinct() %>% mutate(median.seg.kb = "na")
  df %<>% as.matrix() 
  for(i in 1:nrow(df)) {
    sample = df[i, 1]
    data = seg %>% filter(SampleName %in% sample)
    median.seg = median(data$Width)/1000
    df[i, 2] = format(round(median.seg, 3), nsmall = 3)
  }
  df %<>% as.data.frame(check.names = F)
  return(df)
}


##########set up the series to be processed:
geo.series <- c("GSE125700", "GSE83916", "GSE80806", "GSE107394" )
mainDir = "/home/rstudio/oncoscan_pipeline/data" #or setup the your path

#single gse
#gse="GSE125700"

##calculate the HRD scores and GI scores for each sample of each gse series
for(i in 1:length(geo.series)) {
  gse=geo.series[i]
  setwd(file.path(mainDir, gse))
  #load the ASCAT segmentation file
  seg.cn <- read.delim(paste0("combined_ASCN_segmenation_", gse, ".txt", sep=""), header = T, stringsAsFactors = F, check.names = F)
  
  #Remove sex chromosomes
  seg.cn %<>% filter(Chr < 23)
  
  #prepare the segmentation in the correct format for the calculation of the TAI, LST and LOH
  seg.cn.hrd = seg.cn %>% select(-Chrom) %>% mutate(contamination = sapply(seg.cn$Celluarity, function(x) {1-x}))
  seg.cn.hrd  %<>% select(SampleName, Chr, Start, End, Probes, TCN, nMajor, nMinor, ploidy.ASCAT, contamination, Celluarity)
  
  chrominfo <- GetChrominfo() # hg19
  
  #get the 'TAI segment' location detail
  tai.seg <- calc.ai(seg.cn.hrd, chrominfo, return.loc = TRUE)
  tai.seg %<>% filter(AI==1) %>% select(SampleName, Chr, Start, End, AI)
  colnames(tai.seg)[which(names(tai.seg)=="AI")] <- "TAI"
  #get tai segment count of each sample
  tai <- calc.ai(seg.cn.hrd, chrominfo)
  tai = data.frame(tai, check.names = F)
  tai %<>% select("Telomeric AI") %>% set_colnames("TAI") %>% mutate(SampleName = rownames(tai))
  
  ######get the 'LST segment' location detail
  lst.seg <- calc.lst(seg.cn.hrd, chrominfo, return.loc = TRUE)
  lst.seg %<>% filter(LST==1) %>% distinct()
  lst.seg %<>% select(SampleName, Chr, Start, End, LST) %>% distinct()
  ######calculate lst segment count of each sample
  lst <- calc.lst(seg.cn.hrd, chrominfo)
  lst <- data.frame(lst, check.names = F)
  lst %<>% set_colnames("LST") %>% mutate(SampleName = rownames(lst))
  
  ##get the 'LOH segment' location detail
  loh.seg <- calc.hrd(seg.cn.hrd, return.loc = TRUE)
  colnames(loh.seg)[which(names(loh.seg)=="HRD breakpoint")] <- "LOH"
  loh.seg %<>% select(SampleName, Chr, Start, End, LOH)
  ##calculate the loh segment count of each sample
  loh <- calc.hrd(seg.cn.hrd)
  loh <- data.frame(loh, check.names = F)
  loh %<>% set_colnames("LOH") %>% mutate(SampleName = rownames(loh))
  
  ##prepare segmentation for calculation of GI scores (GIS=Genomic instability, GID=Genomic index, GII=Genomic instability index)
  seg.cn$Log2Ratio <- as.numeric(as.character(seg.cn$Log2Ratio))
  #define the CNA in terms of gain and loss as numeric, gain = 1, loss=2, others=0
  seg.cn.gi = seg.cn %>% mutate(CNA= sapply(Log2Ratio, function(x){ifelse(x>0.1, 1, ifelse(x < -0.1, 2, 0))}))
  seg.cn.gi %<>% left_join(tai.seg, by=c("SampleName", "Chr", "Start", "End")) %>% left_join(loh.seg, by=c("SampleName", "Chr", "Start", "End")) %>% left_join(lst.seg, by=c("SampleName", "Chr", "Start", "End"))  
  
  seg.cn.gi$CNA <- as.numeric(as.character(seg.cn.gi$CNA))
  seg.cn.gi$TAI <- as.numeric(as.character(seg.cn.gi$TAI))
  seg.cn.gi$LST <- as.numeric(as.character(seg.cn.gi$LST))
  seg.cn.gi$LOH <- as.numeric(as.character(seg.cn.gi$LOH))
  
  ##calculate Genomic instability =GIS
  gis <- GIS.cal(seg.cn.gi)
  ##calculate Genomic index =GID
  gid <- GID.cal(seg.cn.gi)
  gii <- GII.cal(seg.cn.gi)
  gi.all = gis %>% left_join(gid) %>% left_join(gii)
  
  ##
  gi.hrd.all = gi.all %>%left_join(tai) %>%left_join(lst) %>%left_join(loh) %>% mutate(GEO=gse)
  
  ###calculate the median segment size in kb
  median.seg.kb <- median.segSize.cal(seg.cn.gi)
  median.seg.kb %<>% mutate(GEO=gse)
  
  #save the data
  write.table(seg.cn.gi, paste0("Samples_seg_CNA-HRD-LOC_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
  write.table(gi.hrd.all, paste0("Samples_GI_HRD_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
  write.table(median.seg.kb, paste0("Samples_median-seg-size-kb_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)  
}
