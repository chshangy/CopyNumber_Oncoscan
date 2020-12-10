## This file shows an example to calculate the GI Scores and HRD scores from the ASCAT segmentation

#load the library
library(dplyr)
library(magrittr)

###HRD functions
source("/home/rstudio/Jason_project/oncoscan_pipeline/HRD_calculation_Functions_multiple_methods_csy.R")
source("/home/rstudio/Jason_project/oncoscan_pipeline/GI_calculation_Functions_multiple_methods_csy.R")
source("/home/rstudio/Jason_project/oncoscan_pipeline/read_segmentation_Functions_csy.R")
  

##########set up the series to be processed:
geo.series <- c("GSE125700", "GSE83916", "GSE80806", "GSE107394" )
mainDir = "/home/rstudio/Jason_project/oncoscan_pipeline/data"

#single gse
#gse="GSE125700"
#i=1
##calculate the
for(i in 1:length(geo.series)) {
  gse=geo.series[i]
  setwd(file.path(mainDir, gse))
  seg.cn <- read.delim(paste0("combined_ASCN_segmenation_", gse, ".txt", sep=""), header = T, stringsAsFactors = F, check.names = F)
  
  #Remove sex chromosomes
  seg.cn %<>% filter(Chr < 23)
  
  #prepare the segmentation in the correct format for the calculation of the TAI, LST and LOH
  seg.cn.hrd = seg.cn %>% select(-Chrom) %>% mutate(contamination = sapply(seg.cn$Celluarity, function(x) {1-x}))
  seg.cn.hrd  %<>% select(SampleName, Chr, Start, End, Probes, TCN, nMajor, nMinor, ploidy.ASCAT, contamination, Celluarity)
  
  chrominfo <- GetChrominfo() # hg19
  
  #get the segment as "TAI"
  tai.seg <- calc.ai(seg.cn.hrd, chrominfo, return.loc = TRUE)
  tai.seg %<>% filter(AI==1) %>% select(SampleName, Chr, Start, End, AI)
  colnames(tai.seg)[which(names(tai.seg)=="AI")] <- "TAI"
  #get tai value
  tai <- calc.ai(seg.cn.hrd, chrominfo)
  tai = data.frame(tai, check.names = F)
  tai %<>% select("Telomeric AI") %>% set_colnames("TAI") %>% mutate(SampleName = rownames(tai))
  
  ######get the segment as LST
  lst.seg <- calc.lst(seg.cn.hrd, chrominfo, return.loc = TRUE)
  lst.seg %<>% filter(LST==1) %>% distinct()
  lst.seg %<>% select(SampleName, Chr, Start, End, LST) %>% distinct()
  ######Calculation of LST
  lst <- calc.lst(seg.cn.hrd, chrominfo)
  lst <- data.frame(lst, check.names = F)
  lst %<>% set_colnames("LST") %>% mutate(SampleName = rownames(lst))
  
  ##get the segment as LOH
  loh.seg <- calc.hrd(seg.cn.hrd, return.loc = TRUE)
  colnames(loh.seg)[which(names(loh.seg)=="HRD breakpoint")] <- "LOH"
  loh.seg %<>% select(SampleName, Chr, Start, End, LOH)
  
  ##Calculation of LOH HRD
  loh <- calc.hrd(seg.cn.hrd)
  loh <- data.frame(loh, check.names = F)
  loh %<>% set_colnames("LOH") %>% mutate(SampleName = rownames(loh))
  
  ##prepare segmentation for calculation of GI scores
  seg.cn$Log2Ratio <- as.numeric(as.character(seg.cn$Log2Ratio))
  #define the CNA as numeric, gain = 1, loss=2, others=0
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
  
  #save the data
  write.table(seg.cn.gi, paste0("Samples_seg_CNA-HRD-LOC_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
  write.table(gi.hrd.all, paste0("Samples_GI_HRD_", gse, ".txt", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
  
}