#this file contains the functions to calculate the GI Scores
#load the library
library(dplyr)
library(magrittr)
##function of calculate the GI from the segmenation output
GIS.cal <- function(seg) {
  GI = as.data.frame(unique(seg$SampleName))
  GI %<>% set_colnames(c("SampleName")) %>% distinct() %>% mutate(GIS = "na")
  GI %<>% as.matrix() 
  for(i in 1:nrow(GI)) {
    sample = GI[i, 1]
    data = seg %>% filter(SampleName %in% sample) %>% filter(TCN != 2)
    total = sum(data$Width)
    gi = total/3300000000
    GI[i, 2] = format(round(gi, 2), nsmall = 2)
  }
  GI %<>% as.data.frame(check.names = F)
  return(GI)
}

GID.cal <- function(seg) {
  GI = as.data.frame(unique(seg$SampleName))
  GI %<>% set_colnames(c("SampleName")) %>% distinct() %>% mutate(GID = "na")
  GI %<>% as.matrix() 
  for(i in 1:nrow(GI)) {
    sample = GI[i, 1]
    data = seg %>% filter(SampleName %in% sample)
    data %<>% mutate(cna.sum = select(., CNA, LOH, TAI) %>% rowSums(na.rm = TRUE))
    data %<>% filter(cna.sum > 0)
    A=nrow(data)
    C=length(unique(data$Chr))
    gi = A*A/C
    GI[i, 2] = format(round(gi, 2), nsmall = 2)
  }
  GI %<>% as.data.frame(check.names = F)
  return(GI)
}

GII.cal <- function(seg) {
  GI = as.data.frame(unique(seg$SampleName))
  GI %<>% set_colnames(c("SampleName")) %>% distinct() %>% mutate(GII= "na")
  GI %<>% as.matrix() 
  for(i in 1:nrow(GI)) {
    sample = GI[i, 1]
    data = seg %>% filter(SampleName %in% sample)
    data %<>% mutate(cna.sum = select(., CNA, LOH, TAI) %>% rowSums(na.rm = TRUE))
    data.cnv = data %>% filter(cna.sum > 0)
    gi = nrow(data.cnv)/nrow(data)
    GI[i, 2] = format(round(gi, 2), nsmall = 2)
  }
  GI %<>% as.data.frame(check.names = F)
  return(GI)
}