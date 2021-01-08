# CopyNumber_Oncoscan
Data processing and analysis of oncoscan and oncoscanCNV data

Background

The OncoScan assay is an accepted cancer diagnostic microarray for detection of CNVs, loss of heterozygosity (LOH), and cancer-related somatic mutations,

How to calculate the Genomic Features from the data

Four methods to calculate HRD 

LOH HRD: calculated a total sum of the number of LOH events (segments with only one allele) in each sample. Then, we normalize the value to be in the range 0–1 and termed it LOH HRD

AIL HRD: counted the sum of regions with allelic imbalance, an unequal allele copy number and extension to a sub-telomere without crossing the centromere. 

LST HRD: counted the total number of breakpoints between regions longer than 10 Mb after filtering out regions shorter than 3 Mb. 

HRD: The fourth method was defined as (LOH HRD + AIL HRD + LST HRD)/3,  The HRD value calculated by this method termed “HRD”.

Three methods to calculate Genomic instability

Genomic Instability (GIS): Ratio of the total length of regions with a copy number (CN) other than 2 to a constant of 3.3 × 10^9. So, what I did is just take the sum of all the segment length with CN other than 2 and divide it by 3.3 × 10^9

Genomic Index (GID) : A^2/C, where A is the total number of alterations and C is the number of chromosomes affected by these alterations

Genomic Instability Index (GII) : Fraction of the genome affected by copy number change. Count the number of segments that have been altered (amplification, deletions, LoH, LST, TAI) and then divide it by the number of segments in that sample. 

What the codes do?

1) download raw data and matrix data from GEO using GEOquery

2) process the raw data to get the ASCAT segment data using EaCoN package

3) Calculate the HRD scores (LOH, LST, TAI and mixed-HRD), GI scores (GIS-Genomic instability, GID-Genomic index, GII-Genomic instability index) and median segment size (in kb) from the segmentation data



