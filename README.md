# LGSM-AIL
This repository contains code and data of the manuscript (submitted/in review), "Replication and discovery of musculoskeletal QTLs in LG/J and SM/J advanced intercross lines"

# Summary of the files in this directory

Folder: LGSM gemma script MS
phenotypes.csv: Muscleskeletal phenotypes data collected on LGSM AIL male mice.
genotypes.csv:  

Folder: LGSM rel script MS
phenotypes.csv: Muscleskeletal phenotypes data collected on LGSM AIL male mice.
RQTLcross.format_data.csv: Cross format file with genotypes and phenotypes. This file format is used by r/qtl.
genotypes.csv:
map.csv:

# Summary of the codes in this directory

Folder: LGSM gemma script MS
map.MuscleM.gemma.PSU.R: R script for mapping muscleskeletal QTLs by chromosome using the software GEMMA 
map.Tibia.gemma.PSU.R: R script for mapping skeletal traits QTLs  using the software GEMMA 
read.gemma.genotypesPSU.R: R function used to run the analyses mapping.
functions.gemma.PSU.R:R script functions needed for mapping QTLs
misc.gemma.PSU.R: extra R functions needed for mapping QTLs

Folder: LGSM rel script MS
map.rel.PSU.R:R script for mapping muscleskeletal QTLs using QTLRel package
map.functions.PSU:R script functions needed for mapping QTLs
map.cross.PSU.R:R script for mapping QTLs using r/qtl package


Credits:
The R code implementing the analysis was based on the work of:
Peter Carbonetto
Department of Human Genetics
University of Chicago
June 2015
and modified by:
Ana I Hernandez Cordero
University of Aberdden
May 2017
