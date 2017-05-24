# This file contains the functions to read in genotype data

# ----------------------------------------------------------------------
# Reads the genotypes of the F1 mice derived from the Mouse Diversity
# Array panel. This function returns a list containing two list
# elements: "geno", the n x p matrix of genotypes, where n is the
# number of strains, and p is the number of SNPs; "map", a data frame
# giving the chromosome, base-pair position and alleles for each SNP.
# SKIP the first column , this column have extra information.
read.genotypes <- function (file, skip = 1) {
  chromosomes <- c(1:19,"X")
  genotypes   <- c("A","T","G","C")
  
  # Read the genotype information from the CSV file. Here I'm assuming
  # that the first 25 lines of the file are comments (lines beginning
  # with #).
  d <- fread(file,sep = ",",header = TRUE,stringsAsFactors = FALSE,
             colClasses = list(character = "chr"),na.strings = "NA",
             skip = skip)
  
  # Discard the data.table attributes.
  class(d) <- "data.frame"
  
  # Set the row names to the SNP ids.
  rownames(d) <- d$id
  d           <- d[-1]
  
  # Split the data frame into an n x p matrix of genotypes (where n is
  # the number of strains, and p is the number of SNPs), and the
  # remaining SNP information ("map").
  map  <- d[c("chr","pos","A1","A2")]
  geno <- t(d[-(1:4)])
  storage.mode(geno) <- "double"
  rm(d)
  
  # I convert chromosome to a factor manually so that I can control
  # the order of the chromosomes in the factor. I also convert the
  # allele columns ("A1" and "A2") to factors.
  map <- transform(map,
                   chr = factor(chr,chromosomes),
                   A1  = factor(A1,genotypes),
                   A2  = factor(A2,genotypes))
  
  # Return a list containing two list elements: the genotype matrix
  # ("geno"), and the SNP information ("map").
  return(list(map = map,geno = geno))
}









