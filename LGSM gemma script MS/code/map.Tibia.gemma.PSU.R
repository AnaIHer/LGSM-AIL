# Maps QTLs for phenotypes measured in wild-type F1 mice using a linear
# mixed model (GEMMA) that corrects for possible confounding to due
# relatedness.
# Load packages.

require(data.table)
require(qtl)

# SCRIPT PARAMETERS
# -----------------
analysis <-Bones
which.analysis <- "Tibia"
choose.samples <- "ALL"
use.kinship    <- TRUE
results.file   <- "results.gemma.muscles.RData"
seed           <- 1
num.perms      <- 1000

cat("Initiating Analysis ",analysis,".\n",sep="")

# Set location of the data files and the gemma executable.
pheno.filename <- "/Users/anahernandez/Desktop/projects/Project1/Project 1GEMMA_test/gemma/phenotypes.csv"
geno.filename  <- "/Users/anahernandez/Desktop/projects/Project1/Project 1GEMMA_test/gemma/genotypes.csv"
gemmadir       <- "/Users/anahernandez/Desktop/gemma_output"
gemma.exe      <- "/Users/anahernandez/Desktop/P3_GEMMA/gemma"

# Print out the status.
cp0("On ",which.analysis," and ",choose.samples,"\n")
phenotype       <- which.analysis
covariates      <-"location"

# LOAD PHENOTYPE DATA
# -------------------
cat("Loading phenotype data.\n")

# Read in the data, transform it if necessary
pheno <- read.csv(pheno.filename,check.names = T)

# Only include samples in the analysis for which the phenotype and all
# covariates are observed.
cols  <- c(phenotype,covariates)
rows  <- pheno$id
rows  <- which(none.missing.row(pheno[cols]))
pheno <- pheno[rows,]

# Print out some more information  
cat("Including all (",nrow(pheno),") samples in analysis.\n",sep="")

# Convert the sex, time of day, and study columns to binary values. Sex variable is needed for GEMMA
levels(pheno$location)         <- 0:1
levels(pheno$sex)         <- 0:1

# LOAD GENOTYPE DATA
# ------------------
cat("Loading genotype data.\n")
d    <- read.genotypes(geno.filename)
map  <- d$map
geno <- d$geno
geno <- geno[order(rownames(geno)), ]
map  <- cbind(data.frame(snp = rownames(map)),map)
rownames(map) <- NULL
rm(d)

# Align the phenotypes and genotypes.
geno <- geno[match(pheno$id,rownames(geno)),]

# Initialize the random number generator.
set.seed(seed)

r <-vector("list",length(which.analysis))
names(r) <- which.analysis
gwscan.gemma <-r
pve.gemma    <-r

for (i in which.analysis) {
  
  if (use.kinship) {
    
    # Map QTLs using a separate kinship matrix for each chromosome
    out <- run.gemma(i,covariates,pheno,geno,map,gemmadir,gemma.exe)
    gwscan.gemma [[i]] <- out$gwscan
    pve.gemma [[i]]<- out$pve
  } 
}

# PERMUTATIONS
# ----------------------------------------------------------------------------
# Create text files containing the mean genotypes and map information for
# all markers in the format used by GEMMA.

cat("Writing SNP and genotype data to separate files.\n")
write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno,map)
write.gemma.map(paste0(gemmadir,"/map.txt"),map)

# Write out the kinship matrix to file.
cat("Writing identity kinship matrix to file.\n");
write.table(diag(nrow(pheno)),paste0(gemmadir,"/kinship.txt"),
            sep = " ",quote = FALSE,row.names = FALSE,
            col.names = FALSE)

# Initialize the table containing the minimum p-values.
perms.gemma           <- matrix(NA,num.perms,1)

# Define variables used to estimated thresholds

phenotypes<-"Tibia"
  
r <-vector("list",length(phenotypes))
names(r) <- phenotypes
permutations.gemma <-r

# Repeat for each permutation.
for (phenotype in phenotypes) {
  
  for (i in 1:num.perms) {
    
    # Permute the rows of the phenotype table.
    cat("Permutation #",i,'\n',sep="")
    rows <- sample(nrow(pheno))
    
    # Write the phenotype and covariate data to separate text files.
    cat(" * Writing phenotype and covariate data to file.\n")
    write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno[rows,])
    write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),
                           covariates,pheno[rows,])
    
    # Now we are ready to calculate the p-values for all markers using
    # GEMMA.
    cat(" * Computing p-values for all markers using GEMMA.\n")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)
    
    # Load the results of the GEMMA association analysis, and get the
    # minimum p-value.
    gwscan <- read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
    perms.gemma[i] <- max(gwscan$log10p)
  }
  permutations.gemma[[phenotype]]<-perms.gemma
  class(permutations.gemma[[phenotype]]) <- c("scanoneperm","matrix")
}


save(list = c("analysis","which.analysis","covariates","map","gwscan.gemma","pve.gemma","seed",
				"phenotypes","permutations.gemma","num.perms"),
     file = results.file)
     
## ------ OPTIONAL ---------

## change -log10p to LOD scores
# copy the output to a new df
gwscan.gemma.LOD <-gwscan.gemma
# Change lop10p to pvalue.
gwscan.gemma.LOD$Tibia$log10p<- log10p2pval(gwscan.gemma.LOD$Tibia$log10p)
# Convert pvalue to LOD scores.  
gwscan.gemma.LOD$Tibia$log10p<- pval2lod(gwscan.gemma.LOD$Tibia$log10p)
# change the name of the columns
names(gwscan.gemma$Tibia$log10p)<-c("chr","pos","lod")







