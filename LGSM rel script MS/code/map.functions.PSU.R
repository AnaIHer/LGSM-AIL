chromosomes <- c(as.character(1:19),"X")

# ----------------------------------------------------------------------
# Convert a factor to a vector of integers.
factor2integer <- function (x)
  match(x,levels(x))

# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s) {
  cat(s,rep("\b",nchar(s)),sep="")
}

# ----------------------------------------------------------------------
# Returns a logical vector such that an entry is TRUE if all the
# entries in the corresponding list element (or column of the table)
# are missing---that is, all the entries are set to NA.
all.missing.col <- function (d)
    sapply(d,function(x) all(is.na(x)))

# READ GENOTYPE FILE FUNCTION
# --------------------------------------------------------------------------------
# Returns a data frame containing the genotype data stored in a CSV
# file. I remove the first two columns of the genotype data frame
# containing the ID and generation, and set the row names for the
# genotype data to the mouse IDs.
read.genotypes <- function (file) {
  geno <- read.csv(file,comment.char="#",as.is = "id",check.names = FALSE)
  row.names(geno) <- geno$id
# 1:1 since I just need the genotypes col 1 = id and row 1=snp id
  geno            <- geno[-(1:1)]
  return(geno)
}

# READ MAP FILE FUNCTION
# --------------------------------------------------------------------------------
# Returns a data frame containing the marker data stored in a CSV file.
# Here I convert the chromosomes to factors manually (as opposed to
# letting the read.csv function handle this) to make sure that the
# chromosomes are ordered properly.
read.map <- function (file, chromosomes)
  within(read.csv(file,comment.char="#",check.names = FALSE,
                  stringsAsFactors = FALSE),
         chr <- factor(chr,chromosomes))

# GENOTYPE to count if genotype file is on ABHN format
# ----------------------------------------------------------------------
# Convert the genotypes to a matrix of allele counts.
# This was not necesary since my file already had genotypes as a count
genotypes2counts <- function (d)
  as.matrix(data.frame(lapply(d,factor2integer),
                       row.names = rownames(d),
                       check.names = FALSE))

# ----------------------------------------------------------------------
# Adjust the map positions (i.e. genetic distances) slightly (by the
# amount j) so that no two markers have the same position.
jitter.gendist <- function (map, j) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)

  # Repeat for each chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers  <- which(map$chr == i)
    n        <- length(markers)

    # Adjust the genetic distances slightly according to scalar j.
    map[markers,"dist"] <- map[markers,"dist"] + cumsum(rep(j,times = n))
  }

  return(map)
}

# ----------------------------------------------------------------------
# Get the genotype probabilities for the specified samples and markers.
subset.genoprob <- function (gp, samples = NULL, markers = NULL) {

  # If the set of markers is the null object, use all the markers.
  if (is.null(markers)) {
    n       <- length(gp$chr)
    markers <- 1:n
  }

  # If the set of samples is the null object, use all the samples.
  if (is.null(samples)) {
    n       <- nrow(gp$pr)
    samples <- 1:n
  }
  
  class(gp) <- c("Pr","list")
  return(within(gp,{
    pr   <- pr[samples,,markers]
    chr  <- chr[markers]
    dist <- dist[markers]
    snp  <- snp[markers]
  }))
}


# ----------------------------------------------------------------------
# Returns the matrix product A*A'.
matrix.square <- function (A)
  return(A %*% t(A))

# ----------------------------------------------------------------------
# Use all the markers to estimate the n x n pairwise relatedness
# matrix, which I define to be 2 times the matrix of kinship
# coefficients (for the definition of kinship coefficients, see Lange,
# 2002). To allow for uncertainty in the genotypes at the markers, I
# compute the expected value of this matrix using the genotype
# probabilities provided as output from the 'genoProb' function in
# QTLRel.
rr.matrix <- function (gp) { 

  # Get the number of samples (n) and the number of markers (p).
  d <- dim(gp$pr)
  n <- d[1]
  p <- d[3]
  
  # Get the genotype probabilities.
  pAA <- gp$pr[,1,]
  pAB <- gp$pr[,2,]
  pBB <- gp$pr[,3,]

  # Get the probability of the genotype AB averaged over all the markers.
  mAB <- matrix(rep(rowMeans(pAB),times = n),n,n)

  # Return the expected value of the pairwise relatedness matrix.
  return(2*(matrix.square(pAA)/p + matrix.square(pBB)/p) 
         + mAB + t(mAB) - matrix.square(pAB)/p)
}

# ----------------------------------------------------------------------
# This function creates a "scanone" object that will be used to store
# QTL mapping results based on the SNP data provided as input. (For
# further info, see the 'scanone' function in the 'qtl' package.) 
empty.scanone <- function (map) {
  d            <- map[c("chr","dist")]
  row.names(d) <- map[["snp"]]
  names(d)     <- c("chr","pos")
  class(d)     <- c("scanone","data.frame")
  return(d) 
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using scanOne from the QTLRel
# library. The function returns a data frame containing the QTL
# mapping results, specifically the LOD scores and the estimated
# additive and dominance effects corresponding to individual markers.
map.cross.rel <- function (pheno, geno, map, gp, vc, phenotype, covariates) {
  out    <- scanOne(pheno[,phenotype],pheno[,covariates],
                    geno,gp,vc,test = "None")
  gwscan <- empty.scanone(map)

  # Get the LOD scores and the proportion of variance explained by each SNP. 
  gwscan$lod <- out$p/(2*log(10))
  gwscan$pve <- out$v/100
  
  # Get the additive and dominance effects of all the SNPs.
  gwscan$add <- sapply(out$parameters,"[[","a")
  gwscan$dom <- sapply(out$parameters,"[[","d")
  
  return(gwscan)
}
