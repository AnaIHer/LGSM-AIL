# SCRIPT PARAMETERS
#---------------------

analysis <-"Muscles"
num.perm <- 1000
num.markers <-7236

samples <- list(alldata = function (pheno) rep(TRUE,nrow(pheno)))


cat("Initiating Analysis ",analysis,".\n",sep="")

#Specific parameters for the analysis of the four muscles
if (analysis == "Muscles") {
  # QTL MAPPING IN ALL SAMPLES.
  resultsfile   <- "gwscan.muscles.RData"
  keep.mice     <- "all"
  drop.X        <- TRUE
  int.covariate <- NULL
  covariates    <- c("Location","Tibia")
  phenotypes <-c("TA","EDL","Gastroc","Soleus")
} else if (analysis == "Bone") {
  resultfile	   <- "gwscan.bone.RData"
  keep.mice	   <- "all"
  drop.X	   <- TRUE
  int.covariate <-NULL
  covariates	   <- "Location"
  phenotypes<-"Tibia"
} else if (analysis=="BMC and BT"){
  resultfile	   <- "gwscan.bone_t.RData"
  keep.mice	   <- "all"
  drop.X	   <- TRUE
  int.covariate <-NULL
  covariates	   <- c("Location","Tibia")
  phenotypes<-c("appBMC","BT")
}

library(qtl)
capture.output(library(QTLRel))

cat("Loading phenotype data")
pheno <- read.csv("/phenotypes/phenotypes.csv",check.names = T)


# Load genotype and marker data for all the analysis. 
# We do not have the RSnumbers
cat("Loading genotype and marker data")
# read.genotypes function was modifies since I have less arguments
geno 		 <- read.genotypes ("/genotypes/genotypes.csv")
map		 <-read.map("/genotypes/map.csv",chromosomes)
rows         <- which(!is.na(map$snp))
map          <- within(map,snp[rows] <- snp[rows])
# Sort map in ascendent order by chr and position 
map           <- map[with(map,order(map$chr,map$position)),]
names(geno)  <- map$snp

# Drop X chromosome from the mapping 
if (drop.X) {
    
  # Drop the X chromosome from the analysis.
  markers     <- which(map$chr != "X")
  geno        <- geno[,markers]
  map         <- transform(map[markers,],chr = droplevels(chr))
  chromosomes <- levels(map$chr)
} else {

  # Retain the X chromosome, but treat it as an autosome in the QTL
  # mapping.
  chromosomes     <- 1:20
  levels(map$chr) <- chromosomes
}

# TRANSFORM TRAITS
# ----------------
# Convert Location to a binary covariate so that PSU = 0 and UoA = 1
pheno<-transform(pheno,Location=factor2integer(Location)-1)

# Adjust the map possitions slightly so that no two markers have the same position
map <- jitter.gendist(map,1e-6)


# Compute genotype probabilities for all the data
# Genotypes are in 1, 2, 3 format already. Missing values have to be replaced with zero.
G<-(geno)
G[is.na(G)] <- 0
gp	   	<- genoProb(G,map,step = Inf,method = "Haldane",verbose=F)



# QTL MAPPING USING QTLRel
# ------------------------
# Initialize the data structures containing the QTLRel mapping results
# (gwscan.rel) and the parameter estimates for the QTLRel variance 
# components (params).
r          <- vector("list",length(phenotypes))
names(r)   <- phenotypes
gwscan.rel <- r
params     <- r

for (phenotype in phenotypes) {

    # Compute LOD scores for all SNPs across the genome using QTLRel.
    out <- map.cross.rr(pheno,geno,map,gp,phenotype,covariates,int.covariate)
    
    # Get the results of the QTL mapping.
    gwscan.rel[[phenotype]] <- out$gwscan
    params[[phenotype]]     <- out$params
  }

  # Merge the parameter estimates of the variance components for all
  # phenotypes into a single table.
  params[[phenotype]] <- do.call(rbind,params[[phenotype]])
}


# THRESHOLD FOR AUTOSOMES AND X WITH R/QTL
# --------------------------------------------------------------


# LOAD genotypes and phenotype DATA WITH THE R/QTL format cross
data<-read.cross("csv",file="/cross.file/RQTL_data.csv")

#CHECK DATA FOR AUTOSOMES AND X CHR CLASS, PERCENTAGE GENOTYPED.
summary(data)

#ADJUST MARKERS POSSITION, SO NO TWO MARKERS HAVE THE SAME POSITION 
data<-jittermap(data,1e-6)


data<-calc.genoprob(data,step=1,error.prob=0.001)


r          <- vector("list",length(phenotypes))
names(r)   <- phenotypes
gwscan.qtl <- r
perms.qtl  <- r

for (phenotype in phenotypes) {
  
  # Map QTLs using a simple linear regression approach that does not
  # correct for possible confounding to due relatedness.
  cat("Mapping QTLs in muscles using 'qtl'.\n")
  suppressWarnings(gwscan.qtl <-
                     scanone(data,pheno.col = phenotypes,model = "normal",method = "hk",
                             addcovar = data.cov,intcovar = data.intcov,use = "all.obs"))
  
  cat("Permutations Mapping QTLs in MUSCLES using 'qtl'.\n")
  suppressWarnings(perms.qtl <-
                     scanone(data,pheno.col = phenotypes,model = "normal",method = "hk",perm.Xsp=T,
                             addcovar = data.cov,intcovar = data.intcov,use = "all.obs", 
                             n.perm= num.perm, verbose=F))
}



# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(file = resultsfile,
     list = c("analysis","phenotypes","samples","covariates","map",
              "int.covariate","num.markers","keep.mice","drop.X",
              "gwscan.qtl","perms.qtl","gwscan.rel","params"))
