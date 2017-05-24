# This file contains miscellaneous functions.

# ----------------------------------------------------------------------
# Combines cat and paste into one function.
cp0 <- function (...)
  cat(paste0(...))
# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (x)
  rowSums(is.na(x)) == 0
# ----------------------------------------------------------------------
# Centers the columns of matrix X so that the entries in each column
# of X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}

# ----------------------------------------------------------------------
# Does the same thing as repmat(A,m,n) in MATLAB.
repmat <- function (A,m,n)
      return(kronecker(matrix(1,m,n),A))

# ----------------------------------------------------------------------
# transform -log10 pvalue to pvalue

log10p2pval <- function(x) {
10^-(x)
}

# ----------------------------------------------------------------------
# transform pvalue to LOD score with one degree of freedom
pval2lod <- function(x) {
(qchisq(x,df=1,lower.tail=F)/(2*log(10)))
}

# ----------------------------------------------------------------------
# transform LOD score with one degree of freedom to pvalue
lod2pval <- function(x) {
pchisq(x*(2*log(10)),df=1,lower.tail=F)
}

# ----------------------------------------------------------------------
# Quantile normalization
# This function transforms the data to a standard normal
# distribution via "quantile normalization". It randomly assigns
# rank to ties.
qt.random.tie <- function (x) {
  x.rank = rank(x,ties.method = "random")
  return(qqnorm(x.rank,plot.it = FALSE)$x)
}




