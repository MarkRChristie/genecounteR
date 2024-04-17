#' Zero Remover
#'
#' This function removes genes that have a count of zero "0" at every individual 
#' This function also removes genes that have a >= user-specified percentage of samples allowed to have a count of zero; applied for each gene separately
#'
#' @param percent.zeros.retained Genes with <= percent.zeros.retained are retained
#' @param dat gene count matrix
#' @return A filtered count matrix
#' @export
ZeroRemover <- function(dat, percent.zeros.retained) {
  # remove rows that have a 0 for every single sample
  dat2   <- dat[, -c(1, 2)]
  retain <- ncol(dat2)*percent.zeros.retained 
  print(paste("Total number of rows =", nrow(dat2), sep = " "))
 
  rsums   <- rowSums(dat2)
  zeros   <- which(rsums == 0)
  print(paste("Total number of rows/genes where all samples have a zero gene count =", length(zeros), sep = " "))

  counts <- rowSums(dat2 == 0)
  zeros  <- which(counts > retain)
  print(paste("Total number of rows/genes where more than", retain, "samples have a zero gene count =", length(zeros), sep = " "))
  hist(counts, col="skyblue1", main="Counts before zeros removed", xlab="Number of samples that have a count of zero")

  dat    <- dat[-zeros, ]
  counts <- rowSums(dat[, -c(1,2)] == 0)
  hist(counts, col="skyblue1", main="Counts after genes with too many zeros removed", xlab="Number of samples that have a count of zero")
  print(paste("Total number of remaing rows =", nrow(dat), sep = " "))
  
  return(dat)  
}