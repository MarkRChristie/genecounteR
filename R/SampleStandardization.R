#' Sample Standardization
#'
#' This function standardizes gene counts for samples (e.g., libraries, individuals) with different total numbers of gene counts
#' It also reports the mean, median gene counts per sample and plots the values.
#' @param dat gene count matrix
#' @return A standardized count matrix
#' @export
SampleStandardization <- function(dat) {
  # standardize across samples with different read counts
  dat2 <- dat[, -c(1:2)]
  
  #plot read depth per sample
  depths <- colSums(dat2)
  depths <- depths[order(depths, decreasing = TRUE)]
  mean.depth   <- mean(depths)/1000000
  median.depth <- median(depths)/1000000
  
  plot(1:length(depths), depths/1000000, ylim = c(0, (max(depths/1000000)+5)), col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Samples", ylab = "Counts (millions)", main=c(paste("Mean depth =",round(mean.depth, 2), "million reads per sample", sep=" ")), cex.main = 0.85)
  abline(h=median.depth, lty=2, lwd=0.9)
  
  print(paste("Mean number of reads per samples =", round(mean.depth, 3), "million reads per sample", sep = " "))
  print(paste("Median number of reads per samples =", round(median.depth, 3), "million reads per sample", sep = " "))
  
  
  #standardize by read counts per sample
  depths <- colSums(dat2) # important to call again because re-ordered above
  ncols  <- ncol(dat2)
  ncols == length(depths)
  OUT <- NULL
  for(n in 1:ncols){
    col <- dat2[, n]
    col <- col/depths[n]
    OUT <- cbind(OUT, col) 
  }
  
  # colSums(OUT) should all equal 1!; build in as an explicit test
  if(sum(colSums(OUT)) == ncol(dat2)) {print("Standardization Succesful!")} else {print ("ERROR: Standardization not Succesful")}
  dat2 <- OUT
  
  # calculate raw gene counts (rows) irrespective of length
  gcounts <- rowSums(dat2)
  which(gcounts == 0) # check that all zeros removed
  
  #gcounts  <- gcounts*median(depths)  # standardize gcounts to median count
  dat2  <- dat2*median(depths)     # standardize dat to median count
  dat   <- cbind(dat[, 1:2], dat2)
  
  return(dat)  
}
