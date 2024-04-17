#' Gene Count Distribution Plots
#'
#' Retain top x percent of data; results illustrated on plot (red dashed line illustrates threshold)
#'
#' @param threshold percent of data to retain (e.g., if threshold set to 0.15, top 15% (mean expression) of data retained); 0 < threshold < 1
#' @param dat gene count matrix
#' @return A filtered count matrix
#' @export
GeneCountDistributionPlot <- function(dat, threshold) {
  lengths <- dat[, 2]
  ids     <- dat[, 1]
  dat     <- dat[, -c(1:2)]
  
  n.genes <- nrow(dat)
  #rsums   <- rowSums(dat) 
  rsums   <- rowMeans(dat) 
  
  
  # calculate nice expression frequency plot
  evalues    <- log10(rsums)
  max.evalue <- max(evalues)
  intervals  <- 1000 # 10000 bins
  intervals  <- max.evalue/intervals
  
  starts <- seq(0, max.evalue, intervals)
  stops  <- starts
  stops  <- stops[-1]
  stops[length(stops)] <- max.evalue
  starts <- starts[-length(starts)]
  
  DAT <- NULL
  for(n in 1:length(starts)){
    start <- starts[n]
    stop  <- stops[n]
    vals  <- which(evalues >= start & evalues < stop)
    nvals <- length(vals)  
    out   <- cbind(n, nvals, mean(evalues[vals]))
    DAT   <- rbind(DAT, out)
  }
  
  plot(DAT[, 3], DAT[, 2], col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, xlab = "Normalized Expression", ylab = "Frequency", xaxt='n')
  lines(DAT[, 3], DAT[, 2], col="gray50", lwd=1.5)
  
  axis.vals  <- seq(0, max(evalues), 1)
  axis.vals2 <- 10^axis.vals
  axis(side =1, at = axis.vals, labels=axis.vals2)
  
  # show top x%
  rsums2 <- sort(rsums, decreasing = TRUE)
  tops   <- threshold * length(rsums2)
  tops   <- rsums2[1:tops]
  cutoff <- tops[length(tops)]
  abline(v=log10(cutoff), lty=2, col="red")
  
  keep <- which(rsums > cutoff)
  dat  <- dat[keep, ]
  lengths <- lengths[keep]
  
  dat <- cbind(ids[keep], lengths[keep], dat)
  
  return(dat)

}