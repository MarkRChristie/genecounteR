#' Gene Standardization
#' 
#' Standardize across genes with different gene lengths
#' Standardize by residuals between gene length and total expression; linear by default
#'
#' @param dat gene count matrix
#' @return A standarized count matrix
#' @export
GeneStandardization <- function(dat) {
  #scales::alpha() # allow for transparency
  # standardize across genes with different gene lengths
  lengths <- dat[, 2]
  ids     <- dat[, 1]
  dat     <- dat[, -c(1:2)]
  
  # calculate row sums
  rsums   <- rowSums(dat)
  
  # expression length plots
  plot(lengths, rsums/1000000, col = "dodgerblue", bty = "l", pch = 19, cex = 0.85, xlab = "Gene/Contig Length", ylab = "Counts (millions)")
  plot(log(lengths), log(rsums), col =  scales::alpha("dodgerblue", 0.8), bty = "l", pch = 19, cex = 0.85, xlab = "Log Gene/Contig Length", ylab = "Log Counts")
  plot(log(lengths), log(rsums), col =  scales::alpha("dodgerblue", 0.01), bty = "l", pch = 19, cex = 0.85, xlab = "Log Gene/Contig Length", ylab = "Log Counts")
  
  # regression
  tp <- lm(log(rsums) ~ log(lengths))
  summary(tp)
  lines(log(lengths), predict(tp), col = 'blue')
  
  # begin standardize to linear
  y <- log(rsums)
  x <- log(lengths)
  xys <- cbind(1:length(x), x, y)
  new <- predict(tp)
  med.count <- median(xys[, 3])
  xys <- cbind(xys, new, med.count)
  res <- xys[, 3] - xys[, 4]
  
  new.y <- med.count[1]+res # this is correct standardization to rowSums, but next need standardization for each individual
  
  xys <- cbind(xys, new.y)         # add new.y to xys
  xys <- cbind(xys, exp(xys[, 6])) # take exponent to convert back to raw counts
  
  # Begin standardization for each individual
  gene1  <- dat
  m.gene <- apply(gene1, 1, mean) # calculate mean value to standardize to mean; could try median
  diff   <- xys[, 7]-m.gene       # difference between new value and current gene value
  
  gene2  <- gene1 + diff # standardize each gene to each new mean value
  
  
  row.sums2 <- rowSums(gene2) # check to see if standardization succesfull
  plot(x, log(row.sums2), col =  scales::alpha("dodgerblue", 0.01), bty = "l", pch = 19, cex = 0.85, xlab = "Log Gene/Contig Length", ylab = "Log Counts")
  abline(h=median(log(row.sums2)))
  
  # standardize to spline (needs reworked via above format; but probably not needed)
  
  #spline
  #y <- log(rsums)
  #x <- log(lengths)
  #xys <- cbind(x, y)
  #xys <- xys[order(xys[, 1]), ]
  #x <- xys[, 1]
  #y <- xys[, 2]
  #s2  <- smooth.spline(x, y, spar=0.98) # spar is the roughness penalty
  #lines(predict(s2, x), col = "grey70", lwd=2)
  #lines(log(lengths), predict(tp), col = 'blue')
  
  #new <- predict(s2, x)
  #xys[, 4] <- new[[2]]
  #res <- xys[, 3] - xys[, 4]
  #gene4  <- gene2 + res  # add in residuals for each gene
  #row.sums2 <- rowSums(gene4) # check to see if standardization successful
  #plot(x, row.sums2, col = alpha("dodgerblue", 0.01), bty = "l", pch = 19, cex = 0.85, xlab = "Log Gene/Contig Length", ylab = "Log Counts", main = "spline")
  
  
  dat <- cbind(ids, lengths, gene2)
  print("Standardization Successful!")
  

  return(dat)  
}
