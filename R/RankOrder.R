#' Rank Order
#'
#' This function calculates the rank-order expression (among samples) within a gene 
#' @param dat gene count matrix
#' @return a delta matrix
#' @export
RankOrder <- function(dat) {
  # calculate rank-ordered expression within each gene
  # currently set to calculate sum of deltas; uncomment to calculatue pairwise too
  lengths <- dat[, 2]
  ids     <- dat[, 1]
  dat     <- dat[, -c(1:2)]
  n.genes <- nrow(dat)
  dat     <- as.matrix(dat) # had to add this due to error when adding gene ids
  
  
  OUT <- NULL
  for(n in 1:n.genes){
    #for(n in sample(1:n.genes, 500, replace = FALSE)){
    
    gene <- dat[n, ]
    gene <- gene[order(gene, decreasing = TRUE)] # order gene counts
    var.gene <- var(gene)
    sum.gene <- sum(gene)
    mean.gene<- mean(gene)
    #plot(1:length(gene), gene, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "Counts (millions)")
    
    genea <- gene[-length(gene)]
    geneb <- gene[-1]
    
    deltaE <- genea-geneb # should always be positive or 0 (no negative)
    any(deltaE < 0)       # should always return false
    
    #non ordered plot
    #plot(1:length(deltaE), deltaE, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "delta Expression (millions)")
    
    # rank order deltaE
    deltaE <- deltaE[order(deltaE, decreasing = TRUE)]
    sum.deltaE <- sum(deltaE)
    #plot(1:length(deltaE), deltaE, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "delta Expression (millions)")
    
    ranges <- max(gene)-min(gene)  # equals sum(deltaE); report both
    standardized.deltaE <- deltaE/ranges
    
    sum(standardized.deltaE) == 1 # must equal 1!!!
    
    standardized.deltaE <- c(0, standardized.deltaE)
    #plot(cumsum(standardized.deltaE), pch=19, col="dodgerblue") 
    
    # create null
    null <- rep(1/length(deltaE), length(deltaE))
    null <- c(0, null)
    #points(1:length(standardized.deltaE), cumsum(null), pch=19, col="orange") 
    
    
    theta <- sum(cumsum(standardized.deltaE))/sum(cumsum(null))
    out   <- c(ids[n], sum.gene, mean.gene, var.gene, ranges, sum.deltaE, sum(standardized.deltaE), theta, cumsum(standardized.deltaE))
    
    OUT <- rbind(OUT, out)
    print(n)
    
  }
  
  colnames(OUT)[1:9] <- c("n", "sum.gene", "mean.gene", "var.gene", "ranges", "sum.deltaE", "sum(standardized.deltaE)", "theta", "cumsum(standardized.deltaE)")
  if(sum(OUT[, 7]) == nrow(OUT)) {print("standardization successful")}  # all genes stanrdized to 1
  if(all(OUT[, 5] == OUT[, 6])) {print("deltaE calculations verified")}
  
  hist(OUT[, 8], breaks=30, main="distribution of theta") # distribution of thetas
  
  thetas <- OUT[, -c(1:8)]
  plot(-1000,-1000, col = "dodgerblue", bty = "l", pch = 19, cex = 0.85, xlab = "Comparison", ylab = "Cumulative Sum deltaE", xlim=c(0, ncol(thetas)), ylim = c(0, 1))
  comps <- ncol(thetas)
  lines(1:comps, cumsum(rep(1/comps, comps)), lty=2)
  for(n in 1:nrow(thetas)){
    rows <- thetas[n, ]
    lines(1:length(rows), rows, col="red", lwd=0.05)
  }
  
  #note (OUT[, 2] * OUT[, 6]) is not useful = high theta (good) * large change in E (bad)
  plot(log(OUT[, 2]), (OUT[, 2]/OUT[, 6])*OUT[, 8], col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, xlab = "log Expression", ylab = expression(Omega))
  #plot(log(OUT[, 2]), OUT[, 6], xlab = "log Expression", ylab = "sum.deltaE")
  #plot(log(OUT[, 2]), log(OUT[, 6]), xlab = "log Expression", ylab = "log sum.deltaE")
  #plot(log(OUT[, 2]), OUT[, 8], xlab = "log Expression", ylab = "theta")
  
  # comparisons to variance
  #plot(log(OUT[, 4]), log(OUT[, 5])) # relationship between variance and deltaE
  #plot(log(OUT[, 4]), (OUT[, 8]))    # relationship between variance and theta
  #abline(a=0, b=1)
  
  return(OUT)
  
}