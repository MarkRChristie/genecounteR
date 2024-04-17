#' Match Back
#'
#' Matches back significant genes to standardized and plot results
#' @param dat gene count matrix
#' @param delta delta matrix
#' @return returns plots for visualization purposes; no matrices returned
#' @export
MatchBack <- function(delta, dat, ref2) {
  # matches back results to standardized and non-standardized data sets and plot results
  # think about what is important; original counts will not be standardized by diffs in sample depths
  
  # match to standardized data set ---------------------------------------------#
  m1  <- match(deltas[, 1], ref2[, 1])
  red <- ref2[m1, ] #reduced data set
  red <- red[, -c(1:2)]
  red <- as.matrix(red)
  
  # check that matching was correct
  redsums <- rowSums(red)
  test    <- deltas[, 2] == redsums
  if(any(test == FALSE)) {print("FAILED MATCHING")} else {print("MATCHING SUCCESSFUL")}
  
  # average sample of genes for comparison--------------------------------------#
  red2    <- ref2[, -c(1:2)]
  red2    <- as.matrix(red2)
  n.genes <- nrow(red2)
  red2    <- as.matrix(red2) # had to add this due to error when adding gene ids
  
  
  # sample of 20,000 genes
  OUT <- NULL
  for(n in sample(1:n.genes, 20000, replace = FALSE)){
    
    gene <- red2[n, ]
    gene <- gene[order(gene, decreasing = TRUE)] # order gene counts
    #var.gene <- var(gene)
    #plot(1:length(gene), gene, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "Counts (millions)")
    
    #genea <- gene[-length(gene)]
    #geneb <- gene[-1]
    
    #deltaE <- genea-geneb # should always be positive or 0 (no negative)
    #any(deltaE < 0)       # should always return false
    
    #non ordered plot
    #plot(1:length(deltaE), deltaE, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "delta Expression (millions)")
    
    # rank order deltaE
    #deltaE <- deltaE[order(deltaE, decreasing = TRUE)]
    #sum.deltaE <- sum(deltaE)
    #plot(1:length(deltaE), deltaE, col = "dodgerblue", bty = "l", pch = 19, cex = 2, xlab = "Individuals", ylab = "delta Expression (millions)")
    
    #out <- cbind(ids[n], E, var.gene, sd.p.deltaE) # using variance
    OUT <- rbind(OUT, gene)
    
    
  }
  
  red3 <- OUT
  
  red4=rowMeans(red3)
  #red4=mean(log(red4))
  red4=mean(red4)
  
  # calculate percent change/decrease from one column to next 
  OUT2 <- NULL
  for(n in 1:(ncol(red3)-1)){
    #col1 <- red3[, n] + abs(min(red3[, n])) + 0.00000000001
    #col2 <- red3[, n+1] + abs(min(red3[, n+1])) + 0.00000000001
    #asp2=mean(log(col1)-log(col2))
    col1 <- red3[, n] 
    col2 <- red3[, n+1] 
    diff <- (col1-col2)/col1
    diff <- mean(diff)
    OUT2 <- c(OUT2, diff)
  }
    
  red33 <- red4-(cumsum(OUT2)*red4)
  
  #red5 <- colMeans(red3)
  
  # average sample of genes for comparison--------------------------------------#
    for(n in 1:length(m1)){
    gene1 <- red[n, ]
    #gene1 <- log(gene1[order(gene1, decreasing = TRUE)])
    gene1 <- (gene1[order(gene1, decreasing = TRUE)])
    gene2 <- red33 # random set of genes
    #gene2 <- #log(gene2[order(gene2, decreasing = TRUE)])
    
    #ranges <- c(gene1, gene2)
    range2 <- range(gene2)[2]-range(gene2)[1]
    range2 <- range2/2
    
    
    #lines(1:length(gene1), gene1[order(gene1, decreasing = TRUE)], col = "dodgerblue", bty = "l", pch = 19, cex = 1.2)
    #plot(1:length(gene1), gene1, col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, ylim = c(min(ranges)-0.1, max(ranges)+0.1))
    par(mfrow=c(1,2))
    plot(1:length(gene1), gene1, col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, ylim = c(mean(gene1)-range2, mean(gene1)+range2), main = n,  xlab = "Sample", ylab = "Standardized gene count")
    #plot(1:length(gene1), gene1, col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, main = n,  xlab = "Sample", ylab = "Standardized gene count")
    
    
    #gene1 <- red2[n, ] # random set of genes
    #lines(1:length(gene1), gene1[order(gene1, decreasing = TRUE)], col = "red", bty = "l", pch = 19, cex = 1.2)
   
    #plot(1:length(gene2), gene2, col = "red", bty = "l", pch = 19, cex = 1.2, xlab = "Sample", ylab = "Standardized gene count")
    plot(1:length(gene2), gene2, col = "red", bty = "l", pch = 19, cex = 1.2, xlab = "Sample", ylab = "Standardized gene count", main="average gene")
    
    
  }
  
  
  # match to standardized data set log version ---------------------------------------------#
 

  
}
