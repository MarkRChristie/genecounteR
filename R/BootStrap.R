#' Boot Strap
#'
#' Calculate significance of each gene via permutation. Lower p-values represent less among-sample variation.
#' Plots "significant" genes as red points
#'
#' @param percent.zeros.retained Genes with <= percent.zeros.retained are retained
#' @param dat gene count matrix
#' @param deltas deltas matrix
#' @param alpha user set threshold for reporting results (e.g., 0.05 returns values <= 0.05)
#' @return An updated deltas matrix (last column has the p-values)
#' @export
BootStrap <- function(dat, deltas, alpha) {
  # plot empirical observations of test statistic
  # plotting deltas[, 3]/deltas[, 2] will find the MOST variable genes (important for future)
  #plot(log(deltas[, 2]), (deltas[, 2]/deltas[, 6]), col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, xlab = "Log Expression", ylab = "E/sum.p.deltaE")
  plot(log(deltas[, 2]), (deltas[, 2]/deltas[, 6])*deltas[, 8], col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, xlab = "log Expression", ylab = "(E/sum.deltaE) * theta")

  
  
  # calculate rank-ordered expression within each gene
  lengths <- dat[, 2]
  ids     <- dat[, 1]
  dat     <- dat[, -c(1:2)]
  
  #n.genes <- nrow(dat)*10
  #n.indiv <- ncol(dat)
  
  head(deltas)
  
  samp1 <- sample(deltas[, 2], 1000000, replace = TRUE)
  samp2 <- sample(deltas[, 6], 1000000, replace = TRUE)
  samp3 <- sample(deltas[, 8], 1000000, replace = TRUE)
  test.stat1 <- (samp1/samp2)*samp3  # idea 1; calculate as sum.deltaE

  emp.stat1  <- (deltas[, 2]/deltas[, 6]) * deltas[, 8] #test statistic from empirical data

  
  OUT <- NULL
  for(n in 1:length(emp.stat1)){
    
    e.stat1 <- emp.stat1[n]
    p.val1  <- length(which(test.stat1 >= e.stat1))
    p.val1  <- p.val1/length(test.stat1)
    
    out <- cbind(n, p.val1)
    OUT <- rbind(OUT, out)
    
  }
  
  deltas <- cbind(deltas, OUT)
  plot(deltas[, 8], deltas[, ncol(deltas)], xlab="theta", ylab="p-value")

  plot(log(deltas[, 2]), (deltas[, 2]/deltas[, 6])*deltas[, 8], col = "dodgerblue", bty = "l", pch = 19, cex = 1.2, xlab = "log Expression", ylab = "(E/sum.deltaE) * theta")
  sigs <- deltas[deltas[, ncol(deltas)] <= alpha,  ]
  points(log(sigs[, 2]), (sigs[, 2]/sigs[, 6])*sigs[, 8], col = "red", pch = 19, cex = 1.2)

  
  return(sigs)
  
}