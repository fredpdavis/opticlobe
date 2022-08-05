library(rstan)
library(mclust)
library(loo)
library(matrixStats)
library(txtplot)
library(caret)

inferState.partitionK <- function(n, k, labels, skipThese) {
  # Purpose: Partitions a set of points for k-fold cross-validation
  #
  # returns: array of k-partition assignment
  
  fullList <- 1:n
  if (!missing(skipThese)) {
    fullList <- setdiff(fullList, skipThese)
  } else {
    skipThese <- c()
  }
  
  kparts <- vector(mode="numeric",length=n)
  
  shuffledI <- sample(fullList, length(fullList))
  medSize <- floor(length(fullList) / k)
  for (i in 1:k) {
    t.start <- ((i - 1) * medSize) + 1
    for (j in t.start:(t.start + medSize - 1)) {
      kparts[shuffledI[j]] <- i
    }
  }
  
  if (length(fullList) %% k > 0) {
    for (j in (k * medSize + 1):length(fullList)) {
      kparts[shuffledI[j]] <- k
    }
  }
  
  
  # If points of a particuilar label are all found in same parition, shuffle with
  # another class point to reblanace
  
  keepLabels <- labels[setdiff(fullList, skipThese)]
  
  allGood <- FALSE
  while (!allGood) {
    allGood <- TRUE
    for (curLabel in unique(keepLabels)) {
      parts <- unique(kparts[labels == curLabel])
      if (length(parts) == 1) {
        
        # swap 2 points: 1 from this label, 2 from another label, in another class.
        i <- sample(which(labels == curLabel))[1]
        j <- sample(which(labels != curLabel & kparts != parts[1]))[1]
        part.t <- kparts[i]
        kparts[i] <- kparts[j]
        kparts[j] <- part.t
        
        allGood <- FALSE
      } 
    }
  }
  
  for (i in skipThese) {
    kparts[i] <- k + 1 }
  
  return(kparts)
  
}