library(ica)

addError <- function(datablock, error) {
    errorM <- replicate(ncol(datablock), rnorm(nrow(datablock)))
    errorM <- SSequal(errorM, datablock)
    errorlevel <- error / (1 - error)

    res <- datablock + (errorM * sqrt(errorlevel))
    return(res)
}

SSequal <- function(m1, m2) {
    res <- (m1 / sqrt(sum(m1^2)) * sqrt(sum(m2^2))) # R
    return(res)
}

# number of modalities: M = 2; 
# number of features/Voxels: Vm = 2500;
# Cluster number: k = 4 -> R; 
# subjects per cluster: Nk = 10;
# number of components: Qm = 5 -> Q; change to a vector

Simulate_CJICA <- function(Nk, Vm, K, Qm, E, M, type = 1, cor = .5) {
  # type 1: Sk %*% Ak
  # type 2: S %*% Ak
  # type 3: Sk %*% A
  # type 4: is type 1 but with pairwise correlated signals
  
  # Qm is now a vector of length K
  # Check if Qm is a vector with length K
  if (length(Qm) != K) {
    stop("Length of Qm vector must equal K (number of clusters)")
  }
  
  dnames <- c("b")
  P <- rep(1:K, each = Nk)
  
  # generate s
  Slist <- list()
  if (type == 1 | type == 3 | type == 4) {
    for (k in 1:K) {
      S <- matrix(data = NA, nrow = Vm * M, ncol = Qm[k])
      for (q in 1:Qm[k]) {
        s <- numeric()
        for (m in 1:M) {
          s <- c(s, s <- icasamp(
            dname = sample(dnames, size = 1),
            query = "rnd", nsamp = Vm
          ))
        }
        S[, q] <- s
      }
      Slist[[k]] <- S
      
      # For type 4, modify the code to handle variable Qm
      if (type == 4) {
        if (K == 2) {
          r <- matrix(c(1, cor, cor, 1), nrow = 2)
          chol <- chol(r)
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[2])
          for (sig in 1:Qm[1]) {
            # Only process up to the minimum number of signals between clusters
            if (sig <= Qm[2]) {
              ss <- cbind(
                Slist[[1]][, sig],
                s <- icasamp(
                  dname = sample(dnames, size = 1),
                  query = "rnd", nsamp = Vm * 2
                )
              )
              ss <- ss %*% chol
              S2[, sig] <- ss[, 2]
            }
          }
          # If the second cluster has more signals, generate additional ones
          if (Qm[2] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[2]) {
              S2[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          Slist[[2]] <- S2
        } else if (K == 3) {
          r <- matrix(c(
            1, cor, cor - .1,
            cor, 1, cor,
            cor - .1, cor, 1
          ), byrow = T, nrow = 3)
          chol <- chol(r)
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[2])
          S3 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[3])
          
          # Process signals up to the minimum number across clusters
          min_sig1_2 <- min(Qm[1], Qm[2])
          min_sig1_3 <- min(Qm[1], Qm[3])
          
          for (sig in 1:Qm[1]) {
            # Generate correlated signals for those that exist in multiple clusters
            if (sig <= min_sig1_2 || sig <= min_sig1_3) {
              ss <- cbind(
                Slist[[1]][, sig],
                icasamp(
                  dname = sample(dnames, size = 1),
                  query = "rnd", nsamp = Vm * 2
                ),
                icasamp(
                  dname = sample(dnames, size = 1),
                  query = "rnd", nsamp = Vm * 2
                )
              )
              ss <- ss %*% chol
              
              # Only assign if within range for that cluster
              if (sig <= Qm[2]) {
                S2[, sig] <- ss[, 2]
              }
              if (sig <= Qm[3]) {
                S3[, sig] <- ss[, 3]
              }
            }
          }
          
          # Generate additional independent signals for cluster 2 if needed
          if (Qm[2] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[2]) {
              S2[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 3 if needed
          if (Qm[3] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[3]) {
              S3[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          Slist[[2]] <- S2
          Slist[[3]] <- S3
        } else if (K == 4) {
          r <- matrix(
            c(
              1, cor, cor - .1, cor - .2,
              cor, 1, cor, cor - .1,
              cor - .1, cor, 1, cor,
              cor - .2, cor - .1, cor, 1
            ),
            byrow = T, nrow = 4
          )
          chol <- chol(r)
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[2])
          S3 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[3])
          S4 <- matrix(data = NA, nrow = Vm * M, ncol = Qm[4])
          
          # Process signals up to the minimum number across clusters
          for (sig in 1:Qm[1]) {
            # Generate correlated signals for those that exist in multiple clusters
            ss <- cbind(
              Slist[[1]][, sig],
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * 2
              ),
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * 2
              ),
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * 2
              )
            )
            ss <- ss %*% chol
            
            # Only assign if within range for that cluster
            if (sig <= Qm[2]) {
              S2[, sig] <- ss[, 2]
            }
            if (sig <= Qm[3]) {
              S3[, sig] <- ss[, 3]
            }
            if (sig <= Qm[4]) {
              S4[, sig] <- ss[, 4]
            }
          }
          
          # Generate additional independent signals for cluster 2 if needed
          if (Qm[2] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[2]) {
              S2[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 3 if needed
          if (Qm[3] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[3]) {
              S3[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 4 if needed
          if (Qm[4] > Qm[1]) {
            for (sig in (Qm[1]+1):Qm[4]) {
              S4[, sig] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          Slist[[2]] <- S2
          Slist[[3]] <- S3
          Slist[[4]] <- S4
        }
      }
    }
  } else if (type == 2) {
    # For type 2, use the maximum Qm value
    max_Qm <- max(Qm)
    S <- matrix(data = NA, nrow = Vm * M, ncol = max_Qm)
    for (q in 1:max_Qm) {
      s <- numeric()
      for (m in 1:M) {
        s <- c(s, s <- icasamp(
          dname = sample(dnames, size = 1),
          query = "rnd", nsamp = Vm * 2
        ))
      }
      S[, q] <- s
    }
    Slist[[1]] <- S
  }
  
  # generate a
  Alist <- list()
  for (k in 1:K) {
    # Use Qm[k] instead of Qm
    A <- matrix(data = NA, nrow = Nk, ncol = Qm[k])
    for (q in 1:Qm[k]) {
      A[, q] <- runif(n = Nk, min = -2, max = 2)
    }
    Alist[[k]] <- A
  }
  
  # Create X with variable-sized matrices
  if (type == 1 | type == 4) {
    X <- lapply(seq_along(Alist), function(anom) {
      Slist[[anom]] %*% t(Alist[[anom]])
    })
  } else if (type == 2) {
    X <- lapply(seq_along(Alist), function(anom) {
      # Use only the relevant columns for each cluster
      Slist[[1]][, 1:Qm[anom], drop = FALSE] %*% t(Alist[[anom]])
    })
  } else {
    X <- lapply(seq_along(Alist), function(anom) {
      Slist[[anom]] %*% t(Alist[[1]])
    })
  }
  
  X <- do.call(cbind, X)
  Xe <- addError(X, error = E)
  
  out <- list()
  out$Xe <- Xe
  out$X <- X
  out$P <- P
  out$S <- Slist
  out$A <- Alist
  
  return(out)
}
