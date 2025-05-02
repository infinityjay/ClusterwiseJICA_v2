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

ssq = function(z) sum(z^2)

# number of modalities: M = 2; 
# number of features/Voxels: Vm = 2500;
# Cluster number: k = 4 -> R; 
# subjects per cluster: Nk = 10;
# number of components: Qm = 5 -> Q; change to a vector

Simulate_CJICA <- function(Nk, Vm, K, Qvect, E, M, type = 1, cor = .5, VAF) {
  # type 1: Sk %*% Ak
  # type 2: S %*% Ak
  # type 3: Sk %*% A
  # type 4: is type 1 but with pairwise correlated signals
  
  # Qvect is now a vector of length K
  # Check if Qvect is a vector with length K
  if (length(Qvect) != K) {
    stop("Length of Qvect vector must equal K (number of clusters)")
  }
  
  dnames <- c("b")
  P <- rep(1:K, each = Nk)
  
  # generate s
  Slist <- list()
  if (type == 1 | type == 3 | type == 4) {
    for (k in 1:K) {
      # Changed matrix dimensions to (Vm*M, Nk)
      S <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
      for (n in 1:Nk) {
        s <- numeric()
        for (m in 1:M) {
          s <- c(s, s <- icasamp(
            dname = sample(dnames, size = 1),
            query = "rnd", nsamp = Vm
          ))
        }
        S[, n] <- s
      }
      Slist[[k]] <- S
      
      # For type 4, modify the code to handle variable Qvect
      if (type == 4) {
        if (K == 2) {
          r <- matrix(c(1, cor, cor, 1), nrow = 2)
          chol <- chol(r)
          # Changed dimensions for S2
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          for (n in 1:min(Nk, Qvect[1])) {
            # Only process up to the minimum number of signals between clusters
            ss <- cbind(
              Slist[[1]][, n],
              s <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            )
            ss <- ss %*% chol
            S2[, n] <- ss[, 2]
          }
          # If the second cluster has more columns, generate additional ones
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S2[, n] <- icasamp(
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
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          S3 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          
          # Process signals up to the minimum number across clusters
          min_cols_1_2 <- min(Qvect[1], Nk)
          min_cols_1_3 <- min(Qvect[1], Nk)
          
          for (n in 1:min(Qvect[1], Nk)) {
            # Generate correlated signals for those that exist in multiple clusters
            ss <- cbind(
              Slist[[1]][, n],
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              ),
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            )
            ss <- ss %*% chol
            
            # Only assign if within range for that cluster
            if (n <= Nk) {
              S2[, n] <- ss[, 2]
            }
            if (n <= Nk) {
              S3[, n] <- ss[, 3]
            }
          }
          
          # Generate additional independent signals for cluster 2 if needed
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S2[, n] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 3 if needed
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S3[, n] <- icasamp(
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
          S2 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          S3 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          S4 <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
          
          # Process signals up to the minimum number across clusters
          for (n in 1:min(Qvect[1], Nk)) {
            # Generate correlated signals for those that exist in multiple clusters
            ss <- cbind(
              Slist[[1]][, n],
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              ),
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              ),
              icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            )
            ss <- ss %*% chol
            
            # Only assign if within range for that cluster
            if (n <= Nk) {
              S2[, n] <- ss[, 2]
            }
            if (n <= Nk) {
              S3[, n] <- ss[, 3]
            }
            if (n <= Nk) {
              S4[, n] <- ss[, 4]
            }
          }
          
          # Generate additional independent signals for cluster 2 if needed
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S2[, n] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 3 if needed
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S3[, n] <- icasamp(
                dname = sample(dnames, size = 1),
                query = "rnd", nsamp = Vm * M
              )
            }
          }
          
          # Generate additional independent signals for cluster 4 if needed
          if (Nk > Qvect[1]) {
            for (n in (Qvect[1]+1):Nk) {
              S4[, n] <- icasamp(
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
    max_Qm <- max(Qvect)
    # Changed dimensions for S
    S <- matrix(data = NA, nrow = Vm * M, ncol = Nk)
    for (n in 1:Nk) {
      s <- numeric()
      for (m in 1:M) {
        s <- c(s, s <- icasamp(
          dname = sample(dnames, size = 1),
          query = "rnd", nsamp = Vm
        ))
      }
      S[, n] <- s
    }
    Slist[[1]] <- S
  }

  # generate a
  Alist <- list()
  for (k in 1:K) {
    A <- matrix(data = NA, nrow = Nk, ncol = Nk)
    for (q in 1:Nk) {
      A[, q] <- runif(n = Nk, min = -2, max = 2)
    }
    
    # Ensure Qvect[k] does not exceed Nk
    Qk <- Qvect[k]  
    if (Qk >= Nk) {
      stop("Qvect[k] should be less than Nk")
    }
    
    # Compute scaling factors
    fac1 <- sqrt(VAF / ssq(A[, 1:Qk]))
    fac2 <- sqrt((1 - VAF) / ssq(A[, (Qk + 1):Nk]))
    
    # Apply scaling
    A[, 1:Qk] <- fac1 * A[, 1:Qk]
    A[, (Qk + 1):Nk] <- fac2 * A[, (Qk + 1):Nk]
    
    # Store matrix in list
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
      Slist[[1]][, 1:Qvect[anom], drop = FALSE] %*% t(Alist[[anom]])
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
