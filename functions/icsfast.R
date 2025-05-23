
icafast_adjust <- function (X, nc, center = TRUE, maxit = 100, tol = 1e-06, Rmat = diag(nc), 
                            alg = c("par", "def"), fun = c("logcosh", "exp", "kur"), 
                            alpha = 1) 
{
  X <- as.matrix(X)
  nobs <- nrow(X)
  nvar <- ncol(X)
  nc <- as.integer(nc[1])
  if (nc < 1) {
    stop("Must set nc>=1 component.")
  }
  maxit <- as.integer(maxit[1])
  if (maxit < 1) {
    stop("Must set maxit>=1 iteration.")
  }
  tol <- tol[1]
  if (tol <= 0) {
    stop("Must set ctol>0.")
  }
  if (nc > min(nobs, nvar)) {
    stop("Too many components. Set nc<=min(dim(X)).")
  }
  alpha <- alpha[1]
  if (alpha < 1 | alpha > 2) {
    stop("Must set 'alpha' between 1 and 2.")
  }
  if (nrow(Rmat) != nc | ncol(Rmat) != nc) {
    stop("Input 'Rmat' must be nc-by-nc rotation matrix.")
  }
  if (center) 
    X <- scale(X, scale = FALSE)
  xeig <- eigen(crossprod(X)/nobs, symmetric = TRUE)
  nze <- sum(xeig$val > xeig$val[1] * .Machine$double.eps)
  if (nze < nc) {
    warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
    nc <- nze
    Rmat <- diag(nc)
  }
  Dmat <- sdiag(sqrt(xeig$val[1:nc]))
  
  if(nc == 1){
    Mprt <- Dmat %*% xeig$vectors[,1]  
  }else{
    Mprt <- tcrossprod(Dmat, xeig$vec[, 1:nc])        
  }
  diag(Dmat) <- 1/diag(Dmat)
  Pmat <- xeig$vec[, 1:nc] %*% Dmat
  Xw <- X %*% Pmat
  if (nc == 1L) {
    return(list(S = Xw, M = t(Mprt), W = t(Pmat), Y = Xw, Q = t(Pmat), 
                R = matrix(1), vafs = (sum(Mprt^2) * nobs)/sum(X^2), 
                iter = NA, alg = alg, fun = fun, alpha = alpha))
  }
  if (fun[1] == "kur") {
    fun1d <- function(x) {
      x^3
    }
    fun2d <- function(x) {
      3 * (x^2)
    }
  }
  else if (fun[1] == "exp") {
    fun1d <- function(x) {
      x * exp(-(x^2)/2)
    }
    fun2d <- function(x) {
      exp(-(x^2)/2) * (1 - x^2)
    }
  }
  else {
    fun1d <- function(x) {
      tanh(alpha * x)
    }
    fun2d <- function(x) {
      alpha * (1 - tanh(alpha * x)^2)
    }
  }
  if (alg[1] == "def") {
    myiters <- rep(NA, nc)
    for (j in 1:nc) {
      if (j < 2) {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs, 
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
      else {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        svec <- matrix(0, nc, 1)
        for (k in 1:(j - 1)) {
          svec <- svec + sum(Rmat[, k] * Rmat[, j]) * 
            Rmat[, k]
        }
        Rmat[, j] <- Rmat[, j] - svec
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs, 
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          svec <- matrix(0, nc, 1)
          for (k in 1:(j - 1)) {
            svec <- svec + sum(Rmat[, k] * rnew) * Rmat[, 
                                                        k]
          }
          rnew <- rnew - svec
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
    }
  }
  else {
    rsvd <- svd(Rmat)
    Rmat <- tcrossprod(rsvd$u, rsvd$v)
    iter <- 0
    vtol <- 1
    while (vtol > tol && iter < maxit) {
      smat <- Xw %*% Rmat
      rnew <- crossprod(Xw, fun1d(smat))/nobs
      rnew <- rnew - Rmat %*% sdiag(colMeans(fun2d(smat)))
      rsvd <- svd(rnew)
      rnew <- tcrossprod(rsvd$u, rsvd$v)
      vtol <- 1 - min(abs(colSums(Rmat * rnew)))
      iter <- iter + 1
      Rmat <- rnew
    }
    myiters <- iter
  }
  M <- crossprod(Rmat, Mprt)
  vafs <- rowSums(M^2)
  ix <- sort(vafs, decreasing = TRUE, index.return = TRUE)$ix
  M <- M[ix, ]
  Rmat <- Rmat[, ix]
  vafs <- (vafs[ix] * nobs)/sum(X^2)
  return(list(S = Xw %*% Rmat, M = t(M), W = t(Pmat %*% Rmat), 
              Y = Xw, Q = t(Pmat), R = Rmat, vafs = vafs, iter = myiters, 
              alg = alg[1], fun = fun, alpha = alpha))
}

# Wed Oct 21 11:38:53 2020
# Author: Jeffrey Durieux, MSc

# What: script to compute Ahats 
# algorithm step 3.1

Ahats <- function(X, icapara){
  
  nClus <- length(icapara$Sr)
  N <- ncol(X)
  
  AhatClus <- list()
  for(outer in 1:nClus){
    
    A <- matrix(data = NA, nrow = N, ncol(icapara$Mr[[outer]]))
    
    for(inner in 1:N ){
      A[inner,] <- X[,inner] %*% icapara$Sr[[outer]] %*% 
        NMFN::mpinv( crossprod(icapara$Sr[[outer]]) )
    }
    AhatClus[[outer]] <- A
  }
  
  return(AhatClus)
}

# Wed Oct 21 11:49:35 2020
# Author: Jeffrey Durieux, MSc


# What: computation of Xhats 
# algorithm step 3.2

# Xhat <- Sr %*% Air

XhatsAndLir <- function(X, Sr, Ahats){
  nClus <- length(Sr)
  N <- ncol(X)
  
  ss <- matrix(data = NA,nrow = N, ncol = nClus)
  for(outer in 1:N){
    
    for(inner in 1:nClus){
      xhat <- Sr[[inner]] %*% Ahats[[inner]][outer,]
      ss[outer,inner] <- sum( (X[,outer] - xhat)^2 )
    }
  }
  
  newp <- apply(ss, MARGIN = 1, which.min)
  lossvec <- apply(ss, MARGIN = 1, min)
  loss <- sum(lossvec)
  vaf <- ( sum(X^2)-loss) / sum(X^2)
  
  out <- list()
  out$newp <- newp
  out$lossvec <- lossvec
  out$loss <- loss
  out$vaf <- vaf
  out$ss <- ss
  return(out)
}

# get AIC
AIC <- function(N, J, ss, Q) {
  return(N*J*log(ss) + 2*N*Q)
}


# search the empty clusters
SearchEmptyClusters <- function(nClus, newcluster, SSminVec) {
  
  OriCluster <- 1:nClus
  
  test <- sapply(OriCluster, FUN = '%in%', newcluster)
  
  #test result = no empty clusters so return original newcluster
  
  if ( all( test == TRUE) ){
    newcluster <- newcluster
  }else{
    
    EmptyClusters <- which(test == FALSE)
    singletonnames <- names(which( (table(newcluster)  == 1) == TRUE))
    singletons <- as.numeric(singletonnames)
    id <- which(newcluster %in% singletons == T)
    
    SSminVec[id] <- 0
    
    #worst <- sort( sapply( SSList, FUN = max), decreasing = TRUE)
    worst <- sort( SSminVec, decreasing = TRUE)
    
    #remove worst of singletons, otherwise empties will occur
    
    Index <- sapply( seq_along(EmptyClusters),
                     function(i) FUN = which( SSminVec == worst[i] ) )
    
    # if ties occur in SSminVec
    if( is.null(ncol(Index)) == FALSE ){
      Index <- Index[,1]
      Index <- sample(Index, size = length(EmptyClusters))
    }
    
    for(i in 1:length(Index)){
      newcluster <- replace(newcluster, Index[i], EmptyClusters[i])
      newcluster
    }
    
    
  }# else some emptyclusters
  if( length(unique(newcluster)) != nClus ){
    cat('SearchEmptyCluster, empty occurred')
  }
  return(newcluster)
}

# Thu Oct 22 11:02:22 2020
# Author: Jeffrey Durieux, MSc

# What: avoid nc < N

Avoid_nc_N <- function(newcluster, SSminVec, nc) {
  
  tab <- table(newcluster)
  
  # Check which clusters have fewer than the required components
  toofew <- which(tab < nc)
  
  if (length(toofew) > 0) {
    
    # Number of elements to pick per cluster
    topick <- nc[toofew] - tab[toofew]
    topick <- as.integer(topick)
    
    id_to_pick_from <- which(!(newcluster %in% as.integer(toofew)))
    
    # Sort based on loss value and get index
    ix <- sort.int(SSminVec, decreasing = TRUE, index.return = TRUE)
    id_to_add_to_small_clus <- ix$ix
    
    # Remove IDs that belong to small clusters
    id_to_add_to_small_clus <- id_to_add_to_small_clus[id_to_add_to_small_clus %in% id_to_pick_from]
    
    # Loop over each underpopulated cluster
    for (pp in seq_along(toofew)) {
      num_pick <- topick[pp]
      
      if (num_pick > length(id_to_add_to_small_clus)) {
        stop("Error: Not enough samples to redistribute.")
      }
      
      worst_fit_pick <- id_to_add_to_small_clus[1:num_pick]
      newcluster[worst_fit_pick] <- toofew[pp]
      
      # Remove assigned IDs from available list
      id_to_add_to_small_clus <- id_to_add_to_small_clus[!id_to_add_to_small_clus %in% worst_fit_pick] 
    }
    
    # Final check if all clusters meet the required size
    tabb <- table(newcluster)
    if (any(tabb < nc)) {
      stop("Error: Some clusters still have fewer subjects than required, even after reassignment. Consider reducing the number of components or clusters.")
    }
    
  }
  
  return(newcluster)
}






