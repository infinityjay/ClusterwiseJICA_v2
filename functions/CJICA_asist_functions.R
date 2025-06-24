# Wed Oct 21 09:18:01 2020
# Author: Jeffrey Durieux, MSc

# What: function that allocates subjects into matrix X based on partitioning
# Algorithm step 1

sortX <- function(X, p){
  ClusL <- length( unique(p) )
  NewList <- list()
  
  for(i in 1:ClusL){
    NewList[[i]] <- X[ ,p == i]
  }
  return(NewList)
}


# Wed Oct 21 09:46:37 2020
# Author: Jeffrey Durieux, MSc


# What: apply ICA on sorted data list with matrices
# algorithm step 2

# use VAF > threshold to determine the component number per cluster
cal_nc <- function(X, threshold = 0.8, useChull = TRUE) {
  nc_vector <- NULL
  
  for (data in X) {
    # transfer the dim to (subjects, features)
    data = t(data)
    Maxnc <- min(dim(data))
    # Skip if trivial
    if (Maxnc == 1 || nrow(data) <= 1 || ncol(data) == 0) {
      nc_vector <- c(nc_vector, 1)
      warning("Skipping due to insufficient data dimensions.")
      next
    }
    
    vaf <- NULL
    sstotal = sum(data^2)
    # use PCA here
    # use CHull, complexity = components number, fit = sse or vaf.
    components <- prcomp(data, center = FALSE)
    for (compnr in 1:Maxnc) {
      sse = sum((data - components$x[,1:compnr] %*% t(components$rotation[,1:compnr]))^2)
      vaf = c(vaf, (sstotal- sse) / sstotal)
    }
    if (useChull) {
      comp.fit <- cbind(c(1:Maxnc), vaf)
      comp.fit <- as.data.frame(comp.fit) 
      colnames(comp.fit) <- c("complexity", "fit")
      
      chull_result <- tryCatch({
        CHull(comp.fit, bound = "upper")
      }, error = function(e) {
        warning("CHull failed: ", e$message)
        NULL
      })
      
      # Check if result is a valid CHull object
      if (inherits(chull_result, "CHull") && !is.null(chull_result$Solution$complexity)) {
        nc_vector <- c(nc_vector, chull_result$Solution$complexity)
      } else {
        nc_vector <- c(nc_vector, 1)
      }
    } else {
      nc_vector <- c(nc_vector, min(which(vaf >= threshold)))
    }
  }
  
  return(nc_vector)
}

# change nc to a vector,corresponds to each item in list
# ICAonList <- function(List, nc){
ICAonList <- function(List, nc_vector){
  
  #if else statement does not really matter. ica_adjust and ica::icafast give same results
  # included here only in testphase of coding
  #if(nc == 1){
  Result = NULL
  Result <- mapply(icafast_adjust, List, nc_vector, SIMPLIFY = FALSE)
  
  #  Result <- lapply(List, FUN = ica::icafast, nc = nc) 
  
  S <- lapply(seq_along(Result), function(anom) Result[[anom]]$S)
  M <- lapply(seq_along(Result), function(anom) Result[[anom]]$M)
  
  params <- list(Sr=S,Mr=M)
  
  return(params)
}

# icafast_adjust function
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

XhatsAndLir <- function(X, nClus, Sr, Ahats, Qvec, complex, Vm){
  # todo: N should be different in clusters
  N <- ncol(X)
  V <- Vm
  ss <- matrix(data = NA,nrow = N, ncol = nClus)
  aic <- matrix(data = NA,nrow = N, ncol = nClus)
  penalty <- NA
  for(outer in 1:N){
    
    for(inner in 1:nClus){
      xhat <- Sr[[inner]] %*% Ahats[[inner]][outer,]
      ss[outer,inner] <- sum( (X[,outer] - xhat)^2 )
      I <- nrow(Ahats[[inner]])
      # define the penalty 
      if(complex == 1) {
        penalty <- Qvec[inner] 
      # } else if(complexity == 2) {
      #   penalty <- Qvec[inner]
      } else if(complex == 3) {
        penalty <- Qvec[inner] * V
      } else if(complex == 4) {
        penalty <- Qvec[inner] * V + I
      # } else if(complexity == 5) {
      #   penalty <- Qvec[inner] * V + N
      } else if(complex == 6) {
        penalty <- Qvec[inner] * V + Qvec[inner] * I
      } else if(complex == 7) {
        penalty <- Qvec[inner] * V  + Qvec[inner] * I + I
      } else {
        penalty <- 0
      }
      aic[outer,inner] <- log(ss[outer,inner])*V + 2 * penalty
    }
  }
  
  
  lossvec <- apply(ss, MARGIN = 1, min)
  loss <- sum(lossvec)
  vaf <- ( sum(X^2)-loss) / sum(X^2)
  aicvec <- apply(aic, MARGIN = 1, min)
  aicSum <- sum(aicvec)
  # newp <- apply(ss, MARGIN = 1, which.min)
  newp <- apply(aic, MARGIN = 1, which.min)
  
  out <- list()
  out$newp <- newp
  out$lossvec <- lossvec
  out$loss <- loss
  out$vaf <- vaf
  out$ss <- ss
  out$aic <- aic
  out$aicSum <- aicSum
  out$aicvec <- aicvec
  return(out)
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

# Tucker check function
TuckCheck <- function(S){
  Tucker <- function(X, Y){
    return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
  }
  K <- length(S)
  
  if(K == 5){
    corMod <- mean(c(Tucker(S[[1]][1:2500,], S[[1]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[2]][1:2500,], S[[2]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[3]][1:2500,], S[[3]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[4]][1:2500,], S[[4]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[5]][1:2500,], S[[5]][2501:5000,]) %>% diag() %>% mean)
    )
    
    #combn(1:5,2)
    corClus <- mean(c(Tucker(S[[1]],S[[2]]) %>% diag() %>% mean,
                      Tucker(S[[1]],S[[3]]) %>% diag() %>% mean,
                      Tucker(S[[1]],S[[4]]) %>% diag() %>% mean,
                      Tucker(S[[1]],S[[5]]) %>% diag() %>% mean,
                      
                      Tucker(S[[2]],S[[3]]) %>% diag() %>% mean,
                      Tucker(S[[2]],S[[4]]) %>% diag() %>% mean,
                      Tucker(S[[2]],S[[5]]) %>% diag() %>% mean,
                      
                      Tucker(S[[3]],S[[4]]) %>% diag() %>% mean,
                      Tucker(S[[3]],S[[5]]) %>% diag() %>% mean,
                      
                      Tucker(S[[4]],S[[5]]) %>% diag() %>% mean
    ))
    
  }else if(K == 3){
    
    corMod <- mean(c(Tucker(S[[1]][1:2500,], S[[1]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[2]][1:2500,], S[[2]][2501:5000,]) %>% diag() %>% mean,
                     Tucker(S[[3]][1:2500,], S[[3]][2501:5000,]) %>% diag() %>% mean))
    
    #combn(1:3,2)
    corClus <- mean(c(Tucker(S[[1]],S[[2]])%>% diag() %>% mean,
                      Tucker(S[[1]],S[[3]])%>% diag() %>% mean,
                      Tucker(S[[2]],S[[3]])%>% diag() %>% mean))
  }
  
  return(list(corMod, corClus))
}




