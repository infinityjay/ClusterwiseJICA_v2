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
        nc_value <- chull_result$Solution$complexity
      } else {
        nc_value <- 1
      }
    } else {
      nc_value <- min(which(vaf >= threshold))
    }
    
    # limit the nc value at 8
    nc_value <- min(nc_value, 8)
    nc_vector <- c(nc_vector, nc_value)
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

XhatsAndLir <- function(X, nClus, Sr, Ahats, Qvec, complex, Vm, useInputQ){
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
  # use aic to constrain
  newp <- apply(aic, MARGIN = 1, which.min)
  if(useInputQ) {
    # original algorithm
    newp <- apply(ss, MARGIN = 1, which.min)
  }
  
  
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

# avoid cluster subjects smaller than 3
SearchSmallClusters <- function(nClus, newcluster, SSminVec, minSize = 3) {
  
  OriCluster <- 1:nClus
  
  # Check how many members are in each cluster
  clusterSizes <- table(factor(newcluster, levels = OriCluster))
  
  # Identify clusters that have fewer than minSize members
  SmallClusters <- as.numeric(names(clusterSizes[clusterSizes < minSize]))
  
  if (length(SmallClusters) == 0) {
    return(newcluster)  # all clusters meet the minimum size
  }
  
  for (smallCluster in SmallClusters) {
    currentSize <- clusterSizes[as.character(smallCluster)]
    needToAdd <- minSize - currentSize
    
    if (needToAdd <= 0) next
    
    cat('Processing small cluster', smallCluster, 'which has', currentSize, 'subjects, needs', needToAdd, 'more.\n')
    
    for (i in 1:needToAdd) {
      # Recompute cluster sizes
      currentClusterSizes <- table(factor(newcluster, levels = OriCluster))
      
      # Find the largest cluster that has > minSize members
      eligibleClusters <- which(currentClusterSizes > minSize)
      
      if (length(eligibleClusters) == 0) {
        cat('SearchSmallClusters: Warning - no eligible clusters to move from\n')
        break
      }
      
      # Select the cluster with the most subjects
      largestCluster <- as.numeric(names(which.max(currentClusterSizes[eligibleClusters])))
      cat('Selected largest cluster:', largestCluster, 'with', currentClusterSizes[as.character(largestCluster)], 'subjects\n')
      
      # Find subjects from the largest cluster
      subjectsInLargest <- which(newcluster == largestCluster)
      
      # Get SS values for those subjects
      SSvals <- SSminVec[subjectsInLargest]
      
      # Select subject with the worst (highest) SS value
      worstIndex <- which.max(SSvals)
      worstSubject <- subjectsInLargest[worstIndex]
      
      # Reassign subject to the small cluster
      newcluster[worstSubject] <- smallCluster
      
      cat('Moved subject', worstSubject, 'with SS value', SSminVec[worstSubject],
          'from cluster', largestCluster, 'to cluster', smallCluster, '\n')
      
      # Update clusterSizes
      clusterSizes <- table(factor(newcluster, levels = OriCluster))
    }
  }
  
  # Final check
  finalSizes <- table(factor(newcluster, levels = OriCluster))
  remainingSmall <- names(finalSizes[finalSizes < minSize])
  
  if (length(remainingSmall) > 0) {
    cat('SearchSmallClusters: Warning - clusters', paste(remainingSmall, collapse = ', '), 
        'still have fewer than', minSize, 'subjects.\n')
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
    corMod <- mean(c(Tucker(S[[1]][1:500,], S[[1]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[2]][1:500,], S[[2]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[3]][1:500,], S[[3]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[4]][1:500,], S[[4]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[5]][1:500,], S[[5]][501:1000,]) %>% diag() %>% mean)
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
    
    corMod <- mean(c(Tucker(S[[1]][1:500,], S[[1]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[2]][1:500,], S[[2]][501:1000,]) %>% diag() %>% mean,
                     Tucker(S[[3]][1:500,], S[[3]][501:1000,]) %>% diag() %>% mean))
    
    #combn(1:3,2)
    corClus <- mean(c(Tucker(S[[1]],S[[2]])%>% diag() %>% mean,
                      Tucker(S[[1]],S[[3]])%>% diag() %>% mean,
                      Tucker(S[[2]],S[[3]])%>% diag() %>% mean))
  }
  
  return(list(corMod, corClus))
}

log_with_time <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat("[", timestamp, "] ", message, "\n", sep = "")
}

# function to select the top sources
select_top_sources <- function(S, n_select = 8, method = "variance") {
  if(ncol(S) <= n_select) {
    log_with_time("Number of sources already within limit")
    return(list(selected_matrix = S, selected_indices = 1:ncol(S)))
  }
  
  log_with_time(paste("Reducing from", ncol(S), "to", n_select, "sources using", method))
  
  if(method == "variance") {
    # Select sources with highest variance
    source_vars <- apply(S, 2, var, na.rm = TRUE)
    selected_idx <- order(source_vars, decreasing = TRUE)[1:n_select]
  } else if(method == "energy") {
    # Select sources with highest energy (sum of squares)
    source_energy <- apply(S^2, 2, sum, na.rm = TRUE)
    selected_idx <- order(source_energy, decreasing = TRUE)[1:n_select]
  } else if(method == "random") {
    # Random selection
    selected_idx <- sample(1:ncol(S), n_select)
  }
  
  log_with_time(paste("Selected sources:", paste(selected_idx, collapse = ", ")))
  
  return(list(
    selected_matrix = S[, selected_idx, drop = FALSE],
    selected_indices = selected_idx
  ))
}

FindOptimalPermutSingle <- function( Sest , Strue, verbose = FALSE, selection_method = "variance")
{
  # code to search the optimal permutation of estimated ICA components for
  # comparing it with simulated components
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  
  #JD: code from Tom, adjusted for matrix vs matrix comparison
  #Sest, Strue (nVoxels x nSources)
  
  log_with_time("Starting FindOptimalPermutSingle function")
  # limit the component number of 2 matrix are same
  n_sources_estimate <- ncol(Sest)
  n_row_est <- nrow(Sest)
  n_sources_true <- ncol(Strue)
  n_row_true <- nrow(Strue)
  log_with_time(paste("dimention of Sest: (", n_row_est,",",n_sources_estimate, "); dimention of Strue: (", n_row_true, ",",n_sources_true, ")"))
  
  if(n_sources_estimate != n_sources_true) {
    log_with_time("Applying source selection to match number of sources")
    
    # Determine the minimum number of sources to keep both matrices comparable
    # Also limit the source number smaller than 8
    n_sources_final <- min(n_sources_estimate, n_sources_true, 8)  
    
    # Reduce Sest if it has more sources than needed
    if(n_sources_estimate > n_sources_final) {
      sest_selected <- select_top_sources(Sest, n_sources_final, selection_method)
      Sest <- sest_selected$selected_matrix
      log_with_time(paste("Reduced Sest from", n_sources_estimate, "to", n_sources_final, "sources"))
    }
    
    # # Reduce Strue if it has more sources than needed
    # if(n_sources_true > n_sources_final) {
    #   strue_selected <- select_top_sources(Strue, n_sources_final, selection_method)
    #   Strue <- strue_selected$selected_matrix
    #   log_with_time(paste("Reduced Strue from", n_sources_true, "to", n_sources_final, "sources"))
    # }
  }
  
  library(gtools)
  N_sources = dim(Sest)[2]
  
  log_with_time(paste("Number of sources:", N_sources))
  
  AllPerms = permutations( n = N_sources , r = N_sources , v = 1:N_sources )
  nPerms = dim(AllPerms)[1]
  
  log_with_time(paste("Total permutations to evaluate:", nPerms))
  
  #Find best permutation
  BestRecov = -9999
  BestPerm = -9999
  log_with_time("Starting permutation search")
  
  for( permtel in 1:nPerms )
  {
    if(verbose == TRUE)
    {
      if( (permtel%%50) == 0)
      {
        print( paste( "perm: " , permtel , "/" , nPerms ) )
      }
    }
    
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tempRecovBlock[1] = mean( abs( diag( Tucker(Strue ,
                                                          Sest[, tp] ) ) ) )
    # niet nodig als het goed is
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  Out$TuckerMatrix = Tucker(Strue , Sest[, BestPerm] )
  
  log_with_time("FindOptimalPermutSingle Function completed")
  return(Out)
}

FindOptimalClusPermut <- function(Pest, Ptrue){
  # find optimal cluster permutation of estimated clustering
  # compared to simulated clustering
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  clus <- length(unique(Pest))
  
  AllPerms = gtools::permutations( n = clus , r = clus)
  nPerms = dim(AllPerms)[1]
  
  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if( (permtel%%50) == 0)
    {
      print( paste( "perm: " , permtel , "/" , nPerms ) )
    }
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tab <- table(Ptrue, Pest)
    
    tempRecovBlock[1] = sum( diag( tab[,tp] ) )
    
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  return(Out)
  
}

Tucker <- function(X, Y){
  return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
}

modRV <- function(X, Y){
  
  if(nrow(X) != nrow(Y)){
    stop('Number of rows of input matrices are not equal')
  }
  
  XXtilde <- ( X %*% t(X) ) - diag (diag( X %*% t(X) ) )
  YYtilde <- ( Y %*% t(Y) ) - diag (diag( Y %*% t(Y) ) )
  
  res <-  ( t(c(XXtilde)) %*% c(YYtilde) ) /
    ( sqrt( ( t(c(XXtilde)) %*% c(XXtilde)) * ( t(c(YYtilde)) %*% c(YYtilde)) ) )
  
  
  return(res)
}

perturbation <- function(p, percentage = 0.1){
  
  clusters <- sort(unique(p))
  sel <- ceiling(length(p) * percentage )
  selected <- sample(1:length(p), size = sel, replace = F)
  
  if(length(selected) == 1){
    # change one cluster
    oriclus <- p[selected]
    newclus <- which(clusters != oriclus)
    
    if(length(newclus) > 1){
      newclus <- sample(newclus, size = 1)
    }
    
    np <- replace(p, selected, newclus)
    
  }else{
    # change multiple clusters
    np <- p
    for(i in 1:length(selected)){
      oriclus <- p[selected[i]]
      newclus <- which(clusters != oriclus)
      
      if(length(newclus) > 1){
        newclus <- sample(newclus, size = 1)
      }
      
      np <- replace(np, selected[i], newclus) # check if this works
    }
  }
  return(np)
}

clusf <- function(nBlocks, nClus) {
  #simplyfied cluster generation function using an equal probability
  clus <- GenerateRandomClustering(nBlocks, nClus, rep(c(1 / nClus), nClus))
  
  return(clus)
}


GenerateRandomClustering <- function(nElement , nClust , Prob = NULL)
{
  ####GenerateRandomClustering = for Random Starts
  
  # Author: Tom F. Wilderjans
  # nElement: number of elements to be clustered
  # nClust: number of clusters
  # Prob (1 x nClust): proportion of elements in each cluster
  
  # Added by Jeffrey Durieux: default Prob = equal cluster prob
  # This done to adjust code later on for potential cluster perbutation?
  
  if(is.null(Prob))
  {
    Prob <- rep(x = (1/nClust) , nClust)
  }
  
  
  BestClust = NULL
  ErrorEncountered = F
  
  if (!(length(Prob) == nClust))
  {
    cat('there should be as much probabilities as clusters')
    ErrorEncountered = T
  }
  
  if ((abs(sum(Prob) - 1) > .000000001) | (any(Prob < 0)))
  {
    cat('probabilities should sum to one (and cannot be negative)')
    ErrorEncountered = T
  }
  
  if (!(any(nClust == 1:nElement)))
  {
    cat("nClus should be a number between 1 and maximal number of datamatrices (length of DataList)")
    ErrorEncountered = T
  }
  
  if (!(ErrorEncountered))
  {
    if (nElement > nClust)
    {
      if (nClust == 1)
      {
        BestClust = rep(1 , times = nElement)
      }
      else
      {
        ProbVV = round(Prob * nElement)
        if (!(sum(ProbVV) == nElement) |
            (any(ProbVV < 1)))
          #not enough elements, or empty clusters
        {
          ProbVV = AdjustProb(ProbVV , nElement)
        }
        
        tempclus = rep(1:length(ProbVV) , ProbVV)
        BestClust = tempclus[sample(1:nElement,size = nElement,replace =
                                      FALSE)]
      }
    }
    else
    {
      BestClust = 1:nClust
    }
  }
  
  if (!(length(unique(BestClust)) == nClust))
  {
    BestClust = NULL
  }
  
  return(BestClust)
}
