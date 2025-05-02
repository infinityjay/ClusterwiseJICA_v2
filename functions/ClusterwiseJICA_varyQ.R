# Mon Oct 26 14:37:26 2020
# Author: Jeffrey Durieux, MSc

# Main function of clusterwise JICA

# X: simulated/real data
# k: cluster size
# nc: component number
# starts: iteration times

# remove nc vector from the params, and calculate it before step 2
ClusterwiseJICA_varyQ <- function(X, k = 4, starts = 10, nc = c(5,5), useInputQ, VAF = 1, scale = F, 
                                  complex, useChull = T, Vm = 500, rational = NULL, verbose = F){
  if (useInputQ) {
    if (length(nc) != k) {
      stop("Error: nc length not equal to k")
    }
  }
  
  if(scale == T){
    f1 <- sqrt((2*Vm)/sum(X[1:Vm,]^2))
    f2 <- sqrt((2*Vm)/sum(X[Vm+1:2*Vm,]^2))
    X1 <- f1*X[1:Vm,]
    X2 <- f2*X[Vm+1:Vm*2,]
    X <- rbind(X1,X2)
  }
  
  ResultsStarts <- list()
  
  for(start in 1:starts){
    
    totalSS <- sum(X^2)
    lossiter <- totalSS + 1
    iter <- 0
    
    repeat{
      iter <- iter + 1
      if(iter >= 2){
        
        if(verbose == TRUE){
          cat('Start: ', start ,'Iteration ', iter, ' loss value: ', lossiter[iter],'VAF:' ,Lir$vaf ,'\n')    
        }
        
      }else{
        if(verbose == TRUE){
          cat('Start: ', start, 'Iteration ', iter, ' loss value: ', lossiter[iter],'\n')   
        }
      }
      
      
      # algo step 1
      if(iter == 1){
        if(!is.null(rational)){
          # rational will be null in test, ignore the nc->vector problem for now.
          p <- rational
          
          t <- 0
          while( any( table(p)  < nc ) & t < 100 ){
            id <- which(table(p) < nc)  
            id_to_take <- which(table(p) > nc)
            id_to_take <- which(p == id_to_take)
            
            s <- sample(id_to_take, size = 1)
            p[s] <- sample(id, size = 1)
            
            t <- t + 1
          }
          
          
        }else{
          p <- CICA:::clusf(ncol(X), nClus = k)
          
        }
      }else{
        p <- Lir$newp
      }
      List <- sortX(X, p)
      
      # algo step 2
      # calculate the component number per cluster
      if (!useInputQ) {
        nc <- cal_nc(X = List, VAF, useChull)
      }
      
      #### add stop warning over here ##### about nc <= k_n 
      icaparam <- ICAonList(List, nc = nc)
      
      # algo step 3
      Ahh <- Ahats(X = X, icapara = icaparam)
      Lir <- XhatsAndLir(X = X, nClus = k, Sr = icaparam$Sr, Ahats = Ahh, Qvec = nc, complex = complex, Vm)
      
      # avoid empty clusters
      if( length(unique(Lir$newp)) < k ){
        Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp, 
                            SSminVec = Lir$aicvec)
      }
      
      # # avoid clus size lower than nc
      # Lir$newp <- Avoid_nc_N(Lir$newp, Lir$lossvec, nc = nc)
      
      lossiter <- c(lossiter, Lir$loss)
      
      if( lossiter[iter] - lossiter[iter + 1]  < .00001){
        break()
      }
    }
    
    
    out <- list()
    out$p <- Lir$newp
    out$ica <- icaparam
    out$lossiter <- lossiter
    out$Lir <- Lir
    ResultsStarts[[start]] <- out
  }
  return(ResultsStarts)
}
