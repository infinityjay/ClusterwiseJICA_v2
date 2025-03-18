# Mon Oct 26 14:37:26 2020
# Author: Jeffrey Durieux, MSc

# Main function of clusterwise JICA

# X: simulated/real data
# k: cluster size
# nc: component number
# starts: iteration times

# remove nc vector from the params, and calculate it before step 2
ClusterwiseJICA <- function(X, k = 4, starts = 10, threshold = 0.8, scale = T, 
                            rational = NULL, verbose = F){
  
  if(scale == T){
    f1 <- sqrt(5000/sum(X[1:2500,]^2))
    f2 <- sqrt(5000/sum(X[2501:5000,]^2))
    X1 <- f1*X[1:2500,]
    X2 <- f2*X[2501:5000,]
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
      nc <- cal_nc(X = List, threshold)
      
      #### add stop warning over here ##### about nc <= k_n 
      icaparam <- ICAonList(List, nc = nc)
      
      # algo step 3
      Ahh <- Ahats(X = X, icapara = icaparam)
      Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ahh)
      
      # avoid empty clusters
      if( length(unique(Lir$p)) < k ){
        Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp, 
                            SSminVec = Lir$lossvec)
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
