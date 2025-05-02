ILS_CJICA <- function(X, k, nc=c(2,2), useInputQ, VAF = 1, scale = TRUE, complex,Vm = 500 , useChull = T, iter, stepsize, titlevec, rational=NA, verbose = F){
  
  if(scale == TRUE){
    f1 <- sqrt(5000/sum(X[1:2500,]^2))
    f2 <- sqrt(5000/sum(X[2501:5000,]^2))
    X1 <- f1*X[1:2500,]
    X2 <- f2*X[2501:5000,]
    X <- rbind(X1,X2)
  }
  x <-  ClusterwiseJICA_varyQ(X, k = k, starts = 1, nc, useInputQ, VAF, scale, complex, useChull, Vm, rational, verbose)
  x <- x[[1]]
  
  lossvault <- x$Lir$aicSum
  pvault <- x$Lir$newp
  Lirvault <- list(x$Lir)
  Lir <- x$Lir
  
  n <- length(x$Lir$newp)
  it <- 0
  itvault <- 0
  losstrack <- x$Lir$aicSum
  temperature <- 1
  tempstep <- stepsize
  newp <- NULL
  ARItrack <- NULL
  
  #it < iter & temperature < n
  while(it < iter & temperature < n){
    
    it <- it + 1
    cat('Temperature equals :',temperature, fill = TRUE)
    cat('Iteration: ', it, fill = TRUE)
    newp <- rankperb(Lir = Lirvault[[which.min(lossvault)]], nobj = temperature)
    
    repeat{
      loss1 <- Lir$aicSum
      List <- sortX(X, newp)
      ## test fixed Q vector
      
      if (!useInputQ) {
        nc <- cal_nc(X = List, VAF, useChull)
      }
      icaparam <- ICAonList(List, nc = nc)
      Ahh <- Ahats(X = X, icapara = icaparam)
      Lir <- XhatsAndLir(X = X, nClus = k, Sr = icaparam$Sr, Ahats = Ahh, Qvec = nc, complex = complex, Vm = Vm)
      loss2 <- Lir$aicSum
      newp <- Lir$newp
      # do some empty cluster check.
      ###
      if( length(unique(Lir$p)) < k ){
        Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp,
                                        SSminVec = Lir$lossvec)
      }
      
      loss1 - loss2
      losstrack <- c(losstrack,loss2)
      
      ARItrack <- c(ARItrack, round(adjustedRandIndex(titlevec, newp), digits = 4))
      
      if(loss1 - loss2 < 10){
        
        if(sign(loss1-loss2) == -1){
          #increase in loss
          cat('increase')
          temperature <- temperature + tempstep
          break()
        }else if(sign(loss1-loss2) == 0){
          # equal loss
          
          if(loss2 %in% lossvault){
            #if loss2 already in lossvault: increase temp
            temperature <- temperature + tempstep
            break()
          }
          
          if(loss2 > min(lossvault)){
            temperature <- temperature + tempstep
            break()
          }
          
          temperature <- 1
          lossvault <- c(lossvault, loss2)
          pvault <- cbind(pvault, newp)
          Lirvault <- c(Lirvault, list(Lir))
          itvault <- c(itvault,it)
          break()
        }
        
        break()
      } # end if
    }# end repeat
  }#end while
  
  plot(losstrack)
  
  out <- list()
  out$lossvault <- lossvault
  out$pvault <- pvault
  out$Lirvault <- Lirvault
  out$itvault <- itvault
  out$losstrack <- losstrack
  out$ARItrack <- ARItrack

  return(out)
}

scaleprob <- function(x){x/sum(x)}

uij2 <- function(x){sum(x^2)}

ssranking <- function(ss){
  ssscale <- apply(ss, 1, FUN = scaleprob)
  partcoef <- apply(ssscale, MARGIN = 2, FUN = uij2)
  sorted <- sort(partcoef, index.return=TRUE)
  return(sorted$ix)
}

rankperb <- function(Lir, nobj = 1){
  
  k <- sort(unique(Lir$newp))
  rank <- ssranking(Lir$ss)
  
  for(i in 1:nobj){
    kold <- Lir$newp[rank[i]]  
    ids <- which(kold != k)
    newm <- min(Lir$ss[rank[i], ids])
    knew <- which(Lir$ss[rank[i],] == newm)
    Lir$newp[rank[i]] <- knew
  }
  return(Lir$newp)
}