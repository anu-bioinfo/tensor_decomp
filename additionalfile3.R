#### Additional File 3

#### Summary: this file contains the R-scripts for implementing the JIVE, CCA, sCCA, iCluster and parafac algorithms as used for the simulation model and application to real data and estimating sensitivity.

#### Libraries needed
library(isva); ### necessary for RMT function which estimates number of significant components
library(r.jive); ### necessary for JIVE
library(multiway); ### necessary for parafac
library(PMA); ### necessary for sCCA
library(impute);
library(iCluster); ### necessary for iCluster
#### INPUT Arguments

#### Some of the input objects are common to all methods, others are unique to each method. The common ones are described here:

#### data: this is the data object, which contains the multi-way data in two different formats. The A entry of data gives the array or data-tensor format in which the first mode defines the tissue or data-type, the second mode defines the samples and the third mode the features (e.g. CpGs or genes). The L or "list" entry gives the data in list format, where each list entry is the data-matrix which rows represent features and columns represent samples for a given tissue or data-type. 

#### dim: this is a vector, which contains the number of significant components of each tissue or data type to search for, and typically we define this to be the RMT estimates for each separate data-or-tissue-type matrix.

#### dJ: the number of significant components of joint matrix. Typically we define this to be the RMT estimates for joint matrix.

#### tpJV.idx: this is an index vector labeling the true positive features associated with a factor which we know drives variation in the data, and which therefore the algorithm should capture. These labels must be in the order of the entries in the 3rd mode of the data-tensor.

#### topN: the number of top-ranked features to select from each inferred component. It must be specified and by default it equals the number of true positives.

#### maxiter:	The maximum number of iterations.

#### npermEstDim: this is the number of permutations to use to determine significance of the amount of variance carried by each of the components.

#### Additional input:

#### rankJ:	An integer giving the joint rank of the data, if known. If not given, this will be calculated using the chosen method. If the method is "given" then the default is 1.

#### rankA:	A vector giving the indvidual ranks of the data, if known. If not given, this will be calculated using the chosen method. If the method is "given" then the default is rep(1, length(data)).

#### method:	A string with the method to use for rank selection. Possible options are "given", "perm", and "bic". The default is "perm". If ranks are known, you should use "given".


#### 

### Estimate dim and dJ by RMT
EstDim <- function(data){
  data.a <- data$A;
  data.l <- data$L;
  nt <- dim(data.a)[1];
  data.m <- data.l[[1]];
  for(i in 1:(nt-1)){
    data.m <- rbind(data.m,data.l[[i+1]]);
  }
  d <- EstDimRMT(data.m = data.m , plot = F)$dim;
  dim <- vector();
  for(j in 1:nt){
    dim[j] <- EstDimRMT(data.m = data.l[[j]] , plot = F)$dim;
  }
  return(list(dJ=d, dim=dim));
}




DoJIVE <- function(data,rankJ=dJ,rankA=dim,method="given",maxiter=1000){
  
  data.a <- data$A;
  data.l <- data$L;
  nt <- dim(data.a)[1];
  ns <- dim(data.a)[2];
  ng <- dim(data.a)[3];
  jive.o <- jive(data.l,rankJ=rankJ,rankA=rankA,method = method, dnames = names(data), conv = "default", maxiter = maxiter, scale = TRUE, center = TRUE, orthIndiv = TRUE, est = TRUE, showProgress=FALSE);
  J <- rbind(jive.o$joint[[1]],jive.o$joint[[2]],jive.o$joint[[3]]);
  svd.o <- svd(J);
  jV <- svd.o$v %*% diag(svd.o$d);
  Z.jive  <- svd.o$u;
  rankJV <- jive.o$rankJ;
  rankIV.v <- jive.o$rankA;
  projS.lm <- list();
  for(t in 1:nt){    
    projS.lm[[t]] <- Z.jive[(1+ng*(t-1)):(t*ng),];
  }
  
  
  return(list(projS.lm=projS.lm,rankJV=rankJV,rankIV.v=rankIV.v,sJV=jV[,1:rankJV],gJV=svd.o$u[,1:rankJV],nt=nt,ng=ng));

}




Doparafac <- function(data,dim,dJ){
  data.a <- data$A;
  data.l <- data$L;
  nt <- dim(data.a)[1];
  ng <- dim(data.a)[3];
  nfac <- sum(dim)-dJ;
  parafac.o <- parafac(data.a,nfac=nfac);
  
  Z.parafac  <- parafac.o$C;
  
  return(list(nfac=nfac,Z=Z.parafac,nt=nt,ng=ng));
  
}


DosCCA <- function(data,dim,maxiter=15,npermEstDim=25){
  data.a <- data$A;
  data.l <- data$L;
  nt <- dim(data.a)[1];
  if(nt>2){
    cat("The number of tissue/data type cannot be larger than 2.")
    break;
  }
  ns <- dim(data.a)[2];
  ng <- dim(data.a)[3];
  X <- data.a[1,,];
  Z <- data.a[2,,];
  ccaperm.o <- CCA.permute(X,Z,standardize=TRUE,trace=FALSE)
  
  best.idx <- which.max(ccaperm.o$zstats)
  bestpenaltyX <- ccaperm.o$penaltyx[best.idx]
  bestpenaltyZ <- ccaperm.o$penaltyz[best.idx]
  
  cca.o <- CCA(X,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE)
  
  dP.m <- matrix(nrow=npermEstDim,ncol=max(dim))
  colnames(dP.m) <- paste("CP-",1:max(dim),sep="")
  for(p in 1:npermEstDim){
    perm.idx <- sample(1:nrow(X),nrow(X),replace=FALSE)
    Xp <- X[perm.idx,]
    ccaP.o <- CCA(Xp,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE)
    dP.m[p,] <- ccaP.o$d
  }
  d.m <- rbind(cca.o$d,dP.m)
  maxdP.v <- apply(dP.m,2,max)
  sigCP.idx <- which(cca.o$d > maxdP.v)
  if(length(sigCP.idx)==0){
    sigCP.idx <- 1
  }
  X %*% cca.o$u -> u.m
  Z %*% cca.o$v -> v.m
  
  return(list(jvX=u.m,jvZ=v.m,u=cca.o$u,v=cca.o$v,sigCP.idx=sigCP.idx,nt=nt,ng=ng));
  
}



DoCCA <- function(data,dim,maxiter=15,npermEstDim=25){
  data.a <- data$A;
  data.l <- data$L;
  nt <- dim(data.a)[1];
  if(nt>2){
    cat("The number of tissue/data type cannot be larger than 2.")
    break;
  }
  ns <- dim(data.a)[2];
  ng <- dim(data.a)[3];
  X <- data.a[1,,];
  Z <- data.a[2,,];
  ### don't estimate best penalty parameters to get the maximal non-sparsity
  bestpenaltyX <- 1
  bestpenaltyZ <- 1
  cca.o <- CCA(X,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE);
  dP.m <- matrix(nrow=npermEstDim,ncol=max(dim))
  colnames(dP.m) <- paste("CP-",1:max(dim),sep="")
  for(p in 1:npermEstDim){
    perm.idx <- sample(1:nrow(X),nrow(X),replace=FALSE)
    Xp <- X[perm.idx,]
    ccaP.o <- CCA(Xp,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE)
    dP.m[p,] <- ccaP.o$d
  }
  d.m <- rbind(cca.o$d,dP.m);
  maxdP.v <- apply(dP.m,2,max);
  sigCP.idx <- which(cca.o$d > maxdP.v)
  if(length(sigCP.idx)==0){
    sigCP.idx <- 1
  }
  X %*% cca.o$u -> u.m
  Z %*% cca.o$v -> v.m
  return(list(jvX=u.m,jvZ=v.m,u=cca.o$u,v=cca.o$v,sigCP.idx=sigCP.idx,nt=nt,ng=ng));
}


DoiCluster <- function(data,dim,maxiter=100){
  data.a <- data$A;
  nt <- dim(data.a)[1];
  ns <- dim(data.a)[2];
  ng <- dim(data.a)[3];
  data.l<- data$L;
  for(i in 1: nt){
    data.l[[i]] <- t(data.l[[i]]);
  }
  fit <- iCluster(datasets = data.l,k = max(dim),lambda = c(0.2,0.2),scalar = F,max.iter = maxiter,epsilon = 1e-3)
  iCluster.l <- list()
  for(i in 1:nt){
    iCluster.l[[i]] <- fit$W[(1+(i-1)*ng):(i*ng),];
  }
  return(list(iCluster.l=iCluster.l,nt=nt,ng=ng,d=max(dim)));
}






EstSE <- function(output.o,tp=tpJV.idx,topN=length(tpJV.idx),method=c("TPCA","FOBI","JADE","JIVE","CCA","sCCA","iCluster","parafac")){
  case <- match(method,c("TPCA","FOBI","JADE","JIVE","CCA","sCCA","iCluster","parafac"))
  switch(case,{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:length(output.o$sigCP.idx)){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS.lm[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
        ### select the topN features and declare these to be the positives for joint variation
        pred.li[[t]] <- tmp.s$ix[1:topN];
      }
      se.v[cp] <- length(intersect(unique(unlist(pred.li)),tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      if(pv < 0.05/length(output.o$sigCP.idx)){
        predALL.idx <- union(predALL.idx,unique(unlist(pred.li)));
      }
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v); 
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:length(output.o$sigCP.idx)){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS.lm[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
        ### select the topN features and declare these to be the positives for joint variation
        pred.li[[t]] <- tmp.s$ix[1:topN];
      }
      ### now estimate sensitivity
      se.v[cp] <- length(intersect(unique(unlist(pred.li)),tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      if(pv < 0.05/length(output.o$sigCP.idx)){
        predALL.idx <- union(predALL.idx,unique(unlist(pred.li)));
      }
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v);
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:length(output.o$sigCP.idx)){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS.lm[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
        ### select the topN features and declare these to be the positives for joint variation
        pred.li[[t]] <- tmp.s$ix[1:topN];
      }
      ### now estimate sensitivity
      se.v[cp] <- length(intersect(unique(unlist(pred.li)),tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      if(pv < 0.05/length(output.o$sigCP.idx)){
        predALL.idx <- union(predALL.idx,unique(unlist(pred.li)));
      }
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v);
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:output.o$rankJV){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS.lm[[t]][,cp]),decreasing=TRUE,index.return=TRUE);
        ### select the topN features and declare these to be the positives for joint variation
        pred.li[[t]] <- tmp.s$ix[1:topN];
      }
      ### now estimate sensitivity
      se.v[cp] <- length(intersect(unique(unlist(pred.li)),tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      
      predALL.idx <- union(predALL.idx,unique(unlist(pred.li)));
      
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v); 
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in output.o$sigCP.idx){
      tmp.s <- sort(abs(output.o$u[,cp]),decreasing=TRUE,index.return=TRUE)
      ### select the topN features and declare these to be the positives for joint variation
      predX.idx <- tmp.s$ix[1:topN]
      tmp.s <- sort(abs(output.o$v[,cp]),decreasing=TRUE,index.return=TRUE)
      ### select the topN features and declare these to be the positives for joint variation
      predZ.idx <- tmp.s$ix[1:topN]
      
      ### now estimate sensitivity
      seX <- length(intersect(predX.idx,tp))/length(tp);
      seZ <- length(intersect(predZ.idx,tp))/length(tp);
      se.v[cp] <- mean(c(seX,seZ))
      predALL.idx <- union(predALL.idx,union(predX.idx,predZ.idx));
    }
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v); 
    return(list(seALL=seALL,seMAX=seMAX,se=se.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in output.o$sigCP.idx){
      tmp.s <- sort(abs(output.o$u[,cp]),decreasing=TRUE,index.return=TRUE)
      ### select the topN features and declare these to be the positives for joint variation
      predX.idx <- tmp.s$ix[1:topN]
      tmp.s <- sort(abs(output.o$v[,cp]),decreasing=TRUE,index.return=TRUE)
      ### select the topN features and declare these to be the positives for joint variation
      predZ.idx <- tmp.s$ix[1:topN]
      
      ### now estimate sensitivity
      seX <- length(intersect(predX.idx,tp))/length(tp);
      seZ <- length(intersect(predZ.idx,tp))/length(tp);
      se.v[cp] <- mean(c(seX,seZ))
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      predALL.idx <- union(predALL.idx,union(predX.idx,predZ.idx));
    }
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v); 
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:output.o$d){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$iCluster.l[[t]][,cp]),decreasing=TRUE,index.return=TRUE);
        ### select the topN features and declare these to be the positives for joint variation
        pred.li[[t]] <- tmp.s$ix[1:topN];
      }
      ### now estimate sensitivity
      se.v[cp] <- length(intersect(unique(unlist(pred.li)),tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v);
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  },{
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:output.o$nfac){
      pred.li <- vector();
      
      tmp.s <- sort(abs(output.o$Z[,cp]),decreasing=TRUE,index.return=TRUE);
      ### select the topN features and declare these to be the positives for joint variation
      pred.li <- tmp.s$ix[1:topN];
      
      ### now estimate sensitivity
      se.v[cp] <- length(intersect(pred.li,tp))/length(tp);
      n <- round(se.v[cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/output.o$ng,lower.tail=FALSE);
      pv.v[cp] <- pv;
      
      predALL.idx <- union(predALL.idx,pred.li);
      
    }
    
    seALL <- length(intersect(predALL.idx,tp))/length(tp);
    seMAX <- max(se.v);  
    return(list(seALL=seALL,seMAX=seMAX,se=se.v,pv=pv.v));
  })
}




