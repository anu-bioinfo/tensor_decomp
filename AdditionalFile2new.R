#### Additional File 2
#### The specific implementations of the tICA methods are made available under a GNU General Public License version 3.

#### Summary: this file contains R-scripts for implementing the tPCA, tWFOBI, tWJADE, JIVE, CCA, sCCA, iCluster and PARAFAC algorithms (as used for the simulation model and application to real omic data). We also provide a function to estimate sensitivity.

#### Libraries needed
library(isva); ### necessary for the RMT function which estimates number of significant components of variation
library(tensorBSS); ### implements tPCA, tFOBI and tJADE
library(r.jive); ### implements JIVE
library(multiway); ### implements PARAFAC
library(PMA); ### implements CCA and sCCA
library(iCluster); ### implements iCluster

#### INPUT Arguments

#### Some of the input objects are common to all methods, others are unique to each method. The common ones are described here:

#### data: this is the data object, which contains the multi-way data in two different formats. The "A" entry of data (data$A) gives the array or data-tensor format, whereas the "L" entry of data (data$L) gives the data in list format. In the former case, and in our specific applications, the first mode defines the tissue, cell or data-type, the second mode defines the samples and the third mode the features (e.g. CpGs or genes). In the latter case, each list entry corresponds to the cell/tissue or data-type and consists of the data-matrix which rows representing features and columns representing samples. 

#### dim: this is a vector which contains the number of significant components of each data matrix to search for, and is typically obtained by applying RMT to each separate data/tissue-type matrix (i.e. to the individual entries of data$L above).

#### dJ: the number of significant components of joint variation across data/tissue types. Typically we define this to be the RMT estimate applied on the joint matrix, obtained by stacking up the matrices in each entry of data$L.

#### tpJV.idx: this is an index vector labeling the true positive features associated with a factor of interest, which we know drives joint variation in the data, and which therefore the algorithm should capture. These indices must be in the order of the entries in the 3rd mode of the data$A object, or alternatively in the same order as the rows of the individual data matrices in data$L.

#### topN: the number of top-ranked features to select from each inferred component. It must be specified and by default it equals the number of true positives.

#### maxiter: The maximum number of iterations used in the algorithms.


#### Additional input:

#### rankJ: An integer giving the number of components of significant joint variation to use when implementing JIVE, if known. If not given, this will be calculated using the chosen method. If the method is "given" then the default is 1.

#### rankA: A vector giving the individual ranks (i.e. number of significant components of individual variation of each data matrix), if known. If not given, this will be calculated using the chosen method. If the method is "given" then the default is rep(1, length(data)).

#### method: A string with the method to use for rank selection. Possible options are "given", "perm", and "bic". The default is "perm". If ranks are known, you should use "given".

#### npermEstDim: this is the number of permutations to use to determine significance of the amount of variance carried by each of the canonical vectors (when using CCA/sCCA).


#### R-SCRIPTS

### Auxiliary function to estimate dimensionality parameters
### Estimate dim and dJ by RMT
EstDim <- function(data){
  data.l <- data$L;
  nt <- length(data.l);
  ### estimate joint variation
  data.m <- data.l[[1]];
  for(i in 1:(nt-1)){
    data.m <- rbind(data.m,data.l[[i+1]]);
  }
  sd.v <- sqrt(apply(data.m,1,var));
  d <- EstDimRMT((data.m - rowMeans(data.m))/sd.v,plot=FALSE)$dim;
  dim <- vector();
  for(j in 1:nt){
    dim[j] <- EstDimRMT( data.l[[j]]-rowMeans(data.l[[j]]), plot = FALSE)$dim;
  }
  return(list(dJ=d, dim=dim));
}

DoTPCA <- function(data,dim){

    data.a <- data$A;
    nt <- dim(data.a)[1];
    ng <- dim(data.a)[3];
    dim.v <- c(nt,max(dim));
    tpca.o <- tPCA(data.a,d=dim.v);

    projS.lm <- list();
    for(t in 1:nt){
     projS.lm[[t]] <- tpca.o$S[t,,];
    }
    
    return(list(projS=projS.lm,S=tpca.o$S,U=tpca.o$U,nt=nt,ng=ng));

}

#### tensorial ICA: tWFOBI and tWJADE

DoTICA <- function(data,dim,method=c("FOBI","JADE")){
    data.a <- data$A;
    nt <- dim(data.a)[1];
    ng <- dim(data.a)[3];
    dim.v <- c(nt,max(dim));
    cdata.a <- tensorCentering(data.a)
    tpca.o <- tPCA(cdata.a,d=dim.v);
    pdata.a <- tensorTransform(cdata.a,t(tpca.o$U[[2]]),2); ## whiten the data

    if(method=="FOBI"){    
      tica.o <- tFOBI(pdata.a);
    }
    else {
      tica.o <- tJADE(pdata.a);
    }
    
    projS.lm <- list();
    for(t in 1:nt){    
      projS.lm[[t]] <- tica.o$S[t,,];
    }

    return(list(projS=projS.lm,S=tica.o$S,U=tica.o$W,nt=nt,ng=ng));

}

DoJIVE <- function(data,rankJ=dJ,rankA=dim,method="given",maxiter=1000){
  
  data.l <- data$L;
  nt <- length(data.l);
  ng <- nrow(data.l[[1]]);
  jive.o <- jive(data.l,rankJ=rankJ,rankA=rankA,method = method, dnames = names(data), conv = "default", maxiter = maxiter, scale = TRUE, center = TRUE, orthIndiv = TRUE, est = TRUE, showProgress=FALSE);
  rankJV <- jive.o$rankJ;
  rankIV.v <- jive.o$rankA;    
  J <- rbind(jive.o$joint[[1]],jive.o$joint[[2]],jive.o$joint[[3]]);
  svd.o <- svd(J);
  jV <- svd.o$v %*% diag(svd.o$d);
  projS.lm <- list();
  for(t in 1:nt){    
    projS.lm[[t]] <- svd.o$u[(1+ng*(t-1)):(t*ng),];
  }
  
  return(list(projS=projS.lm,rankJV=rankJV,rankIV.v=rankIV.v,sJV=jV[,1:rankJV],gJV=svd.o$u[,1:rankJV],nt=nt,ng=ng));

}


DoPARAFAC <- function(data,dim,dJ,maxiter=500){
  data.a <- data$A;
  nt <- dim(data.a)[1];
  ng <- dim(data.a)[3];
  nfac <- sum(dim)-dJ;
  parafac.o <- parafac(data.a,nfac=nfac,nstart=10,maxit=maxiter,ctol=10^-4,parallel=FALSE,output=c("best"));
  projS.m  <- parafac.o$C;
  return(list(projS=projS.m,nfac=nfac,nt=nt,ng=ng));
}


DoSCCA <- function(data,dJ,maxiter=15,npermEstDim=25){
  data.l <- data$L;
  nt <- length(data.l);
  if(nt>2){
    cat("The number of tissue/data type cannot be larger than 2.");
    break;
  }
  ng <- nrow(data.l[[1]]);
  X <- data.l[[1]];
  Z <- data.l[[2]];
  ccaperm.o <- CCA.permute(X,Z,standardize=TRUE,trace=FALSE);
  
  best.idx <- which.max(ccaperm.o$zstats);
  bestpenaltyX <- ccaperm.o$penaltyx[best.idx];
  bestpenaltyZ <- ccaperm.o$penaltyz[best.idx];
  
  cca.o <- CCA(X,Z,K=dJ,niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE);
  
  dP.m <- matrix(nrow=npermEstDim,ncol=max(dim));
  colnames(dP.m) <- paste("CP-",1:max(dim),sep="");
  for(p in 1:npermEstDim){
    perm.idx <- sample(1:nrow(X),nrow(X),replace=FALSE);
    Xp <- X[perm.idx,];
    ccaP.o <- CCA(Xp,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE);
    dP.m[p,] <- ccaP.o$d;
  }
  d.m <- rbind(cca.o$d,dP.m);
  maxdP.v <- apply(dP.m,2,max);
  sigCP.idx <- which(cca.o$d > maxdP.v);
  if(length(sigCP.idx)==0){
    sigCP.idx <- 1;
  }
  X %*% cca.o$u -> u.m;
  Z %*% cca.o$v -> v.m;
  
  return(list(jvX=u.m,jvZ=v.m,u=cca.o$u,v=cca.o$v,sigCP.idx=sigCP.idx,nt=nt,ng=ng));
  
}



DoCCA <- function(data,dJ,maxiter=15,npermEstDim=25){
  data.l <- data$L;
  nt <- length(data.l);
  if(nt>2){
    cat("The number of tissue/data type cannot be larger than 2.");
    break;
  }
  ng <- nrow(data.l[[1]]);
  X <- data.l[[1]];
  Z <- data.l[[2]];
  ### don't estimate best penalty parameters to get the maximal non-sparsity
  bestpenaltyX <- 1;
  bestpenaltyZ <- 1;
  cca.o <- CCA(X,Z,K=dJ,niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE);
  dP.m <- matrix(nrow=npermEstDim,ncol=max(dim));
  colnames(dP.m) <- paste("CP-",1:max(dim),sep="");
  for(p in 1:npermEstDim){
    perm.idx <- sample(1:nrow(X),nrow(X),replace=FALSE);
    Xp <- X[perm.idx,];
    ccaP.o <- CCA(Xp,Z,K=max(dim),niter=maxiter,penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,standardize=TRUE,trace=FALSE)
    dP.m[p,] <- ccaP.o$d
  }
  d.m <- rbind(cca.o$d,dP.m);
  maxdP.v <- apply(dP.m,2,max);
  sigCP.idx <- which(cca.o$d > maxdP.v);
  if(length(sigCP.idx)==0){
    sigCP.idx <- 1;
  }
  X %*% cca.o$u -> u.m;
  Z %*% cca.o$v -> v.m;
  return(list(jvX=u.m,jvZ=v.m,u=cca.o$u,v=cca.o$v,sigCP.idx=sigCP.idx,nt=nt,ng=ng));
}


DoICluster <- function(data,dim,maxiter=100){
  nt <- length(data.l);
  ng <- nrow(data.l[[1]]);
  data.l<- data$L;
  for(i in 1:nt){
    data.l[[i]] <- t(data.l[[i]]);
  }
  fit <- iCluster(datasets = data.l,k = max(dim),lambda = c(0.2,0.2),scalar = F,max.iter = maxiter,epsilon = 1e-3)
  iCluster.l <- list()
  for(i in 1:nt){
    iCluster.l[[i]] <- fit$W[(1+(i-1)*ng):(i*ng),];
  }
  return(list(iCluster.l=iCluster.l,nt=nt,ng=ng,d=max(dim)));
}



EstSE <- function(output.o,tp=tpJV.idx,topN=length(tpJV.idx),method=c("TPCA","FOBI","JADE","JIVE","CCA","sCCA","iCluster","PARAFAC")){
  case <- match(method,c("TPCA","FOBI","JADE","JIVE","CCA","sCCA","iCluster","PARAFAC"))
  switch(case,{### TPCA
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:nrow(output.o$projS[[1]])){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
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
  },{### TWFOBI
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:nrow(output.o$projS[[1]])){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
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
  },{### TWJADE
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:nrow(output.o$projS[[1]])){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS[[t]][cp,]),decreasing=TRUE,index.return=TRUE);
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
  },{### JIVE
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:output.o$rankJV){
      pred.li <- list();
      for(t in 1:output.o$nt){
        tmp.s <- sort(abs(output.o$projS[[t]][,cp]),decreasing=TRUE,index.return=TRUE);
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
  },{### SCCA
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
  },{## CCA
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
  },{### iCluster
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
  },{ ### PARAFAC
    predALL.idx <- vector();
    se.v <- vector(); pv.v <- vector();
    for(cp in 1:output.o$nfac){
      pred.li <- vector();
      
      tmp.s <- sort(abs(output.o$projS[,cp]),decreasing=TRUE,index.return=TRUE);
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


