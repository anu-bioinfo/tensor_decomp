#### Additional File 2

#### Summary: this file contains the R-scripts for implementing the tPCA, tWFOBI and tWJADE algorithms as used for the simulation model and application to real data.

#### Libraries needed
library(isva); ### necessary for RMT function which estimates number of significant components
library(tensorBSS); ### necessary for tPCA, tFOBI, tJADE

#### INPUT Arguments
#### Some of the input objects are common to all methods, others are unique to each method. The common ones are described here:

#### data: this is the data object, which contains the multi-way data in two different formats. The A entry of data gives the array or data-tensor format in which the first mode defines the tissue or data-type, the second mode defines the samples and the third mode the features (e.g. CpGs or genes). The L or "list" entry gives the data in list format, where each list entry is the data-matrix for a given tissue or data-type. For tPCA and tICA, only the array format is required.

#### dim: this is the number of significant components of each tissue or data type to search for, and typically we define this to be the RMT estimates for each separate data-or-tissue-type matrix.

#### tpJV.idx: this is an index vector labeling the true positive features associated with a factor which we know drives variation in the data, and which therefore the algorithm should capture. These labels must be in the order of the entries in the 3rd mode of the data-tensor.

#### topN: the number of top-ranked features to select from each inferred component. It must be specified and by default it equals the number of true positives.

#### npermEstDim: this is the number of permutations to use to determine significance of the amount of variance carried by each of the components inferred using tPCA.

DoTPCA <- function(data,dim,npermEstDim=25){

    
    data.a <- data$A;
    nt <- dim(data.a)[1];
    ns <- dim(data.a)[2];
    ng <- dim(data.a)[3];
    dim.v <- c(nt,max(dim));
    tpca.o <- tPCA(data.a,d=dim.v);
    obsE.v <- tpca.o$D[[2]][1:dim.v[2]];

    
    permE.m <- matrix(nrow=npermEstDim,ncol=dim.v[2]);
    for(p in 1:npermEstDim){
     perm.idx <- sample(1:ns,ns,replace=FALSE);
     dataP.a <- data.a;
     dataP.a[1,,] <- data.a[1,perm.idx,];
     tpcaP.o <- tPCA(dataP.a,d=dim.v);
     permE.m[p,] <- tpcaP.o$D[[2]][1:dim.v[2]];
    }
    maxpermE.v <- apply(permE.m,2,max);
    sigCP.idx <- which(obsE.v > maxpermE.v);

    if(length(sigCP.idx)==0){ ### if none significant, still choose top one
      sigCP.idx <- 1;
    }
    
    projS.lm <- list();
    for(t in 1:nt){
     projS.lm[[t]] <- matrix(tpca.o$S[t,sigCP.idx,],nrow=length(sigCP.idx));
    }
    
     return(list(projS.lm=projS.lm,S=tpca.o$S,U=tpca.o$U,sigCP=sigCP.idx,nt=nt,ng=ng));

}

#### tensorial ICA: tWFOBI and tWJADE

DoTICA <- function(data,dim,method=c("FOBI","JADE")){

    dim.v <- c(nt,max(dim));
    selKurt <- c(max(dim),0); ### this ranks all possible components by decreasing kurtosis values
    data.a <- data$A;
    nt <- dim(data.a)[1];
    ns <- dim(data.a)[2];
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
    selCP.m <- selectComponents(tica.o$S,first=selKurt[1],last=selKurt[2]);
    tmp.v <- unlist(strsplit(colnames(selCP.m),split=","));
    sigKCP.idx <- unique(as.numeric(tmp.v[seq(2,length(tmp.v),2)]));

    
    projS.lm <- list();
    for(t in 1:nt){    
      projS.lm[[t]] <- tica.o$S[t,sigKCP.idx,];
    }
    


     return(list(projS.lm=projS.lm,S=tica.o$S,U=tica.o$W,sigCP.idx=sigKCP.idx,nt=nt,ng=ng));

}



