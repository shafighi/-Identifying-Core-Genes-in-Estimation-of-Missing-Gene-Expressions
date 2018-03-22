library(knitr)
library(GEOquery)
library(foreign)
library("corpcor")
library("ggplot2")
library("numDeriv")
library("genlasso")
library("dfoptim")
#library("matrixStats")


akhar <- function(mat){
  mat[nrow(mat),(ncol(mat)-10):ncol(mat)]
}
abaad <- function(mat){
  rbind(nrow(mat),ncol(mat))
  
}


#simlple matrix prediction
prediction <- function(df,cs,rs,all_percent_cols) {
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  print(nrow(df))
  meanS <- NULL;
  for (iteration in 1:9){
    c <- cs[iteration]
    r <- rs[iteration]
    cols <- all_percent_cols[iteration,];
    C <- df[((r+1):nrow(df)),cols]
    X <- df[((r+1):nrow(df)),-cols]
    A <- df[(1:r),cols]
    B <- df[(1:r),-cols]
    main <- rbind(cbind(A,B),cbind(C,X));
    est <- (C %*% pseudoinverse(A)) %*% B;
    err <-sum(abs(est-X)^2);
    err
    meanErr <- sqrt(err / (nrow(X)*ncol(X)))
    meanS <- cbind(meanS,meanErr)
  }
  return(t(meanS) )
}

KNN_impute <- function(df,cs,rs,all_percent_cols,k) {
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  print(nrow(df))
  print(abaad(df))
  meanS <- NULL;
  for (iteration in 1:9){
    c <- cs[iteration]
    r <- rs[iteration]
    cols <- all_percent_cols[iteration,];
    print(iteration)
    print(r)
    print(cols[cols>nrow(df)])
    C <- df[((r+1):nrow(df)),cols]
    X <- df[((r+1):nrow(df)),-cols]
    A <- df[(1:r),cols]
    B <- df[(1:r),-cols]
    NAp <- rep(NA, ncol(X)*nrow(X)); 
    dim(NAp) <- c(nrow(x),ncol(X));
    main <- rbind(cbind(A,B),cbind(C,NAp));
    #est <-  kNNImpute(main, k);
    err <-sum(abs(est-X)^2);
    err
    meanErr <- sqrt(err / (nrow(X)*ncol(X)))
    meanS <- cbind(meanS,meanErr)
  }
  return(t(meanS) )
}

fr <- function(x,A,B,betha) { 
  #return(norm(A%*%x-B[,1],"2") + betha* norm(x,"1"))
  #print(abaad(A))
  #print((B))
  #print(x)
  dim(x) <- c(length(x),1)
  dim(B) <- c(length(B),1)
  #print(betha)
  #print(abaad(x))
  #print(abaad(B))
  #print(betha)
  return (((norm_two(A%*%x-B[,1])^2 ) + betha* norm_one(x)^2))
  #return (norm_two((1/2)*x%*%t(x)%*%t(A)-x%*%B[,1]))
}
gr <- function(fr,x,A,B,betha){
  return(grad(fr,x,A=A,B=B,betha=betha, method="Richardson"))
}
norm_two <- function(x) sqrt(sum(x^2))
norm_one <- function(x) (sum(abs(x)))
optimization<- function() {   ## Rosenbrock Banana function
  betha=0.4;
  eset6 <- GDS2eSet(getGEO("GDS4971"),do.log2=TRUE)
  expr6 <- exprs(eset6)
  expr<-expr6[complete.cases(expr6), ]
  #expr <- data.matrix(expr, rownames.force = NA )
  expr <- expr[,colSums(is.na(expr))<nrow(expr)]
  expr<- (t(expr))  
  c_df<- ncol(expr)
  r_df<-nrow(expr)
  all_percent_cols<-NULL
  cs<-NULL
  rs<-NULL
  df<-expr
  for(percents in 1:9){
    percent <- percents*10;# darsad e missing 
    temp <- integer(c_df);
    cs <- cbind(cs,floor(c_df-sqrt((1-percent/100)*c_df)));
    rs <- cbind(rs,floor(r_df-sqrt(1-percent/100)*r_df));
    #print(c_df)
    #print(c[percents])
    temp[1:cs[percents]] <- sample.int(c_df,cs[percents]);
    all_percent_cols<-rbind(all_percent_cols,temp);
  }
  r<-rs[1]
  c<-cs[1]
  cols<-all_percent_cols[1,]
  C <- df[(r+1):r_df,cols];
  X <- df[(r+1):r_df,-cols];
  A <- df[c(1:r),cols];
  B <- df[c(1:r),-cols]; 
  main <- rbind(cbind(A,B),cbind(C,X));
  alpha <- matrix(runif((c*(c_df-c)),0.00001, 0.0001));
  dim(alpha) <- c(c,(c_df-c)) 
  
  falpha <- matrix(runif(c,-0.001,0.001));
  dim(falpha) <- c(c,1) 
  falpha <- pseudoinverse(A) %*% B[,1]+falpha
  bound <- rep(c(-1,1),c)
  dim(bound) <- c(2,c)
  #o<-optimize(fr,initial ,A=A,B=B,betha=betha)
  #opt4<-optim(falpha,fr,NULL,A=A,B=B,betha=betha,control = list(reltol = 0.0000000001,maxit = 20000))
  #lasso<-fusedlasso(B[,1], A, diag(ncol(A)), gamma = 0, approx = FALSE, maxsteps = 2000,
  #          minlam = 0, rtol = 1e-07, btol = 1e-07, eps = 1e-4,
  #         verbose = FALSE)
  opt <- hjkb(falpha, fr, A=A,B=B[,1],betha=betha, lower = bound[1,], upper = bound[2,])
  op2<-A%*%opt$par-B[,1]
  err<-C%*%opt$par-X[,1]
  errB <-sum(abs(err)^2);
  errB
  meanErrB <- sqrt(errB /nrow(X));
  meanErrB
  
  estB <- (C %*% pseudoinverse(A)) %*% B;
  #main <- rbind(cbind(AA,BB),cbind(CC,estB));
  eB <-sum(abs(estB[,1]-X[,1])^2);
  eB
  meanEB <- sqrt(eB /nrow(X));
  meanEB
  
}


#Itrative way of matrix prediction


iterative <- function(df,cs,rs,all_percent_cols){
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]  
  #parameter setting
  r_df<-nrow(df);
  c_df<-ncol(df);
  all_iter_means<-NULL;
  
  for(percents in 1:9){
    percent <- percents*10;# darsad e missing 
    c <- cs[percents]
    r <- rs[percents]
    cols <- all_percent_cols[percents,];
    
    #initialization
    meanS <- NULL;
    minS <- NULL;
    maxS <- NULL;
    percent=percent/10;
    
    #first whole estimation
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    #print(abaad(XX))
    #print(percents)
    estB <- (CC %*% pseudoinverse(AA)) %*% BB;
    main <- rbind(cbind(AA,BB),cbind(CC,estB));
    errB <-sum(abs(estB-XX)^2);
    errB
    meanErrB <- sqrt(errB / (nrow(XX)*ncol(XX)));
    meanS <- cbind(meanS,meanErrB);
    minS <- cbind(minS,min(estB-XX));
    maxS <- cbind(maxS,max(estB-XX));
    bestMain<- main;
    bestMean<-meanErrB;
    bestMeanS<-meanErrB;
    for(iteri in 1:1000) {
      
      
      #ms<- meanS
      
      #for(iteri in 1:100) {
      #obtimization by row fraction [bellow columns]
      exr <- sample.int(r_df-r-4,1)+r;#floor((r_df+r)/2);
      exr
      part1r<-exr-r;
      part2r<-r_df-exr;
      
      #main<-bestMain;
      Cr <- main[(exr+1):r_df,1:c];
      Xr <- main[(exr+1):r_df,(c+1):c_df];
      Ar <- main[1:exr,1:c];
      Br <- main[1:exr,(c+1):c_df];  
      estr1 <- (Cr %*% pseudoinverse(Ar)) %*% Br;
      main <- rbind(cbind(Ar,Br),cbind(Cr,estr1));
      estr <- main[c((r+1):r_df),(c+1):c_df];
      err <-sum(abs(estr-XX)^2);  
      meanErr <- sqrt(err / (nrow(XX)*ncol(XX)));
      meanS <- cbind(meanS,meanErr);
      if(meanErr<bestMean){# this condition can be on the min and max
        bestMeanS <- cbind(bestMeanS,meanErr);
        minS <- cbind(minS,min(estr-XX));
        maxS <- cbind(maxS,max(estr-XX));
        bestMain<-main;
        bestMean<-meanErr;
      }
      
      #moving double estimated block
      #main <- bestMain[c(1:r,(exr+1):r_df,(r+1):exr),];
      main <- main[c(1:r,(exr+1):r_df,(r+1):exr),];
      
      Crn <- main[(r+part2r+1):r_df,1:c];
      Xrn <- main[(r+part2r+1):r_df,(c+1):c_df];
      Arn <- main[1:(r+part2r),1:c];
      Brn <- main[1:(r+part2r),(c+1):c_df];  
      estrn <- (Crn %*% pseudoinverse(Arn)) %*% Brn;
      
      #moving back the block and make estimated block in order it should be
      main <- cbind(rbind(AA,Cr,Crn),rbind(BB,estrn,estr1))#AA,B,...,,,,CC,est,Xp
      #est <- main[c((r+1):nrow(main)),(c+1):ncol(main)];
      
      est <- rbind(estrn,estr1);
      #print("third")
      print(abaad(est))
      print(abaad(XX))
      print(abaad(estrn))
      print(abaad(estr1))
      err <-sum(abs(est-XX)^2);  
      meanErr <- sqrt(err / (nrow(XX)*ncol(XX)));
      
      meanS <- cbind(meanS,meanErr);
      
      if(meanErr<bestMean){# this condition can be on the min and max
        bestMeanS <- cbind(bestMeanS,meanErr);
        minS <- cbind(minS,min(est-XX));
        maxS <- cbind(maxS,max(est-XX));
        bestMain<-main;
        bestMean<-meanErr;
      }
    }
    all_iter_means <- rbind(all_iter_means,meanS);
  }
  return(all_iter_means)
}


#clustering

clustering <- function(expr3795,clusters,cs,rs,all_percent_cols){
  
  df<-data.frame(t(expr3795),check.names = FALSE);
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  #dfc <- scale(df[-1]) 
  dft<-data.frame((expr3795),check.names = FALSE);
  dft <- data.matrix(dft, rownames.force = NA )
  dft <- dft[,colSums(is.na(dft))<nrow(dft)]
  
  res <- kmeans(dft, centers=clusters, nstart=10)
  #res1 <- kcca(dft,3,dist=cov())
  #res <- kcca2df(res1)
  
  #initialization
  meanS <- NULL;
  minS <- NULL;
  maxS <- NULL;
  meanSP <- NULL;
  minSP <- NULL;
  maxSP <- NULL;
  
  
  r_df<-nrow(df);
  c_df<-ncol(df);
  for(i in 1:9){
    c <- cs[i]
    r <- rs[i]
    cols <- all_percent_cols[i,];
    
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    estB <- (CC %*% pseudoinverse(AA)) %*% BB;
    main <- rbind(cbind(AA,BB),cbind(CC,estB));
    errB <-sum(abs(estB-XX)^2);
    errB
    meanErrB <- sqrt(errB / (nrow(XX)*ncol(XX)));
    meanS <- rbind(meanS,meanErrB);
    minS <- rbind(minS,min(estB-XX));
    maxS <- rbind(maxS,max(estB-XX));
    minP<-NULL;
    maxP<-NULL;
    errP<-NULL;
    
    for(clus in 1:clusters){
      truecol1<-names(which(res$cluster==clus));
      idx_first <- match(truecol1,colnames(CC));
      idx_second<-match(truecol1,colnames(BB));
      f<- idx_first[is.na(idx_first)==FALSE];
      s<- idx_second[is.na(idx_second)==FALSE];
      if(length(f)<2){
        CP<-CC;
        AP<-AA;
        BP<-BB[,s];
        XP1<-XX[,s];
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }else if(length(s)==0){
        meanErrP1<-0;
      }else{
        CP<-CC[,f];
        AP<-AA[,f];
        BP<-BB[,s];
        XP1<-XX[,s];
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }
    }
    
    meanErrP <- sqrt(sum(errP) / (nrow(XX)*ncol(XX)));
    meanSP <- rbind(meanSP,meanErrP);
    minSP <- rbind(minSP,min(minP));
    maxSP <- rbind(maxSP,max(maxP));
  }
  return(meanSP)
}


#regression
regres <-function(expr3795,cs,rs,all_percent_cols,per,rep,iter,betha,first_gamma, last_gamma,first_alpha, last_alpha){
  
  df<-data.frame(t(expr3795),check.names = FALSE);
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  
  r_df<-nrow(df);
  c_df<-ncol(df);
  all_percents<-NULL;
  for(percent in 1:per){
    c <- cs[percents]
    r <- rs[percents]
    cols <- all_percent_cols[percents,];
    
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    estB <- (CC %*% pseudoinverse(AA)) %*% BB;
    main <- rbind(cbind(AA,BB),cbind(CC,estB));
    errB <-sum(abs(estB-XX)^2);
    errB
    meanErrB <- sqrt(errB / (nrow(XX)*ncol(XX)));
    all_meanS<- NULL;
    
    for(j in 1:rep){
      #initialize
      alpha <- matrix(runif((c*(c_df-c)),first_alpha, last_alpha));
      dim(alpha) <- c(c,(c_df-c))
      #alpha <- pseudoinverse(AA) %*% BB
      alpha_pre <- alpha
      grad_pre <- 0;
      #gamma <- matrix(runif(((c_df-c)*(c_df-c)), 0, 0.004));
      gamma <- matrix(runif(((c_df-c)*(c_df-c)),first_gamma, last_gamma));
      dim(gamma) <- c((c_df-c),(c_df-c))
      
      #initialization
      meanS <- integer(iter);
      minS <- integer(iter);
      maxS <- integer(iter);
      for(i in 1:iter){
        #if(i>1){
        #gamma <- (t(alpha-alpha_pre) %*% (grad-grad_pre))/sum((grad-grad_pre)^2);
        alpha_pre <- alpha;
        
        #}
        
        grad <- t((1/c)*(t((AA%*%alpha)-BB)%*%AA));
        grad_pre<-grad;
        alpha <- alpha_pre - ((grad-betha) %*%gamma)
        
        err_element <- XX-CC%*%alpha;
        err <- sum(abs(err_element));
        meanErr <- sqrt(err / ((r_df-r)*(c_df-c)));
        if(meanErr>100 ){
          break
        }
        meanS[i]<-meanErr;
        minS[i]<-min(err_element);
        maxS[i]<-max(err_element);
        meanErr
        
      }
      all_meanS <- cbind(all_meanS,meanS)
    }
    all_percents<- rbind(all_percents,all_meanS)
  }
  return(all_percents)
  
}

nearMatrix <- function(expr3795,cs,rs,all_percent_cols){
  
  df<-data.frame(t(expr3795),check.names = FALSE);
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  
  #initialization
  meanS <- NULL;
  minS <- NULL;
  maxS <- NULL;
  
  r_df<-nrow(df);
  c_df<-ncol(df);
  for(i in 1:9){
    c <- cs[i]
    r <- rs[i]
    cols <- all_percent_cols[i,];
    
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    Bd <-NULL;
    Ad <-NULL;
    
    print("loop")
    print(i)
    for(element in 1:(r_df-r)){
      print(element)
      ref<-rep(CC[element,],(r_df-r))
      dim(ref) <- c((r_df-r),(c))
      distance<-apply(AA,1,function(x)sqrt(sum((x-ref)^2)))
      index<-match(min(distance),distance)
      Ad <- rbind(Ad,AA[index[1],]);
      Bd <- rbind(Bd,BB[index[1],]);
      print(" another loop")
    }
    
    #print(abaad(Ad))
    #print(abaad(Bd))
    #print(abaad(CC))
    estB <- (CC %*% pseudoinverse(Ad)) %*% Bd;
    main <- rbind(cbind(AA,BB),cbind(CC,estB));
    errB <-sum(abs(estB-XX)^2);
    errB
    meanErrB <- sqrt(errB / (nrow(XX)*ncol(XX)));
    meanS <- rbind(meanS,meanErrB);
    minS <- rbind(minS,min(estB-XX));
    maxS <- rbind(maxS,max(estB-XX));
  }
  return(meanS)
}

clusterNearMatrix <- function(expr3795,clusters,cs,rs,all_percent_cols){
  
  df<-data.frame(t(expr3795),check.names = FALSE);
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  #dfc <- scale(df[-1]) 
  dft<-data.frame((expr3795),check.names = FALSE);
  dft <- data.matrix(dft, rownames.force = NA )
  dft <- dft[,colSums(is.na(dft))<nrow(dft)]
  
  res <- kmeans(dft, centers=clusters, nstart=10)
  #res1 <- kcca(dft,3,dist=cov())
  #res <- kcca2df(res1)
  
  #initialization
  meanSP <- NULL;
  minSP <- NULL;
  maxSP <- NULL;
  
  
  r_df<-nrow(df);
  c_df<-ncol(df);
  for(i in 1:9){
    print(i)
    c <- cs[i]
    r <- rs[i]
    cols <- all_percent_cols[i,];
    
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    
    
    minP<-NULL;
    maxP<-NULL;
    errP<-NULL;
    
    for(clus in 1:clusters){
      truecol1<-names(which(res$cluster==clus));
      idx_first <- match(truecol1,colnames(CC));
      idx_second<-match(truecol1,colnames(BB));
      f<- idx_first[is.na(idx_first)==FALSE];
      s<- idx_second[is.na(idx_second)==FALSE];
      if(length(f)<2){
        CP<-CC;
        AP<-AA;
        BP<-BB[,s];
        XP1<-XX[,s];
        
        Bd <-NULL;
        Ad <-NULL;
        
        for(element in 1:(r_df-r)){
          ref<-rep(CP[element,],(r_df-r))
          dim(ref) <- c((r_df-r),(c))
          distance<-apply(AP,1,function(x)sqrt(sum((x-ref)^2)))
          index<-match(min(distance),distance)
          Ad <- rbind(Ad,AP[index[1],]);
          Bd <- rbind(Bd,BP[index[1],]);
        }
        AP<-Ad;
        BP<-Bd;
        
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }else if(length(s)==0){
        meanErrP1<-0;
      }else{
        CP<-CC[,f];
        AP<-AA[,f];
        BP<-BB[,s];
        XP1<-XX[,s];
        
        Bd <-NULL;
        Ad <-NULL;
        
        for(element in 1:(r_df-r)){
          ref<-rep(CP[element,],(r_df-r))
          dim(ref) <- c((r_df-r),(ncol(CP)))
          distance<-apply(AP,1,function(x)sqrt(sum((x-ref)^2)))
          index<-match(min(distance),distance)
          Ad <- rbind(Ad,AP[index[1],]);
          Bd <- rbind(Bd,BP[index[1],]);
        }
        AP<-Ad;
        BP<-Bd;
        
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }
    }
    
    meanErrP <- sqrt(sum(errP) / (nrow(XX)*ncol(XX)));
    meanSP <- rbind(meanSP,meanErrP);
    minSP <- rbind(minSP,min(minP));
    maxSP <- rbind(maxSP,max(maxP));
  }
  return(meanSP)
}

nearMatrixClustering <- function(expr3795,clusters,cs,rs,all_percent_cols){
  
  df<-data.frame(t(expr3795),check.names = FALSE);
  df <- data.matrix(df, rownames.force = NA )
  df <- df[,colSums(is.na(df))<nrow(df)]
  #dfc <- scale(df[-1]) 
  dft<-data.frame((expr3795),check.names = FALSE);
  dft <- data.matrix(dft, rownames.force = NA )
  dft <- dft[,colSums(is.na(dft))<nrow(dft)]
  
  res <- kmeans(dft, centers=clusters, nstart=10)
  #res1 <- kcca(dft,3,dist=cov())
  #res <- kcca2df(res1)
  
  #initialization
  meanSP <- NULL;
  minSP <- NULL;
  maxSP <- NULL;
  
  
  r_df<-nrow(df);
  c_df<-ncol(df);
  for(i in 1:9){
    print(i)
    c <- cs[i]
    r <- rs[i]
    cols <- all_percent_cols[i,];
    
    CC <- df[c((r+1):r_df),cols];
    XX <- df[c((r+1):r_df),-cols];
    AA <- df[c(1:r),cols];
    BB <- df[c(1:r),-cols];
    
    Bd <-NULL;
    Ad <-NULL;
    
    for(element in 1:(r_df-r)){
      ref<-rep(CC[element,],(r_df-r))
      dim(ref) <- c((r_df-r),c)
      distance<-apply(AA,1,function(x)sqrt(sum((x-ref)^2)))
      index<-match(min(distance),distance)
      Ad <- rbind(Ad,AA[index[1],]);
      Bd <- rbind(Bd,BB[index[1],]);
    }
    AA<-Ad;
    BB<-Bd;    
    
    minP<-NULL;
    maxP<-NULL;
    errP<-NULL;
    
    for(clus in 1:clusters){
      truecol1<-names(which(res$cluster==clus));
      idx_first <- match(truecol1,colnames(CC));
      idx_second<-match(truecol1,colnames(BB));
      f<- idx_first[is.na(idx_first)==FALSE];
      s<- idx_second[is.na(idx_second)==FALSE];
      if(length(f)<2){
        CP<-CC;
        AP<-AA;
        BP<-BB[,s];
        XP1<-XX[,s];
        
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }else if(length(s)==0){
        meanErrP1<-0;
      }else{
        CP<-CC[,f];
        AP<-AA[,f];
        BP<-BB[,s];
        XP1<-XX[,s];
        
        
        estP1 <- (CP %*% pseudoinverse(AP)) %*% BP;
        maxP <- cbind(maxP,max(estP1-XP1));
        minP <- cbind(minP,min(estP1-XP1));
        errP1 <-sum(abs(estP1-XP1)^2);
        errP <- cbind(errP,errP1);
        meanErrP1 <- sqrt(errP1 / (nrow(XP1)*ncol(XP1)));
      }
    }
    
    meanErrP <- sqrt(sum(errP) / (nrow(XX)*ncol(XX)));
    meanSP <- rbind(meanSP,meanErrP);
    minSP <- rbind(minSP,min(minP));
    maxSP <- rbind(maxSP,max(maxP));
  }
  return(meanSP)
}


expr<-expr11;
#predictions
all_predictions <- function(expr){
  
  #random missing
  all_percent_cols<-NULL;
  r<-NULL;
  c<-NULL;
  r_df <-ncol(expr); 
  c_df <-nrow(expr);
  k=10;
  for(percents in 1:9){
    percent <- percents*10;# darsad e missing 
    temp <- integer(c_df);
    c <- cbind(c,floor(c_df-sqrt((1-percent/100)*c_df)));
    r <- cbind(r,floor(r_df-sqrt(1-percent/100)*r_df));
    temp[1:c[percents]] <- sample.int(c_df,c[percents]);
    all_percent_cols<-rbind(all_percent_cols,temp);
  }
  
  mean6<-prediction(data.frame(t(expr)),c,r,all_percent_cols);
  knnimpute <- KNN_impute(data.frame(t(expr)),c,r,all_percent_cols,k); 
  #iter_mean6<-iterative(data.frame(t(expr)),c,r,all_percent_cols);
  cluster_mean6_3<-clustering(expr,3,c,r,all_percent_cols);
  cluster_mean6_2<-clustering(expr,2,c,r,all_percent_cols);
  cluster_mean6_4<-clustering(expr,4,c,r,all_percent_cols);
  cluster_mean6_5<-clustering(expr,5,c,r,all_percent_cols);
  #regres_means<- regres(expr,c,r,all_percent_cols,3,10,500,0.0001,0.0000001, 0.000001,0.0000001,0.000001);
  nearMatrix_means <- nearMatrix(expr,c,r,all_percent_cols);
  clusterNearMatrix_means6_5<-clusterNearMatrix(expr,2,c,r,all_percent_cols);
  clusterNearMatrix_means6_2<-clusterNearMatrix(expr,3,c,r,all_percent_cols);
  clusterNearMatrix_means6_3<-clusterNearMatrix(expr,4,c,r,all_percent_cols);
  clusterNearMatrix_means6_4<-clusterNearMatrix(expr,5,c,r,all_percent_cols);
  nearMatrixClustering_means6_2<-nearMatrixClustering(expr,2,c,r,all_percent_cols);
  nearMatrixClustering_means6_3<-nearMatrixClustering(expr,3,c,r,all_percent_cols);
  nearMatrixClustering_means6_4<-nearMatrixClustering(expr,4,c,r,all_percent_cols);
  nearMatrixClustering_means6_5<-nearMatrixClustering(expr,5,c,r,all_percent_cols);
  
  #ploting
  observability <- rbind(10,20,30,40,50,60,70,80,90)
  data <- cbind(observability,mean6,nearMatrix_means,cluster_mean6_2,cluster_mean6_3,cluster_mean6_4,cluster_mean6_5,clusterNearMatrix_means6_2,clusterNearMatrix_means6_3,clusterNearMatrix_means6_4,clusterNearMatrix_means6_5,nearMatrixClustering_means6_2,nearMatrixClustering_means6_3,nearMatrixClustering_means6_4,nearMatrixClustering_means6_5)
  colnames(data) <- c("observability","simple","nearMatrix","two_cluster","tree_cluster","four_cluster","five_cluster","two_clusterAndNearMatrix","tree_clusterAndNearMatrix","four_clusterAndNearMatrix","five_clusterAndNearMatrix","two_nearMatrixClustering","tree_nearMatrixClustering","four_nearMatrixClustering","five_nearMatrixClustering")
  data <- data.frame(data);
  return(data)
}


#Data
eset6 <- GDS2eSet(getGEO("GDS4971"),do.log2=TRUE)
expr6 <- exprs(eset6)
expr6<-expr6[complete.cases(expr6), ]

eset7 <- GDS2eSet(getGEO("GDS4296"),do.log2=TRUE)
expr7 <- exprs(eset7)

eset8 <- GDS2eSet(getGEO("GDS3312"),do.log2=TRUE)
expr8 <- exprs(eset8)

eset10 <- GDS2eSet(getGEO("GSE55347"),do.log2=TRUE)
expr10 <- exprs(eset8)

eset11 <- GDS2eSet(getGEO("GSM2284083"),do.log2=TRUE)
expr11 <- exprs(eset8)

load("Desktop/oneToOne.RData")
missingDatas<-which(is.na(RNAOneToOneMatrix),arr.ind = T);
uMissingDatas <- unique(missingDatas[,2]);
micro <- oneToOneRMAMatrix[,-uMissingDatas];
rna <- RNAOneToOneMatrix[,-uMissingDatas];

expr9<-t(rna)
expr<- expr9;
mean9<-prediction(data.frame(t(expr)),c,r,all_percent_cols);
iter_mean9<-iterative(data.frame(t(expr)),c,r,all_percent_cols);

data6<-all_predictions(expr6)
data7<-all_predictions(expr7)
data8<-all_predictions(expr8)
data9<-all_predictions(expr9)

percent <- 10;# darsad e missing 
df <- t(df)
r_df <- abaad(df)[1,1]
c_df <- abaad(df)[2,1]
c <- (c_df-sqrt((1-percent/100)*c_df));
r <- (r_df-sqrt(1-percent/100)*r_df);
c<-ceiling(c)-1
r<-ceiling(r)-1
C <- df[c((r+1):r_df),c(1:c)];
X <- df[c((r+1):r_df),((c+1):c_df)];
A <- df[c(1:r),c(1:c)];
B <- df[c(1:r),c((c+1):c_df)];

est <- (C %*% pseudoinverse(A)) %*% B[,1];
err <-sum(abs(est-X[,1])^2);
err
meanErr <- sqrt(err / (nrow(X)*ncol(X)))

require(glmnet)
set.seed(999)
cv.lasso <- cv.glmnet(A, B[,1], alpha=1, parallel=TRUE, standardize=TRUE, type.measure='mse')
lasso.pred <- predict(cv.lasso, s = cv.lasso$lambda.min, newx = C)
mean((lasso.pred-X[,1])^2)
cv.lasso$glmnet.fit

NAp <- rep(NA, ncol(X)*nrow(X)); 
dim(NAp) <- c(nrow(X),ncol(X));
abaad(NAp)
main <- rbind(cbind(A,B),cbind(C,NAp));
knn <- cv.kNNImpute(main, k, x.dist = NULL, impute.fn, verbose = T)

#ploting

#simple
means <- (cbind(mean6,mean7,mean8,mean9,observability));
rownames(means) <- c("10%","20%","30%","40%","50%","60%","70%","80%","90%");
colnames(means) <- c("GDS4971","GDS4296","GDS3312","sample","observability");
meanS<-data.frame(means)
p1<-ggplot()+
  geom_line(data=meanS,aes(x=observability,y=GDS4971,colour="GDS4971"))+
  geom_line(data=meanS,aes(x=observability,y=GDS4296,colour="GDS4296"))+
  geom_line(data=meanS,aes(x=observability,y=GDS3312,colour="GDS3312"))+
  geom_line(data=meanS,aes(x=observability,y=sample,colour="sample"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 0.9),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p1




#simple-iter
best <- apply(iter_mean9, 1, min) 
data <- cbind(observability,mean9,iter_mean9[,ncol(iter_mean9)],best)
colnames(data) <- c("observability","simple","iterative","best")
data <- data.frame(data)
p2<-ggplot()+
  geom_line(data=data,aes(x=observability,y=best,colour="iterative"))+
  geom_line(data=data,aes(x=observability,y=simple,colour="simple"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p2

#iters
data_iter <- t(iter_mean9)
colnames(data_iter) <- c("one","two","three","four","five","six","seven","eight","nine")
data_iter <- data.frame(data_iter);
num <- 1:2001 
data_iter <- cbind(num,data_iter)
p3<-ggplot()+
  geom_line(data=data_iter,aes(x=num,y=one,colour="10"))+
  geom_line(data=data_iter,aes(x=num,y=two,colour="20"))+
  geom_line(data=data_iter,aes(x=num,y=three,colour="30"))+
  geom_line(data=data_iter,aes(x=num,y=four,colour="40"))+
  geom_line(data=data_iter,aes(x=num,y=five,colour="50"))+
  geom_line(data=data_iter,aes(x=num,y=six,colour="60"))+
  geom_line(data=data_iter,aes(x=num,y=seven,colour="70"))+
  geom_line(data=data_iter,aes(x=num,y=eight,colour="80"))+
  geom_line(data=data_iter,aes(x=num,y=nine,colour="90"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p3


#simple-cluster
data_cluster <- cbind(observability,mean6,cluster_mean6_2,cluster_mean6_3,cluster_mean6_4,cluster_mean6_5)
colnames(data_cluster) <- c("observability","simple","two_cluster","tree_cluster","four_cluster","five_cluster")
data_cluster <- data.frame(data_cluster);
p4<-ggplot()+
  geom_line(data=data_cluster,aes(x=observability,y=simple,colour="simple"))+
  geom_line(data=data_cluster,aes(x=observability,y=two_cluster,colour="two_cluster"))+
  geom_line(data=data_cluster,aes(x=observability,y=tree_cluster,colour="tree_cluster"))+
  geom_line(data=data_cluster,aes(x=observability,y=four_cluster,colour="four_cluster"))+
  geom_line(data=data_cluster,aes(x=observability,y=five_cluster,colour="five_cluster"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p4


#simple-nearMatrix
means <- cbind(observability,nearMatrix_means,(mean6[2:length(mean6)]))
colnames(means) <- c("observability","nearMatrix","simple");
meanS<-data.frame(means)
p5<-ggplot()+
  geom_line(data=meanS,aes(x=observability,y=nearMatrix,colour="nearMatrix"))+
  geom_line(data=meanS,aes(x=observability,y=simple,colour="simple"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p5

#simple-nearMatrix-cluster
data_cluster <- cbind(observability,mean6,nearMatrix_means,cluster_mean6_2,cluster_mean6_3,cluster_mean6_4,cluster_mean6_5)
colnames(data_cluster) <- c("observability","simple","nearMatrix","two_cluster","tree_cluster","four_cluster","five_cluster")
data_cluster <- data.frame(data_cluster);
p6<-ggplot()+
  geom_line(data=meanS,aes(x=observability,y=nearMatrix,colour="nearMatrix"))+
  geom_line(data=data_cluster,aes(x=observability,y=simple,colour="simple"))+
  geom_line(data=data_cluster,aes(x=observability,y=five_cluster,colour="five_cluster"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p6


#both simple-nearMatrix-cluster
data <- cbind(observability,mean6,nearMatrix_means,cluster_mean6_2,cluster_mean6_3,cluster_mean6_4,cluster_mean6_5,clusterNearMatrix_means6_2,clusterNearMatrix_means6_3,clusterNearMatrix_means6_4,clusterNearMatrix_means6_5,nearMatrixClustering_means6_4)
colnames(data) <- c("observability","simple","nearMatrix","two_cluster","tree_cluster","four_cluster","five_cluster","two_clusterAndNearMatrix","tree_clusterAndNearMatrix","four_clusterAndNearMatrix","five_clusterAndNearMatrix","nearMatrixClustering_means6_4")
data <- data.frame(data);
p7<-ggplot()+
  geom_line(data=data,aes(x=observability,y=nearMatrix,colour="nearMatrix"))+
  geom_line(data=data,aes(x=observability,y=simple,colour="simple"))+
  geom_line(data=data,aes(x=observability,y=five_cluster,colour="five_cluster"))+
  geom_line(data=data,aes(x=observability,y=tree_clusterAndNearMatrix,colour="tree_clusterAndNearMatrix"))+
  geom_line(data=data,aes(x=observability,y=four_clusterAndNearMatrix,colour="four_clusterAndNearMatrix"))+
  geom_line(data=data,aes(x=observability,y=five_clusterAndNearMatrix,colour="five_clusterAndNearMatrix"))+
  geom_line(data=data,aes(x=observability,y=nearMatrixClustering_means6_4,colour="nearMatrixClustering_means6_4"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous(expand = c(0, 0))+
  ylab("mean error")+
  scale_color_discrete(name ="datasets")
p7
