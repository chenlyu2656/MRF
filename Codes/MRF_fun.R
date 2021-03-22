#Y: trait,n length vector; Z: genotype matrix n*m ; X: covariate matrix, n*p
MRF<-function(Y,Z,X=NULL,kernel='GR',centerS=F,weights=1,type='normal',od=3) 
{
  #Preliminary calculation
  n<-length(Y);m<-ncol(Z);X1<-cbind(rep(1,n),X)
  #Define genetic similarity
  if(kernel=='GR'){S<-getGR(Z,weights)}
  if(kernel=='IBS'){S<-getIBS(Z,weights)}
  if(kernel=='EDS'){S<-getEDS(Z,weights)}
  if(kernel=='DPS'){S<-getDPS(Z,weights,od=od)}
  if(kernel=='LIN'){S<-getLIN(Z,weights)}
  #Non-centered or Centered similarity
  #Non-centered: enhanced power for one direction effect, high causal proportion
  #Centered: more robust performance to bi-direction effect for rare variants
  if(centerS==T){J<-diag(rep(1,n))-matrix(1/n,n,n);S<-J%*%S%*%J;}	
  diag(S)<-0
  #Quantitative trait
  if(type=='normal')
  {
    #calculate residual
    B=diag(1,n)-X1%*%solve(t(X1)%*%X1)%*%t(X1)
    Y.res<-B%*%matrix(Y,ncol=1);
    E<-S%*%Y.res
    #observed test statistic
    r.GRF<-as.numeric(solve(t(E)%*%E)%*%t(E)%*%Y.res)
    #null distribution
    A<-S-S%*%S*r.GRF
    M1<-B%*%A%*%B
  }
  #Binary trait
  if(type=='binary')
  {
    #calculate residual
    fit<-glm(Y~X1,family=binomial(link="logit"))
    mu<-fit$fitted
    Y.res<-as.matrix(Y-mu)
    E<-S%*%Y.res;
    #observed test statistic
    r.GenRF<-as.numeric(solve(t(E)%*%E)%*%t(E)%*%Y.res)
    #null distribution
    A<-S-S%*%S*r.GenRF;
    V<-diag(mu*(1-mu));
    P0<-V-V%*%X1%*%solve(t(X1)%*%V%*%X1)%*%t(X1)%*%V
    K<-Get.sqrt(P0)
    M1<-K%*%A%*%K
  }
  # beta distributed trait
  if(type=='beta')
  {
    #calculate residual
    fit<-betareg(Y~X,link = "logit")
    mu<-fit$fitted
    phi<-summary(fit)$coef$precision[,1]
    Y.res<-as.matrix(Y-mu)
    E<-S%*%Y.res;
    #observed test statistic
    r.GenRF<-as.numeric(solve(t(E)%*%E)%*%t(E)%*%Y.res)
    #null distribution
    A<-S-S%*%S*r.GenRF;
    V<-diag(mu*(1-mu)/(phi+1));
    P0<-V-V%*%X1%*%solve(t(X1)%*%V%*%X1)%*%t(X1)%*%V
    K<-Get.sqrt(P0)
    M1<-K%*%A%*%K
  }
  #Cauculate p-value
  p.value<-davies(0,eigen(M1,symmetric=TRUE,only.value=TRUE)$values,acc = 10^(-8))$Qq
  return(list(pvalue=p.value))
}



###Calculate genetic distance : IBS

getIBS<-function(geno,weights=1)
{
  n<-nrow(geno);m<-ncol(geno);
  gtemp1<-geno;gtemp2<-geno;
  gtemp1[geno==2]<-1;gtemp2[geno==1]<-0;gtemp2[geno==2]<-1;
  gtemp<-cbind(gtemp1,gtemp2);
  Inner<-gtemp%*%diag(weights,nrow=2*m,ncol=2*m)%*%t(gtemp);
  X2<-matrix(diag(Inner),nrow=n,ncol=n);
  Y2<-matrix(diag(Inner),nrow=n,ncol=n,byrow=T);
  Bound<-sum(matrix(weights,nrow=1,ncol=m)*2);
  Dis<-(X2+Y2-2*Inner);
  IBS<-Bound-Dis;
  IBS;
}


getGR<-function(geno,weights=1){
  n<-nrow(geno);
  m<-ncol(geno);
  cm<-colMeans(geno);	
  gtemp<-geno-matrix(cm,nrow=n,ncol=m,byrow=T);
  W<-diag(weights,m,m);
  GR<-gtemp%*%W%*%t(gtemp)/m;
  GR;
}


getLIN<-function(geno,weights=1)
{
  n<-nrow(geno);m<-ncol(geno);maf<-colMeans(geno)/2;	
  gtemp<-geno-matrix(maf*2,nrow=n,ncol=m,byrow=T);
  W<-diag(weights,m,m);
  LIN<-gtemp%*%W%*%t(gtemp)/m;
  LIN;
}

getEDS<-function(geno,weights=1)
{
  n<-nrow(geno);m<-ncol(geno);
  gtemp<-geno;gtemp[is.na(geno)]<-0;
  Inner<-(gtemp)%*%diag(weights,nrow=m,ncol=m)%*%t(gtemp);
  X2<-(gtemp)^2%*%diag(weights,nrow=m,ncol=m)%*%(is.na(t(geno))==0); 
  Y2<-t(X2);
  Wtemp<-matrix(weights,nrow=n,ncol=m,byrow=T);
  Wtemp[is.na(geno)]<-0;
  Bound<-(is.na(geno)==0)%*%diag(weights,nrow=m,ncol=m)%*%(is.na(t(geno))==0)*4;
  ED<-X2+Y2-2*Inner;
  ED[ED<0]<-0;
  EDS<-sqrt(Bound)-sqrt(ED);
  EDS;
}

getDPS<-function(geno,weights=1,od=3)
{
  n<-nrow(geno);m<-ncol(geno);
  DPS<-matrix(0,nrow=n,ncol=n);
  for(i in 1:n)
  {
    gtemp<-matrix(geno[i,],nrow=n,ncol=m,byrow=T);
    temp<-(abs(gtemp-geno))^od;
    score<-rowSums(temp*matrix(weights,nrow=n,ncol=m,byrow=T),na.rm=T);
    D<-score^(1/od);
    Bound<-rowSums((matrix(weights,nrow=n,ncol=m,byrow=T))*(1-is.na(geno))*(1-is.na(gtemp)))^(1/od)*2;
    DPS[,i]<-Bound-D;
  }
  DPS;
}

Get.sqrt<-function(A)
{
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
