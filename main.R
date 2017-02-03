################################################
#Non-Null distribution is estimated by B-splines
#CM-Alpha: 
################################################

getwd()
setwd("/home/yourFolder") #put code in certain directory

rm(list=ls())
 
load("toyData.RData") #Only information needed from the data is X and Z (others can be ignored)

source("LFDR_CMBsplines.R")
X=data.matrix(X);
N=length(Z)

#For test run, take small iterations, such as nIter=1100, burnIn=100, thin=5, etc.
thin=10;
nIter=18000;
burnIn=4000;
useCov=1;


if(useCov==1){
 	M=dim(X)[2]; #M: number of the covariates including intercept
 	SSA=c(20,rep(10,M-1)) #change step-size for alpha draw; length(SSA)=M 
 	SSAP=rep(0.5,M); # prob that step-size of alpha will be changed to SSA; length(SSAP)=M;
 	MA=rep(5,M); #num of multiple try for alpha; length(MA)=M;
 }else{
 	M=1;
 	SSA=20;SSAP=0.4;MA=5;
 }
 
SSG=2 # change step-size for gamma draw
MG=3 #num of multiple try for gamma

K=5;
mu=0.68;
initNULL=0.95;


#process begins;
ptm <- proc.time() 

MCMCfit<-cmlFDR_CMBsplinesDist(Z,X,K,nIter,burnIn,thin,initNULL,SSA=SSA,SSAP=SSAP,SSG=SSG,MA=MA,MG=MG,mu=mu,useCov=useCov)

#process ends;
ptm2<-proc.time()
#process time needed in hour
(ptm2 - ptm)/3600 



#save file;
save(file="result.R",MCMCfit,X,Z,N,nIter,burnIn,thin,mu,K,useCov)
	




if(useCov==1){
	M=dim(X)[2]
}else{
	M=1;
	X=matrix(1,nrow=N,ncol=1);
}
#Alpha
AMN<-matrix(0,nrow=K,ncol=M);
print("Mean of Alpha");
for (m in 1:M){		
	for(k in 1:K){   		 
		AMN[k,m]=mean(sapply(MCMCfit[[1]],'[',m,k))
   	}
}
print(AMN)


ASD<-matrix(0,nrow=K,ncol=M);
print("Alpha SD:");
for (m in 1:M){
	for(k in 1:K){
		ASD[k,m]=sd(sapply(MCMCfit[[1]],'[',m,k))
			
	}
}
print(ASD)

AMD<-matrix(0,nrow=K,ncol=M);
print("Median of Alpha");
for (m in 1:M){		
	for(k in 1:K){
   		 
		AMD[k,m]=median(sapply(MCMCfit[[1]],'[',m,k))
   	}
}
AMD


#percnetile of alpha
L<-matrix(0,nrow=K,ncol=M);
U<-matrix(0,nrow=K,ncol=M);
for (m in 1:M){
		for(k in 1:K){
  			L[k,m]=quantile(sapply(MCMCfit[[1]],'[',m,k),prob=0.025)
 			U[k,m]=quantile(sapply(MCMCfit[[1]],'[',m,k),prob=0.975)
		}
}
print("posterior quantile 0.025:");
print(L);
print("posterior quantile 0.975:");
print(U);


#Tau_sq
print("Tau_sq mean:");
TMN<-NULL
for(i in 1:M){
    	TMN<-rbind(TMN,mean(as.numeric(unlist(lapply(MCMCfit[[3]], function(x) x[i])))))
	
	}
TMN

print("Tau_sq SD:");
for(i in 1:M){
    	print(sd(as.numeric(unlist(lapply(MCMCfit[[3]], function(x) x[i])))))
	}

print("Tau_sq median:");
TMD<-NULL
for(i in 1:M){
    	TMD<-rbind(TMD,median(as.numeric(unlist(lapply(MCMCfit[[3]], function(x) x[i])))))

}
TMD

#percnetile Tau-sq
print("Tau_sq 95% posterior credible region:");
for(i in 1:M){
  
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[3]], function(x) x[i]))),prob=0.025))
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[3]], function(x) x[i]))),prob=0.975))
}

#####
#a
#####
print("a mean:");
aMN<-NULL
for(i in 1:M){
    	aMN<-rbind(aMN,mean(as.numeric(unlist(lapply(MCMCfit[[15]], function(x) x[i])))))
	
	}
aMN

print("a SD:");
for(i in 1:M){
    	print(sd(as.numeric(unlist(lapply(MCMCfit[[15]], function(x) x[i])))))
	}

print("a median:");
aMD<-NULL
for(i in 1:M){
    	aMD<-rbind(aMD,median(as.numeric(unlist(lapply(MCMCfit[[15]], function(x) x[i])))))

}
aMD

#percnetile a
print("a 95% posterior credible region:");
for(i in 1:M){
  
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[15]], function(x) x[i]))),prob=0.025))
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[15]], function(x) x[i]))),prob=0.975))
}

#sigma_sq
print(paste("sigma_sq mean:",mean(as.numeric(MCMCfit[[13]]))))
print(paste("sigma_sq SD:",sd(as.numeric(MCMCfit[[13]]))))
print(paste("sigma_sq median:",median(as.numeric(MCMCfit[[13]]))))

print("sigma_sq 95% posterior credible region:");
quantile(as.numeric(MCMCfit[[13]]),prob=0.025)
quantile(as.numeric(MCMCfit[[13]]),prob=0.975)


#Gamma
print("Gamma mean:");
GMN<-NULL
for(i in 1:M){
    	GMN<-rbind(GMN,mean(as.numeric(unlist(lapply(MCMCfit[[2]], function(x) x[i])))))
	
	}
GMN

print("Gamma SD:");
for(i in 1:M){
    	print(sd(as.numeric(unlist(lapply(MCMCfit[[2]], function(x) x[i])))))
	}

print("Gamma median:");
GMD<-NULL
for(i in 1:M){
    	GMD<-rbind(GMD,median(as.numeric(unlist(lapply(MCMCfit[[2]], function(x) x[i])))))

}
GMD

#percnetile gamma
print("Gamma 95% posterior credible region:");
for(i in 1:M){
  
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[2]], function(x) x[i]))),prob=0.025))
  print(quantile(as.numeric(unlist(lapply(MCMCfit[[2]], function(x) x[i]))),prob=0.975))
}

X[1:5,]

#total number of runs after thinning;
rn=length((MCMCfit)[[1]])


##########################################
#CM-LFDR
##########################################
mu=0.68
f0<-2*dnorm(abs(Z),mean=0,sd=sqrt(median(as.numeric(MCMCfit[[13]]))))
f1<-f0 #just declare a vector with the same lenth as f0
f1=rowSums(exp((X%*%t(AMD) + log(MCMCfit[[12]])))/rowSums(exp(X%*%t(AMD))))
f1[abs(Z)<=mu]=0

p0=1;
p1=exp(X%*%GMD)


pi0=p0/(p0+p1)
pi1=1-pi0



CM.lfdr=pi0*f0/(pi0*f0+pi1*f1)


alpha=seq(0,0.2,0.005)
totalNonNull<-NULL
NonNULL<-matrix(0,nrow=N,ncol=length(alpha));
j=1;
#cutoff
for (i in alpha){
	NonNULL[,j]=ifelse(CM.lfdr<=i,1,0)
	totalNonNull[j]=sum(NonNULL[,j])
	j=j+1;
}
cbind(alpha,totalNonNull);
