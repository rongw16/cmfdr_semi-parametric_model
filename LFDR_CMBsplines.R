
cmlFDR_CMBsplinesDist=function (Z,X,K=20,nIter=1100,burnIn=100,thin=5,initNULL=0.95,simulate=FALSE,SSA,SSAP,SSG=1,MA,MG=3,theoNULL=2,mu=0.68,return_Bden=1,Liby=1,useCov=1,const=100,nu=10,A=10,initYN=2,saveWhileRun=2,initAlpha=NULL,initGamma=NULL,initTauSq=NULL,initSigmaSq=NULL,init.a=NULL) 
{
#SSA:scale the diagnal of the covariance matrix in MH of Alpha draw to increase/decrease the step size
#SSG:scale the diagnal of the covariance matrix in MH of Gamma draw to increase/decrease the step size

#initAlpha: M(row) by K (col) matrix where M is the number of covariates, includign intercept, K is the number of the B-spline component; 1st col all 0;
#initGamma: M by 1 vector
#initTauSq: scalar
#initSigmaSq: scalar
#saveWhileRun: 1: save 4 times during run; 2: don't save during the run; defulat is 2, not saving.

###########################
#Load packages/functions
###########################
	if(Liby==1){
		library(tmvtnorm)
		library(mnormt)
		library(magic)
		library(arm)
		library(pscl)
		library(splines)
		library(plyr)

		source("fcns.R")
	}


	if(K<=4){stop("K needs to be >=5.");return; }
	
	N=length(Z)


	if (useCov==1){
		if(length(which(X[,1]==1))!= length(Z)) {X=cbind(Intcpt=1,X)}
		
		M=dim(X)[2]

	}else{
		X=matrix(1,nrow=N,ncol=1)

		M=1;
	}
	
	# hyperparameters	
	a0=0.001;b0=0.001; #Prior: P(sigma^2) ~ IVG(a0,b0) #b0 is scale;
	Sigma_Gamma=matrix(0,nrow=M,ncol=M)
	diag(Sigma_Gamma)=10000; #Prior: P(Gamma) ~ N(0,Sigma_Gamma)

	
	D<-matrix(0,nrow=K-3,ncol=K-1)
	for (i in 1:(K-3)){
		D[i,i:(i+2)]=c(1,-2,1);
	}
	Omega=t(D)%*%D;
	Omega_star=Omega;
	Omega_star[1,1]=Omega[1,1]+1/const;
	Omega_star[2,2]=Omega[2,2]+1/const;
	
	#P(tau_m^2|a_m) ~ IVG(nu/2,nu/a_m);
	#P(a_m) ~ IVG(1/2,1/A^2)
	#nu and A are from input or default
	
	

	df=4;
	

	# data

	if(length(SSA) != M | length(SSAP) != M | length(MA)!= M){
		stop("length(SSA),length(SSAP) and length(MA) should be the number of the covariates (including intercept)."); 
		return;
	}


	# Parameter arrays	
	
	ALPHA_array=list()
	Tau_sq_array=list()
	a_array=list()
	if(theoNULL !=1){
		SIGMA_SQ_array=list()
	}
	GAMMA_array=list()
	Accp_Rate_array_a=list() 	#Alpha draw accept rate
	Accp_Rate_array_g=list()	#Gamma draw accept rate
	
	array_ind=0

	# Initialize global parameters 	
	pi0=initNULL 
	
	if (initYN!=1){
		gamma0=log((1-pi0)/pi0) 
		Gamma=array(0,dim=c(M,1)); Gamma[1]=gamma0
		Sigma_sq=1;
		Tau_sq=array(1,dim=c(M,1)) ; 
		a=array(1,dim=c(M,1));
	}else{
		Gamma=initGamma;
		Tau_sq=initTauSq;
		a=init.a;
		Sigma_sq=initSigmaSq;
		
	}
	Gamma_mean=array(0,dim=c(M,1));
	Tau_sq_mean=array(0,dim=c(M,1));
	a_mean=array(0,dim=c(M,1));
	
	Phi=1-as.numeric(abs(Z)< sort(abs(Z))[round(pi0*N)]);
	Phi_mean=0*Phi	
	PHI_match_rate=NULL

	#Initialze Alpha;
	#Constrain on C: Sum(C)=1, 0<C[i]<1, exp(Alpha[,1])=1;
	if (initYN!=1){
	
		Alpha=matrix(0,nrow=M,ncol=K);

		for (m in 1:M){
			Alpha[m,2:3]=mvrnorm(n = 1, c(0,0), matrix(c(Tau_sq[m],0,0,Tau_sq[m]),nrow=2,ncol=2,byrow=T))
			for (i in 4:K){
				Alpha[m,i]=rnorm(1,mean=2*Alpha[m,i-1]-Alpha[m,i-2],sd=sqrt(Tau_sq[m]))	
			}
		}
	}else{
		Alpha=initAlpha;
		
	}
	Alpha_mean=matrix(0,nrow=M,ncol=K);
	Alpha_accp=rep(0,M); #draw alpha[m,] from m=1 to M;

	cnt=0;
	
	Gden<-matrix(0,nrow=N,ncol=K);  
	byStep=0.001;
	
	
	source("otherFcns.R")

	Gden<-matrix(0,nrow=N,ncol=K);
	Gden<-Bdensity(byStep,mu,K,N,Z,Gden)



	#beginning of the iterations
	for(iter in 1:nIter){
		
		print(iter)
		
		Z1=abs(Z[Phi==1])
			
		if(useCov==1){
			X1=X[Phi==1,]
		}else{
			X1=X[Phi==1]
		}

		nn=(1:N)[Phi==1]
		
		###############################
		#get B-spline densities for Z1;
		###############################
		gden<-matrix(0,nrow=length(Z1),ncol=K);
		gden<-Gden[nn,];
		

		## Draw Eta: local indicator
		Eta<-rep(0,N);
		XA=X1%*%Alpha;

		j=1;
		for (i in nn){
			p=exp(XA[j,])*gden[j,];
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			

			Eta[i]=try(sample(c(1:K),size=1,replace=T,prob=p),silent=TRUE);
			if(length(grep("Error", Eta[i], fixed=TRUE))>0){
				Eta[i]=sample(c(1:K),size=1,replace=T,prob=rep(1/K,K))


			}
			j=j+1;
		}


		## Draw ALPHA (multiple-try MH)
		for (m in 1:M){
			
			obj=Draw_Alpha_MMH(Alpha[m,2:K],Alpha[,2:K],X1,Tau_sq[m],df,MA[m],Omega_star,SSA[m],SSAP[m],Eta[!Eta ==0],K,m,iter,burnIn)			
			Alpha[m,]=c(0,obj$par);					
			Alpha_accp[m]=obj$accp;
			print(paste("Alpha_m=",m));
			print(Alpha[m,]);

			
		}
						
			
		## Draw Tau_sq and a;
			for (m in 1:M){
				Tau_sq[m]=rigamma(1,alpha=0.5*(K+nu-1),beta=0.5*(Alpha[m,2:K]%*%Omega_star%*%Alpha[m,2:K])+nu/a[m])
				a[m]=rigamma(1,alpha=0.5*(nu+1),beta=nu/Tau_sq[m] + 1/A^2)
				
			}
			
			print("Tau_sq");print(Tau_sq);
			print("a");print(a);
			
		## Draw GAMMA			
			#Gamma draw: Multiple try MH
			objg=Draw_Gamma_log_M(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,MG,useCov)
			Gamma=objg$par;
			print("Gamma");
			print(Gamma);
		

		if(theoNULL !=1){
		## Draw SIGMA_SQ;
		   	Z0=abs(Z[Phi==0])
			Sigma_sq=rigamma(1,alpha=a0+length(Z0)/2,beta=b0+(Z0%*%Z0)/2);
			print(paste("Sigma_sq",Sigma_sq));
		}


		## Draw PHI : global indicator	
			f1<-NULL;
			f1=rowSums(exp((X%*%Alpha + log(Gden)))/rowSums(exp(X%*%Alpha)))
			
			log_P_phi=cbind(X%*%Gamma,0)
			log_P_phi[,1]=log_P_phi[,1]+ log(f1)

			if(theoNULL ==1){
				log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi)-0.5*Z^2;
			}else{
				log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi*Sigma_sq)-(1/(2*Sigma_sq))*Z^2;
			}


			P_phi=exp(log_P_phi)/apply(exp(log_P_phi),1,sum) 
			
			#none of the non-null SNPs should be <= 0.68;
			P_phi[abs(Z)<=mu,1]=0;P_phi[abs(Z)<=mu,2]=1;
			
			Phi_new=Phi #declare variable Phi_new, create a vector;
			for(i in 1:N){
				#in case NaN happens in P_phi[i,];
				P_phi[i,is.na(P_phi[i,])]=(1-sum(P_phi[i,!is.na(P_phi[i,])]))/sum(is.na(P_phi[i,]))

				Phi_new[i]=sample(c(1,0),size=1,replace=TRUE,prob=P_phi[i,])
			}
			Phi=Phi_new
			
			
			if(simulate == TRUE) print(paste("Phi matches Phi_true at current iteration",sum(Phi == Phi_true)/N))
			print(paste("Proprotion of the non-null cases at current iteration",sum(Phi)/N));
			

			

		## Save results after thin
				
			if(iter%%thin==0 & iter>=burnIn){
				array_ind=array_ind+1
				ALPHA_array[[array_ind]]=Alpha
				
				GAMMA_array[[array_ind]]=Gamma

				Tau_sq_array[[array_ind]]=Tau_sq
				a_array[[array_ind]]=a;
				
				if(theoNULL !=1){
					SIGMA_SQ_array[[array_ind]]=Sigma_sq
				}

			
				Accp_Rate_array_a[[array_ind]]=Alpha_accp;		
				Accp_Rate_array_g[[array_ind]]=objg$accp;
				
				for (m in 1:M){
					print(paste("Mean of Alpha m =",m));
					Alpha_mean[m,]=((array_ind-1)*Alpha_mean[m,]+ALPHA_array[[array_ind]][m,])/array_ind
					if(simulate==TRUE) {
						print(cbind(Alpha_mean[m,],Alpha_true[m,]))
					}else{
						print(Alpha_mean[m,])
					}
					print(paste("Average Multiple-try MH Accept Rate for Alpha m =",m,":",mean(sapply(Accp_Rate_array_a, '[', m))));
				}
				

				print("Tau_sq mean:");
				Tau_sq_mean=((array_ind-1)*Tau_sq_mean+Tau_sq_array[[array_ind]])/array_ind
				if(simulate==TRUE) {
					print(cbind(Tau_sq_mean,TAU_SQ_true))
				}else{
					print(Tau_sq_mean)
				}
				
				
				print("a mean:");
				a_mean=((array_ind-1)*a_mean+a_array[[array_ind]])/array_ind
				if(simulate==TRUE) {
					print(cbind(a_mean,a_true))
				}else{
					print(a_mean)
				}
				
				

				if(theoNULL !=1){
					print("Sigma_sq mean:");
					if(simulate==TRUE) {
						print(cbind(mean(as.numeric(SIGMA_SQ_array)),SIGMA_SQ_true))
					}else{
						print(mean(as.numeric(SIGMA_SQ_array)))
					}

				}

				print("Gamma mean:");
				Gamma_mean=((array_ind-1)*Gamma_mean+GAMMA_array[[array_ind]])/array_ind
				if(simulate==TRUE) {
					print(cbind(Gamma_mean,Gamma_true))
				}else{
					print(Gamma_mean)
				}
				print(paste("Multiple-try MH Accept Rate for Gamma (mean):",mean(as.numeric(Accp_Rate_array_g))));				
			
				if(simulate==TRUE) {
					PHI_match_rate=rbind(PHI_match_rate,sum(Phi==Phi_true)/N);
					print(paste("Phi matching rate (mean):",mean(PHI_match_rate)))
				}
				
				#probability of each SNP being Non-NULL, average of Phi over the iterations saved;
				Phi_mean=((array_ind-1)*Phi_mean+Phi)/array_ind;

			}

			
			if(iter%%(floor(nIter/4))==0 & iter>=burnIn & saveWhileRun==1){
				cnt=cnt+1;
				results_tmp=list()
				results_tmp[[1]]=ALPHA_array
				results_tmp[[2]]=GAMMA_array
				results_tmp[[3]]=Tau_sq_array
				if(simulate==TRUE) {

					results_tmp[[4]]=PHI_match_rate
				}
				results_tmp[[5]]=Accp_Rate_array_g
				results_tmp[[6]]=Alpha_mean
				results_tmp[[7]]=Gamma_mean
				results_tmp[[8]]=P_phi
				results_tmp[[9]]=Accp_Rate_array_a
				results_tmp[[10]]=Phi #last iteration Phi value
				results_tmp[[11]]=Phi_mean
				if (return_Bden ==1){
					results_tmp[[12]]=Gden
				}
				if(theoNULL !=1){
					results_tmp[[13]]=SIGMA_SQ_array
				}else{
					results_tmp[[13]]=1
				}
				results_tmp[[14]] = Eta #last iteration Eta value
				results_tmp[[15]]=a_array;
					
				save(file=paste("_",cnt,".R",sep=""),results_tmp,X,Z,N,nIter,burnIn,thin,simulate,mu,K,useCov)
				
				
										
			}

		

	}	

	#Calculate SD;
	#Alpah	
	for (m in 1:M){
		print(paste("SD of Alpha m = ",m));
		for(k in 1:K){
			print(sd(sapply(ALPHA_array,'[',m,k)))
			
		}
	}

	#Tau_sq
	print("Tau_sq SD:");
	for(i in 1:dim(X)[2]){
    		print(sd(as.numeric(unlist(lapply(Tau_sq_array, function(x) x[i])))))
	}	
	
	
	#a
	for(i in 1:dim(X)[2]){
    		print(sd(as.numeric(unlist(lapply(a_array, function(x) x[i])))))
	}

	if(theoNULL !=1){
		#Sigma_sq
		print(paste("Sigma SD:",sd(as.numeric(SIGMA_SQ_array))))
	}

	#Gamma
	print("Gamma SD:");
	for(i in 1:dim(X)[2]){
    		print(sd(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i])))))
	}


	#return results;
	results=list()
	results[[1]]=ALPHA_array
	results[[2]]=GAMMA_array
	results[[3]]=Tau_sq_array
	if(simulate==TRUE) {

		results[[4]]=PHI_match_rate
	}
	results[[5]]=Accp_Rate_array_g
	results[[6]]=Alpha_mean
	results[[7]]=Gamma_mean
	results[[8]]=P_phi
	results[[9]]=Accp_Rate_array_a
	results[[10]]=Phi #last iteration Phi value
	results[[11]]=Phi_mean
	if (return_Bden ==1){
		results[[12]]=Gden
	}
	if(theoNULL !=1){
		results[[13]]=SIGMA_SQ_array
	}else{
		results[[13]]=1
	}
	results[[14]] = Eta #last iteration Eta value
	results[[15]] = a_array;

	return(results)
}

