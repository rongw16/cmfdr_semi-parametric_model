rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

log_p_alpha=function(Alpha,A,X,Tsq,O,E,K,m){
	

	Ind=matrix(0,nrow=length(E),ncol=K)
	
	for (i in 1: length(E)){
		Ind[i,E[i]]=1;
	}
	
	if (length(A)==length(Alpha)){ #intercept only
		log_p=sum(rowSums(Ind*(X%*%t(c(0,Alpha))- log(rowSums(exp(X%*%t(c(0,Alpha))))))))-0.5*(t(Alpha)%*%O%*%Alpha)/Tsq
	}else{
		 Xm=X[,m]; X1=X[,-m];

		if(dim(A)[1]==2){ #int and X1 only
			AA=cbind(0,t(A[-m,]));
		}else{ #more covariates;		
			AA=cbind(rep(0,(dim(A)[1]-1)),A[-m,]);
		}		
		
		log_p=sum(rowSums(Ind*((Xm%*%t(c(0,Alpha))+(X1%*%AA))- log(rowSums(exp(Xm%*%t(c(0,Alpha))+(X1%*%AA)))))))-0.5*(t(Alpha)%*%O%*%Alpha)/Tsq
		
	}
	return(log_p)
}



## Draw Alpha,log-scale MH alg:multiple-try;
Draw_Alpha_MMH=function(Alpha,A,X1,Tau_sq,df,Multiple,Omega_star,SSA,SSAP,Eta,K,m,iter,burnIn)
{					
		
	
		log_p_Alpha_star=rep(0,Multiple)
		log_p_Alpha_2star=rep(0,Multiple)
		p=rep(0,Multiple)
		den=0;num=0;

			#if(H==1){
				alpha_opt <-try(optim(Alpha,A=A,X=X1,Tsq=Tau_sq,O=Omega_star,E=Eta,K=K,m=m,log_p_alpha,method="Nelder-Mead",
					hessian=TRUE,control=list(maxit=10,fnscale=-1)),silent=TRUE)

				if(length(grep("Error", alpha_opt[[1]], fixed=TRUE))>0){
					return(list(par=rnorm(K-1,0,1), accp=0))

				}
				sigma=solve(-alpha_opt$hessian)

				if (sum(  diag(sigma) < 0  )> 0){
					sigma=matrix(0.1,nrow=K-1,ncol=K-1);
				}

				#print(paste("MH draw candidate: rmvt ( sigma=solve(-hessian)), Alpha_m=",m));	
				#print(sigma);
			#}else{
				#alpha_opt <-optim(Alpha,A=A,X=X1,Tsq=Tau_sq,O=Omega_star,E=Eta,K=K,m=m,log_p_alpha,method="Nelder-Mead",
				#hessian=FALSE,control=list(maxit=10,fnscale=-1))
				#sigma=diag(K-1)
				#sigma=matrix(0,nrow=K-1,ncol=K-1);
				#diag(sigma)=0.001;
			#}
			
	

			alpha=alpha_opt$par	
		
			if(iter > burnIn & sample(c(0,1),1,prob=c(1-SSAP,SSAP))==1){
				diag(sigma)=diag(sigma)*SSA
			}



			Alpha_star=t(rmvt(n=Multiple,alpha,sigma=sigma,df=df))
			
			for (i in 1:Multiple){
				log_p_Alpha_star[i]=log_p_alpha(Alpha_star[,i],A,X1,Tau_sq,Omega_star,Eta,K,m)
			}
			

			#control overfloat, -max(log_p_Alpha_star);
			p=exp(log_p_Alpha_star-max(log_p_Alpha_star))/sum(exp(log_p_Alpha_star - max(log_p_Alpha_star)))
			#in case there is still overfloat;
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			j=sample(c(1:Multiple),1,prob=p);							
			
			Alpha_2star=t(rmvt(n=Multiple-1,Alpha_star[,j],sigma=sigma,df=df))
			Alpha_2star <-cbind(Alpha_2star,Alpha)
			
			for (i in 1:Multiple){
				log_p_Alpha_2star[i]=log_p_alpha(Alpha_2star[,i],A,X1,Tau_sq,Omega_star,Eta,K,m)
			}
			
			#control overfloat
			num=sum(exp(log_p_Alpha_star -max(log_p_Alpha_star)))
			den=sum(exp(log_p_Alpha_2star -max(log_p_Alpha_star)))
											
			rho=min(1,num/den)
			
			#in case overfloat again
			if(is.na(rho)) {rho=0.5};
			
			accp=0;
			u=runif(1)
			if(u<rho){
				Alpha=Alpha_star[,j]
				accp=1;
			}
		
		
		return(list(par=Alpha, accp=accp))
}



log_p_gamma=function(Gamma,Z,X,SG,phi){
	
	log_p=sum(phi*(X%*%Gamma)-log(1+exp(X%*%Gamma)))-0.5*t(Gamma)%*%solve(SG)%*%Gamma

	
	return(log_p)
}

##Draw Gamma, log scale, multiple-try MH;
Draw_Gamma_log_M=function(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,Multiple,useCov)

{

		log_p_Gamma_star=rep(0,Multiple)
		log_p_Gamma_2star=rep(0,Multiple)
		p=rep(0,Multiple)
		
		
		
		if (useCov==1){
			gamma_opt <-optim(Gamma,Z=Z,X=X,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Nelder-Mead",
				hessian=TRUE,control=list(maxit=10,fnscale=-1))
		}else{
			gamma_opt <-optim(Gamma,Z=Z,X=X,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Brent",
				hessian=TRUE,lower=-1000-abs(Gamma),upper=1000+abs(Gamma),control=list(maxit=10,fnscale=-1))

		}
		gamma=gamma_opt$par
		sigma=solve(-gamma_opt$hessian)
		diag(sigma)=diag(sigma)*SSG
		

		Gamma_star=t(rmvt(n=Multiple,gamma,sigma=sigma,df=df))
			
			for (i in 1:Multiple){
				log_p_Gamma_star[i]=log_p_gamma(Gamma_star[,i],Z,X,Sigma_Gamma,Phi)

			}
			
			
			#control overfloat, -max(log_p_Gamma_star);
			p=exp(log_p_Gamma_star-max(log_p_Gamma_star))/sum(exp(log_p_Gamma_star - max(log_p_Gamma_star)))
			
			#in case there is still overfloat;
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			j=sample(c(1:Multiple),1,prob=p);	
							
			
			Gamma_2star=t(rmvt(n=Multiple-1,Gamma_star[,j],sigma=sigma,df=df))
			Gamma_2star <-cbind(Gamma_2star,Gamma)
			
			for (i in 1:Multiple){
				log_p_Gamma_2star[i]=log_p_gamma(Gamma_2star[,i],Z,X,Sigma_Gamma,Phi)
			}
			
											
			#control overfloat
			num=sum(exp(log_p_Gamma_star -max(log_p_Gamma_star)))
			den=sum(exp(log_p_Gamma_2star -max(log_p_Gamma_star)))
											
			rho=min(1,num/den)
			
			#in case overfloat again
			if(is.na(rho)) {rho=0.5};
			
			
			accp=0;
			u=runif(1)
			if(u<rho){
				Gamma=Gamma_star[,j]
				accp=1;
			
			}
		
		
		return(list(par=Gamma, accp=accp))
}

