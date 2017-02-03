library(splines)

#if K,degree=3,beginning and end points are fixed, the B-splines would be the same
	Bdensity<-function(byStep=0.001,mu=0.68,K,N,Z,Gden){
	
	grid=seq(mu-byStep,max(abs(Z))+byStep,by=byStep)
	knots=seq(mu-byStep,max(abs(Z))+byStep,len=K) #equal space;

	Basis.mat=matrix(0,nrow=length(grid),ncol=K);

	Basis.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=3,intercept=F,
	Boundary.knots=c(knots[1],knots[length(knots)]))
				
	phi.mat=Basis.mat[,1:K]

	#normalization
	phi.mat[,3:(K-2)]=phi.mat[,3:(K-2)]/(knots[2]-knots[1])
	phi.mat[,1]=phi.mat[,1]*2/(knots[2]-knots[1])
	phi.mat[,2]=phi.mat[,2]*4/(3*(knots[2]-knots[1]))
	phi.mat[,(K-1)]=phi.mat[,(K-1)]*4/(3*(knots[2]-knots[1]))
	phi.mat[,K]=phi.mat[,K]*2/(knots[2]-knots[1]);

	
	#interpolate;
	for(i in 1:K){				
		Gden[,i]=approx(grid,phi.mat[,i],abs(Z),rule=2,ties="ordered")$y
	}
	return (Gden);
}
