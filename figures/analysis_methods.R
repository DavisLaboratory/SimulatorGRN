# liblist=c('diggit','stringr','plyr','EBcoexpress','foreach')

z.score<-function(emat,conditions){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#apply the Fisher transformation
	r1=cor(t(expr1))
	r2=cor(t(expr2))
	z1=atanh(r1)
	z2=atanh(r2)
	z=(z1-z2)/sqrt(1/(sum(conditions==1)-3)+1/(sum(conditions==2)-3))
	z[is.nan(z)]=0
	
	return(z)
}

z.score.sp<-function(emat,conditions){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#apply the Fisher transformation
	r1=cor(t(expr1), method = 'spearman')
	r2=cor(t(expr2), method = 'spearman')
	z1=atanh(r1)
	z2=atanh(r2)
	z=(z1-z2)/sqrt(1.06/(sum(conditions==1)-3)+1.06/(sum(conditions==2)-3))
	z[is.nan(z)]=0
	
	return(z)
}

magic.score<-function(emat,conditions){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#apply the Fisher transformation
	r1=cor(t(expr1))
	r2=cor(t(expr2))
	z1=atanh(r1)/sqrt(1/(sum(conditions==1)-3))
	z2=atanh(r2)/sqrt(1/(sum(conditions==2)-3))
	
	#compute magic score
	tgtsize=round(ncol(emat)/2) - 3
	expz1=exp(1/sqrt(tgtsize)*2*z1)
	expz2=exp(1/sqrt(tgtsize)*2*z2)
	IM=abs((expz1-1)/(expz1+1))-abs((expz2-1)/(expz2+1))
	IM[is.nan(IM)]=0
	
	return(IM)
}

ftgi.score<-function(emat,conditions){
	score<-foreach(i=rownames(emat),.combine=cbind) %:%
		foreach(j=rownames(emat),.combine=rbind) %dopar% {
			e1=as.numeric(emat[i,])
			e2=as.numeric(emat[j,])
			m1=glm(as.factor(conditions)~e1+e2,family=binomial(link='logit'))
			m2=glm(as.factor(conditions)~e1*e2,family=binomial(link='logit'))
			sc=-log10(anova(m1,m2)[2,4])
			return(sc)
		}
	rownames(score) = colnames(score) = rownames(emat)
	return(score)
}

diffcoex.score<-function(emat,conditions,beta=1,cor.method='pearson'){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#apply the Fisher transformation
	r1=cor(t(expr1),method=cor.method)
	r2=cor(t(expr2),method=cor.method)
	D=sqrt(0.5*abs(sign(r1)*r1^2-sign(r2)*r2^2))
	D=D^beta
	T=D%*%D+ncol(D)*D
	
	mins=matrix(rep(rowSums(D),ncol(D)),nrow=ncol(D))
	mins=pmin(mins,matrix(rep(colSums(D),each=ncol(D)),nrow=ncol(D)))
	T=1-(T/(mins+1-D))
	return(1-T)
}

ebcoexpress.score<-function(emat,conditions,rand.seed=36,plot=F){
	set.seed(rand.seed)
	pat=ebPatterns(c("1,1","1,2"))
	D=makeMyD(emat, conditions, useBWMC=FALSE)
	initHP=initializeHP(D, conditions)
	oout=ebCoexpressOneStep(D, conditions, pat, initHP)
	result1=oout$POSTPROBS
	ppbDC1=result1[,2]
	
	#diagnostic plots if required
	if(plot){
		par(mfrow=c(1,2))
		priorDiagnostic(D, conditions, oout, 1)
		priorDiagnostic(D, conditions, oout, 2)
		par(mfrow=c(1,1))
	}
	
	#convert to matrix
	scoremat=matrix(rep(0,nrow(emat)^2),nrow=nrow(emat))
	colnames(scoremat)=rownames(scoremat)=rownames(emat)
	corpairs=ldply(str_split(names(ppbDC1),'~'))
	corpairs['prob']=ppbDC1
	for(i in 1:nrow(corpairs)){
		g1=as.character(corpairs[i,1])
		g2=as.character(corpairs[i,2])
		scoremat[g1,g2]=scoremat[g2,g1]=corpairs[i,3]
	}
	
	return(scoremat)
}

cai.score<-function(emat,conditions,tau=1,thresh.method=c('soft','lasso','hard','score')[1],eta=1){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#compute correlation matrices
	cov1=cov(t(expr1))
	cov2=cov(t(expr2))
	r1=cov2cor(cov1)
	r2=cov2cor(cov2)
	
	#calculate xi parameter and lambda as in paper
	#condition 1
	ctrdexpr1=as.matrix(expr1-rowMeans(expr1)%*%t(rep(1,ncol(expr1))))
	xi1=(ctrdexpr1^2)%*%t(ctrdexpr1^2)
	xi1=xi1-2*cov1*(ctrdexpr1 %*% t(ctrdexpr1))
	xi1=xi1+ncol(expr1)*cov1^2
	xi1=xi1/(ncol(expr1)*(diag(cov1)%*%t(diag(cov1))))
	lambda1=tau*sqrt((log(nrow(expr1))/ncol(expr1)))
	diagsqrtxi1=sqrt(diag(xi1))%*%t(rep(1,ncol(xi1)))
	lambda1=lambda1*(sqrt(xi1)+abs(r1/2)*(diagsqrtxi1+t(diagsqrtxi1)))
	
	#condition 2
	ctrdexpr2=as.matrix(expr2-rowMeans(expr2)%*%t(rep(1,ncol(expr2))))
	xi2=(ctrdexpr2^2)%*%t(ctrdexpr2^2)
	xi2=xi2-2*cov2*(ctrdexpr2 %*% t(ctrdexpr2))
	xi2=xi2+ncol(expr2)*cov2^2
	xi2=xi2/(ncol(expr2)*(diag(cov2)%*%t(diag(cov2))))
	lambda2=tau*sqrt((log(nrow(expr2))/ncol(expr2)))
	diagsqrtxi2=sqrt(diag(xi2))%*%t(rep(1,ncol(xi2)))
	lambda2=lambda2*(sqrt(xi2)+abs(r2/2)*(diagsqrtxi2+t(diagsqrtxi2)))
	
	#apply data derived thresholds calculated above
	z=r1-r2
	lambda=lambda1+lambda2
	if(thresh.method %in% 'soft'){
		return(sign(z)*pmax(abs(z)-lambda,0))
	}
	if(thresh.method %in% 'hard'){
		return(z*(abs(z)>lambda))
	}
	if(thresh.method %in% 'lasso'){
		return(z*pmax(1-abs(lambda/z)^eta,0))
	}
	if(thresh.method %in% 'score'){
		return(abs(z)-lambda)
	}
}

mindy.score<-function(emat,conditions,ncores=1,rand.seed=36){
	set.seed(rand.seed)
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	m1=mutualInfo(t(expr1),cores=ncores)
	m2=mutualInfo(t(expr2),cores=ncores)
	
	scoremat=m1-m2
	scoremat[is.na(scoremat)]=0
	return(scoremat)
}

dicer.score<-function(emat,conditions){
	expr1=emat[,conditions==1]
	expr2=emat[,conditions==2]
	
	#apply the Fisher transformation
	r1=cor(t(expr1))
	r2=cor(t(expr2))
	mu1=mean(r1)
	mu2=mean(r2)
	var1=var(as.numeric(r1))
	var2=var(as.numeric(r2))
	tscore=((r1-r2)-(mu1-mu2))/sqrt(var1+var2)
	return(tscore)
}

