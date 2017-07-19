library(nleqslv)
library(stringr)
library(plyr)
library(igraph)

fAct<-function(TF,EC50=0.5,n=1.39){
	B=(EC50^n-1)/(2*EC50^n-1)
	K_n=(B-1)
	act=B*TF^n/(K_n+TF^n)
	return(act)
}

model.init<-function(grn){
	
	#maximum expression values: checks and defaults
	if(is.null(grn$nodes$max)){
		grn$nodes['max']=1
	}else{
		if(sum(grn$nodes$max<0)!=0){
			stop('Maximum expression value of a gene must be >0')
		}
		
		if(sum(grn$nodes$max>1)!=0){
			warning('Maximum expression value >1 found, simulating overexpression data')
		}
	}
	
	#degradation rates: checks and defaults
	if(is.null(grn$nodes$dg)){
		grn$nodes['dg']=1
	}else{
		if(sum(grn$nodes$dg<0)!=0){
			stop('Degradation rate of a transcript must be >0')
		}
	}
	
	#time constant: checks and defaults
	if(is.null(grn$nodes$tau)){
		grn$nodes['tau']=1
	}else{
		if(sum(grn$nodes$tau<0)!=0){
			stop('Time constant must be >0')
		}
	}
	
	#identify and split/re-represent AND interactions
	andInts=which(grepl(',',grn$edges$regulator))
	cnames=colnames(grn$edges)
	cnames=cnames[cnames %in% c('regulator','interaction','EC50','n')]#colnames where splits are to be performed
	for (i in andInts){
		edge=sapply(grn$edges[i,],as.vector)
		splt=str_split(edge[cnames],',')
		names(splt)=cnames
		
		#check whether the format for AND interactions network is correct
		if(length(unique(laply(splt,length)))!=1){
			stop('Missing parameters for the AND interaction:\n\t',paste(edge,collapse=' | '))
		}
		
		#create temporary node and connections associated with it
		newnode=paste(c('AND_',splt$regulator),collapse='')
		splt=ldply(splt)
		rownames(splt)=splt$.id
		splt=as.data.frame(t(splt[,-1]),stringsAsFactors=F)
		# splt[colnames(grn$edges)[!colnames(grn$edges) %in% cnames]]=NA
		splt['target']=newnode
		splt=merge(splt,as.data.frame(t(c('regulator'=newnode,'interaction'='and',edge['target']))),all=T)
		if('weight' %in% colnames(grn$edges)){
			splt['weight']=NA
			splt[splt$regulator %in% newnode,'weight']=edge['weight']
		}
		
		#add new nodes and edges to existing data structure
		grn$nodes=rbind(grn$nodes,c('gene'=newnode,NA,NA,NA))
		grn$edges=rbind(grn$edges,splt)
	}
	#remove comma-formated entries of AND interactions
	grn$edges=grn$edges[-andInts,]
	rownames(grn$edges)=NULL
	#correct types of cols in data frame
	for(n in setdiff(colnames(grn$edges),c('regulator','interaction','target'))){
		grn$edges[,n]=as.numeric(grn$edges[,n])
	}
	
	#EC50: checks and defaults
	if(is.null(grn$edges$EC50)){
		grn$edges['EC50']=0.5
	}else{
		if(sum(grn$edges$EC50<0 | grn$edges$EC50>1,na.rm=T)!=0){
			stop('EC50 of an interaction must be in interval [0,1]')
		}
	}
	
	#Hill constant: checks and defaults
	if(is.null(grn$edges$n)){
		grn$edges['n']=1.39
	}else{
		if(sum(grn$edges$n==1,na.rm=T)!=0){
			stop('Hill constant of an interaction can not be 1')
		}
	}
	
	#Weights: checks and defaults
	if(is.null(grn$edges$weight)){
		grn$edges['weight']=1
	}else{
		if(sum(grn$edges$weight<0 | grn$edges$weight>1,na.rm=T)!=0){
			stop('Interaction weights must be in the interval [0,1]')
		}
	}
	
	#normalize weights
	warning('Weights will be normalized')
	
	#correct the parameters for AND interactions
	grn$edges[grepl('^AND_',grn$edges$regulator),c('EC50','n')]=NA
	grn$edges[grepl('^AND_',grn$edges$target),c('weight')]=NA
	
	#initialize formulae
	grn$nodes['formula']='0'
	
	return(grn)
}

grn.to.igrn<-function(grn){
	nodes=grn$nodes
	nodes$max=as.numeric(nodes$max)
	nodes$dg=as.numeric(nodes$dg)
	nodes$tau=as.numeric(nodes$tau)
	
	edges=grn$edges
	edges=edges[,c(1,3,2,4:ncol(edges))]
	map=c('act'=1,'rep'=-1,'and'=NA)
	edges$interaction=map[edges$interaction]
	
	igrn=graph_from_data_frame(edges,directed=TRUE,vertices=nodes)
	
	#as_data_frame(igrn,what='both')
	#x=adjacent_vertices(igrn,'AND_CD',mode='in')
	#x[[1]]$name
	
	return(igrn)
}

#need to remap interactions
igrn.to.grn<-function(igrn){
	grn=as_data_frame(igrn,what='both')
	names(grn)=c('nodes','edges')
	colnames(grn$edges)[1:2]=c('regulator','target')
	grn$edges=grn$edges[,c(1,3,2,4:ncol(grn$edges))]
	
	return(grn)
}

model.generate.equations<-function(igrn){
	nodes=V(igrn)
	nodes=nodes[!grepl('AND_',nodes$name)]
	for(i in 1:length(nodes)){
		v=nodes[i]
		fn=v$formula
		#if external condition, rate of change is zero therefore return 0
		if(degree(igrn,v,mode='in')!=0){
			vReg=adjacent_vertices(igrn,v,mode='in')[[1]]#get list of regulators
			andInt=grepl('AND_',vReg$name)
			
			#calculate the activation function
			for(j in 1:length(vReg)){
				vR=vReg[j]#retrieve regulator
				e=E(igrn,c(vR$name,v$name))#retrieve interaction (edge)
				if(andInt[j]){
					#in the AND case
					vComplex=adjacent_vertices(igrn,vR,mode='in')[[1]]#get list of regulators in regulatory complex
					fn=paste(c(fn,'+W_',sort(vComplex$name),v$name),collapse='')
					for(k in 1:length(vComplex)){
						vC=vComplex[k]#retrieve regulator
						eC=E(igrn,c(vC$name,vR$name))#retrieve interaction (edge)
						
						#in the OR case
						pname=paste(vC$name,v$name,sep='')
						act=paste('fAct(',vC$name,',EC50_',pname,',n_',pname,')',sep='')
						#if inhibitory interaction
						if(eC$interaction==-1)
							act=paste('(1-',act,')',sep='')
						fn=paste(fn,'*',act,sep='')
					}
				}else{
					#in the OR case
					pname=paste(vR$name,v$name,sep='')
					act=paste('fAct(',vR$name,',EC50_',pname,',n_',pname,')',sep='')
					#if inhibitory interaction
					if(e$interaction==-1)
						act=paste('(1-',act,')',sep='')
					fn=paste(fn,'+W_',pname,
							 '*',act,
							 sep='')
				}
			}
			
			fn=paste('min(',fn,',1)',sep='')
			fn=paste(fn,'*',v$name,'_MAX-','r_',v$name,'*',v$name,sep='')
			fn=paste('(',fn,')/t_',v$name,sep='')
		}
		V(igrn)[v$name]$formula=fn
	}
	return(igrn)
}

params.vectorize<-function(igrn){
	nodes=V(igrn)
	edges=E(igrn)
	nodes=nodes[!grepl('AND_',nodes$name)]
	parms=list()
	
	#node parameters
	for(i in 1:length(nodes)){
		v=nodes[i]
		parms[paste(v$name,'MAX',sep='_')]=v$max
		parms[paste('r',v$name,sep='_')]=v$dg
		parms[paste('t',v$name,sep='_')]=v$tau
		
		#interactions
		#if external condition, rate of change is zero therefore return 0
		if(degree(igrn,v,mode='in')!=0){
			vReg=adjacent_vertices(igrn,v,mode='in')[[1]]#get list of regulators
			andInt=grepl('AND_',vReg$name)
			
			#calculate the activation function
			for(j in 1:length(vReg)){
				vR=vReg[j]#retrieve regulator
				e=E(igrn,c(vR$name,v$name))#retrieve interaction (edge)
				if(andInt[j]){
					#in the AND case
					vComplex=adjacent_vertices(igrn,vR,mode='in')[[1]]#get list of regulators in regulatory complex
					parms[paste(c('W_',sort(vComplex$name),v$name),collapse='')]=e$weight
					for(k in 1:length(vComplex)){
						vC=vComplex[k]#retrieve regulator
						eC=E(igrn,c(vC$name,vR$name))#retrieve interaction (edge)
						
						#in the OR case
						pname=paste(vC$name,v$name,sep='')
						parms[paste('EC50_',pname,sep='')]=eC$EC50
						parms[paste('n_',pname,sep='')]=eC$n
					}
				}else{
					#in the OR case
					pname=paste(vR$name,v$name,sep='')
					parms[paste('W_',pname,sep='')]=e$weight
					parms[paste('EC50_',pname,sep='')]=e$EC50
					parms[paste('n_',pname,sep='')]=e$n
				}
			}
		}
	}
	return(parms)
}

model.get.function<-function(exprs,ext,igrn){
	#start function
	fn='function(exprs,parms){'
	
	#function body
	fn=paste(fn,'\tparms=c(parms,exprs)',sep='\n')
	fn=paste(fn,'\trates=exprs*0',sep='\n')
	fn=paste(fn,'\trates=with(parms,{',sep='\n')
	#start with: create equations
	nodes=V(igrn)
	for(i in names(exprs)){
		eqn=paste('\t\t','rates[\"',i,'\"]<-',nodes[i]$formula,sep='')
		fn=paste(fn,eqn,sep='\n')
	}
	#end with
	fn=paste(fn,'\t\treturn(rates)',sep='\n')
	fn=paste(fn,'\t})',sep='\n')
	fn=paste(fn,'\n\treturn(rates)',sep='\n')
	#end function
	fn=paste(fn,'}',sep='\n')
	return(eval(parse(text = fn)))
}

model.solve<-function(igrn,ext){
	#validity of ext input
	indegree=degree(igrn,mode='in')
	missinginputs=setdiff(names(indegree)[indegree==0],names(ext))
	if(length(missinginputs!=0)){
		stop('missing external input values: ',missinginputs)
	}
	
	#generate expr vector and initialize
	genes=names(indegree[indegree!=0])
	genes=genes[!grepl('^AND_',genes)]
	exprs=numeric(length(genes))+0.5
	names(exprs)=genes
	
	#generate params vector
	parms=model.vectorize(igrn)
	parms=c(parms,ext)
	
	#generate equations and get funcion
	igrn=model.generate.equations(igrn)
	fn=model.get.function(exprs,ext,igrn)
	
	#solve equations using nleqslv
	soln=nleqslv(exprs,fn,jac=NULL,parms,jacobian=T)
	exprs=c(soln$x,ext)
	exprs=exprs[sort(names(exprs))]
	soln=c('exprs'=list(exprs),soln)
	
	return(soln)
}

model.addNoise<-function(soln){
	
}

# grn=model.init(grn1)
# igrn=grn.to.igrn(grn)
# ext=c(0.5,0.3)
# names(ext)=LETTERS[1:2]
# res=model.solve(igrn,ext)















# model.function<-function(exprs,igrn){
# 	newexprs=exprs*0
# 	if(is.null(names(exprs))){
# 		stop('Input vector must be named with gene names')
# 	}
# 	edgeList=as_data_frame(igrn,what='edges')
#
# 	nodes=V(igrn)
# 	nodes=nodes[!grepl('AND_',nodes$name)]
# 	for(i in 1:length(nodes)){
# 		v=nodes[i]
# 		#if external condition, rate of change is zero therefore return 0
# 		if(degree(igrn,v,mode='in')==0){
# 			newexprs[v$name]=0
# 		}else{
# 			vReg=adjacent_vertices(igrn,v,mode='in')[[1]]#get list of regulators
# 			andInt=grepl('AND_',vReg$name)
#
# 			#calculate the activation function
# 			act=0
# 			for(j in 1:length(vReg)){
# 				vR=vReg[j]#retrieve regulator
# 				e=E(igrn,c(vR$name,v$name))#retrieve interaction (edge)
# 				if(andInt[j]){
# 					#in the AND case
# 					tmp=e$weight
# 					vComplex=adjacent_vertices(igrn,vR,mode='in')[[1]]#get list of regulators in regulatory complex
# 					for(k in 1:length(vComplex)){
# 						vC=vComplex[k]#retrieve regulator
# 						eC=E(igrn,c(vC$name,vR$name))#retrieve interaction (edge)
# 						tmp=tmp*(fAct(exprs[vC$name],eC$EC50,eC$n)*eC$interaction+
# 								 	(eC$interaction==-1))
# 					}
# 					act=act+tmp
# 				}else{
# 					#in the OR case
# 					tmp=e$weight*fAct(exprs[vR$name],e$EC50,e$n)*e$interaction+
# 						(e$interaction==-1)
# 					act=act+tmp
# 				}
# 			}
# 			act=min(act,1)
# 			newexprs[v$name]=(act*v$max-v$dg*exprs[v$name])/v$tau
# 		}
# 	}
# 	return(newexprs)
# }