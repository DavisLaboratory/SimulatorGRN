#----graphGRN: functions----
validGRN<-function(object){
	if(length(nodes(object))!=0){
		#RNA maximum in range
		nodeRNAmax = unlist(nodeData(object,nodes(object),'rnamax'))
		if(sum(nodeRNAmax<0 | nodeRNAmax>1,na.rm=T)!=0){
			stop('RNA maximum expression has to be between 0 and 1')
		}
		
		#RNA degradation in range
		nodeRNAdeg = unlist(nodeData(object,nodes(object),'rnadg'))
		if(sum(nodeRNAdeg<0 | nodeRNAdeg>1,na.rm=T)!=0){
			stop('RNA degradation rate has to be between 0 and 1')
		}
		
		#Time constant in range
		nodetau = unlist(nodeData(object,nodes(object),'rnadg'))
		if(sum(nodetau<=0,na.rm=T)!=0){
			stop('Time constant must be positive')
		}
		
		#Interaction weight in range
		edgeweight = unlist(edgeData(object,from=nodes(object),attr='weight'))
		if(sum(edgeweight<0 | edgeweight>1,na.rm=T)!=0){
			stop('Interaction weight has to be between 0 and 1')
		}
		
		#EC50 in range
		edgeEC50 = unlist(edgeData(object,from=nodes(object),attr='EC50'))
		if(sum(edgeEC50<0 | edgeEC50>1,na.rm=T)!=0){
			stop('EC50 has to be between 0 and 1')
		}
		
		#Hill constant in range
		edgen = unlist(edgeData(object,from=nodes(object),attr='n'))
		if(sum(edgen==1,na.rm=T)!=0){
			stop('Hill constant cannot be 1')
		}
		
		#check for missing values
		nodetype = unlist(nodeData(object,nodes(object),'type'))
		if(sum(nodetype %in% 'OR' &
			   is.na(c(
			   	nodeRNAmax, nodeRNAdeg, nodetau, edgen, edgeEC50, edgeweight
			   )))!=0) {
			stop('Missing parameters not allowed for OR interactions')
		}
	}
	return(T)
}

initGRN<-function(.Object, ...){
	.Object = callNextMethod(.Object, ...)
	.Object@nodeParams = c('rnamax','rnadg','tau')
	.Object@edgeParams = c('weight','EC50','n')
	
	#set defaults for node and edge attributes
	nodeDataDefaults(.Object,'rnamax')=1
	nodeDataDefaults(.Object,'rnadg')=1
	nodeDataDefaults(.Object,'tau')=1
	nodeDataDefaults(.Object,'type')='OR'
	edgeDataDefaults(.Object,'weight')=1
	edgeDataDefaults(.Object,'EC50')=0.5
	edgeDataDefaults(.Object,'n')=1.39
	
	validObject(.Object)
	return(.Object)
}

nodeDataSetter<-function(self, n, attr, value) {
	self = callNextMethod()
	validObject(self)
	return(self)
}


edgeDataSetter<-function(self, from, to, attr, value) {
	self = callNextMethod()
	validObject(self)
	return(self)
}

#----graphGRN_AND: functions----
validGRN_AND<-function(object){
	if(length(nodes(object))!=0){
		edgeweight = unlist(edgeData(object,from=nodes(object),attr='weight'))
		edgeEC50 = unlist(edgeData(object,from=nodes(object),attr='EC50'))
		edgen = unlist(edgeData(object,from=nodes(object),attr='n'))
		
		#check for missing values
		nodetype = unlist(nodeData(object,nodes(object),'type'))
		if(sum(nodetype %in% 'AND' &
			   is.na(c(
			   	edgen, edgeEC50, edgeweight
			   )))!=0) {
			stop('Missing parameters not allowed for OR interactions')
		}
	}
	return(T)
}

initGRN_AND<-function(.Object, ...){
	.Object = callNextMethod(.Object, ...)
	
	#set defaults for node and edge attributes
	edgeDataDefaults(.Object,'type')='AND'
	
	validObject(.Object)
	return(.Object)
}

