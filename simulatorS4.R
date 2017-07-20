library(graph)
source('simulatorS4-methods.R')

#----create graphGRN class----
setClass(
	Class = 'graphGRN',
	contains = 'graphNEL',
	slots = list(nodeParams = 'character',
				 edgeParams = 'character')
)

setValidity('graphGRN', validGRN)

setMethod(
	f = 'initialize',
	signature = 'graphGRN',
	definition = initGRN
)

setMethod(
	f = 'nodeData<-',
	signature = c(self = 'graphGRN', n = 'character', attr = 'character'),
	definition = nodeDataSetter
)

setMethod(
	f = 'edgeData<-',
	signature = c(
		self = 'graphGRN',
		from = 'character',
		to = 'character',
		attr = 'character'
	),
	definition = edgeDataSetter
)

setGeneric(
	name = 'addComplexInteraction',
	def = function(node) {
		standardGeneric('addComplexInteraction')
	}
)

setGeneric(
	name = 'generateEqn',
	def = function(node) {
		standardGeneric('generateEqn')
	}
)

#----create AND interaction class----
setClass(
	Class = 'graphGRN_AND',
	contains = 'graphGRN'
)

setValidity('graphGRN_AND', validGRN_AND)

setMethod(
	f = 'initialize',
	signature = 'graphGRN_AND',
	definition = initGRN
)

setMethod(
	f = 'nodeData<-',
	signature = c(self = 'graphGRN_AND', n = 'character', attr = 'character'),
	definition = nodeDataSetter
)

setMethod(
	f = 'edgeData<-',
	signature = c(
		self = 'graphGRN_AND',
		from = 'character',
		to = 'character',
		attr = 'character'
	),
	definition = edgeDataSetter
)





graphGRN<-function(nodes = character(), edgeL = list()){
	grn=new(Class = 'graphGRN',nodes=nodes, edgeL=edgeL, edgemode='directed')
	return(grn)
}

# graphGRN<-function(nodedf = data.frame(),edgedf = data.frame()){
# 	nodeL=nodedf$gene
# 	#add and nodes
# 	andnodes=edgedf$regulator[edges$regulator]
# 	
# 	grn=new(Class = 'graphGRN',nodes=nodedf$gene, edgeL=edgeL, graphData=list('edgemode'='directed'))
# }

#----use cases----
grn=graphGRN()
grn1=addNode(LETTERS[1:4],grn)
grn1=addEdge('A','B',grn1)
grn1=addEdge('A','C',grn1)

# setClass(Class = 'A', slots = list(n='numeric'))
# setValidity('A',function(object){if(object@n>100)stop('Error @ A')})
# inc=function(object){object@n=object@n+1;validObject(object);return(object)}
# setGeneric('inc')
# setClass(Class = 'B', contains = 'A')
# setValidity('B',function(object){if(object@n>100)stop('Error @ B')})
