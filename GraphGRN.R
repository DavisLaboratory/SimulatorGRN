source('simulatorS4-methods.R')

#----Node----
setClass(
	Class = 'Node',
	slots = list(
		name = 'character',
		tau = 'numeric',
		rnamax = 'numeric',
		rnadeg = 'numeric',
		inedges = 'list'
	)
)

setValidity('Node',validNode)

setMethod(
	f = 'initialize',
	signature = 'Node',
	definition = initNode
)

setMethod(
	f = '$',
	signature = 'Node',
	definition = function(x, name) {
		return(slot(x, name))
	}
)

setMethod(
	f = '$<-',
	signature = 'Node',
	definition = function(x, name, value) {
		slot(x, name)<-value
		validObject(x)
		return(x)
	}
)

setMethod(
	f = 'show',
	signature = 'Node',
	definition = function(object) {
		print(paste('Name: ', object@name, sep = ''))
		print(paste('Time const: ', object@tau, sep = ''))
		print(paste('RNA max: ', object@rnamax, sep = ''))
		print(paste('RNA degradation rate: ', object@rnadeg, sep = ''))
	}
)

setGeneric(
	name = 'generateEqn',
	def = function(object){
		standardGeneric('generateEqn')
	}
)

setMethod(
	f = 'generateEqn',
	signature = 'Node',
	definition = generateRateEqn
)

#----Edge----
setClass(
	Class = 'Edge',
	slots = list(
		from = 'list',
		to = 'Node',
		weight = 'numeric',
		EC50 = 'numeric',
		n = 'numeric',
		activation = 'logical',
		name = 'character'
	)
)

setValidity('Edge',validEdge)

setMethod(
	f = 'initialize',
	signature = 'Edge',
	definition = initEdge
)

setMethod(
	f = '$',
	signature = 'Edge',
	definition = function(x, name) {
		return(slot(x, name))
	}
)

setMethod(
	f = '$<-',
	signature = 'Edge',
	definition = function(x, name, value) {
		if (name %in% 'name') {
			stop('Edge name is generated automatically. It cannot be modified.')
		}
		
		slot(x, name)<-value
		validObject(x)
		return(x)
	}
)

setMethod(
	f = 'show',
	signature = 'Edge',
	definition = function(object) {
		return(NULL)
	}
)

setGeneric(
	name = 'generateActivationEqn',
	def = function(object){
		standardGeneric('generateActivationEqn')
	}
)

#----EdgeOr----
setClass(
	Class = 'EdgeOr',
	contains = 'Edge'
)

setValidity('EdgeOr',validEdgeOr)

setMethod(
	f = 'generateActivationEqn',
	signature = 'EdgeOr',
	definition = generateActivationEqnOr
)

#----EdgeAnd----
setClass(
	Class = 'EdgeAnd',
	contains = 'Edge'
)

setValidity('EdgeAnd',validEdgeAnd)

setMethod(
	f = 'initialize',
	signature = 'EdgeAnd',
	definition = initEdgeAnd
)

setMethod(
	f = 'generateActivationEqn',
	signature = 'EdgeAnd',
	definition = generateActivationEqnAnd
)

#----GraphGRN----
setClass(
	Class = 'GraphGRN',
	slots = list(
		nodeset = 'list',
		edgeset = 'list'
	)
)

setValidity('GraphGRN', validGraphGRN)

setMethod(
	f = 'initialize',
	signature = 'GraphGRN',
	definition = initGraphGRN
)

#----GraphGRN:addNode----
setGeneric(
	name = 'addNode',
	def = function(graph, node, tau, rnamax, rnadeg, inedges) {
		standardGeneric('addNode')
	}
)

setMethod(
	f = 'addNode',
	signature = c('GraphGRN', 'Node', 'missing', 'missing', 'missing', 'missing'),
	definition = function(graph, node, tau, rnamax, rnadeg, inedges) {
		graph@nodeset = c(graph@nodeset, node)
		
		#named entry to graph structure
		names(graph@nodeset)[length(graph@nodeset)] = node$name
		validObject(graph)
		return(graph)
	}
)

setMethod(
	f = 'addNode',
	signature = c('GraphGRN', 'character', 'ANY', 'ANY', 'ANY', 'missing'),
	definition = function(graph, node, tau, rnamax, rnadeg, inedges) {
		#create default node
		nodeObj = new('Node', name = node)
		
		#modify default node with provided parameters
		if(!missing(tau))
			edgeObj$tau = tau
		if(!missing(rnamax))
			edgeObj$rnamax = rnamax
		if(!missing(rnadeg))
			edgeObj$rnadeg = rnadeg
		
		graph = addNode(graph, nodeObj)
		
		return(graph)
	}
)

#----GraphGRN:addEdge----
setGeneric(
	name = 'addEdge',
	def = function(graph, from, to, edgetype, activation, weight, EC50, n) {
		standardGeneric('addEdge')
	}
)

setMethod(
	f = 'addEdge',
	signature = c('GraphGRN', 'character', 'character', 'ANY', 'ANY', 'ANY', 'ANY', 'ANY'),
	definition = function(graph, from, to, edgetype, activation, weight, EC50, n) {
		#retrieve nodes from graph
		tonode = getNode(graph, to)
		fromnode = list()
		for (i in 1:length(from)){
			fromnode = c(fromnode, getNode(graph, from[i]))
		}
		
		#get class
		if(missing(edgetype) || edgetype %in% 'or')
			edgeclass = 'EdgeOr'
		else if(edgetype %in% 'and')
			edgeclass = 'EdgeAnd'
		else
			stop('Unrecognised edgetype')
		
		#create default edge
		edgeObj = new(edgeclass, from = c(fromnode), to = tonode)
		#modify default edge with provided parameters
		if(!missing(activation))
			edgeObj$activation = activation
		if(!missing(weight))
			edgeObj$weight = weight
		if(!missing(EC50))
			edgeObj$EC50 = EC50
		if(!missing(n))
			edgeObj$n = n
		
		#add edge to graph
		graph@edgeset = c(graph@edgeset, edgeObj)
		
		#update node inedges information
		tonode$inedges = c(tonode$inedges, edgeObj)
		graph@nodeset[[tonode$name]] = tonode
		
		#named entry to graph structure
		names(graph@edgeset)[length(graph@edgeset)] = edgeObj$name
		validObject(graph)
		
		return(graph)
	}
)

#----GraphGRN:getEdge----
setGeneric(
	name = 'getEdge',
	def = function(graph, from, to) {
		standardGeneric('getEdge')
	}
)

setGeneric(
	name = 'getEdge<-',
	def = function(graph, from, to, value) {
		standardGeneric('getEdge<-')
	}
)

setMethod(
	f = 'getEdge',
	signature = c('GraphGRN', 'character', 'character'),
	definition = function(graph, from, to) {
		#generate name
		edgename = paste(from, collapse = '')
		edgename = paste(edgename, to, sep = '->')
		
		edgeObj = graph@edgeset[[edgename]]
		return(edgeObj)
	}
)

setReplaceMethod(
	f = 'getEdge',
	signature = c('GraphGRN', 'character', 'character', 'Edge'),
	definition = function(graph, from, to, value) {
		#generate name
		edgename = paste(from, collapse = '')
		edgename = paste(edgename, to, sep = '->')
		
		graph@edgeset[[edgename]] = value
		return(graph)
	}
)

#----GraphGRN:getInputNode----
setGeneric(
	name = 'getInputNodes',
	def = function(graph, nodename) {
		standardGeneric('getInputNodes')
	}
)

setMethod(
	f = 'getInputNodes',
	signature = c('GraphGRN'),
	definition = function(graph) {
		inputnodes = sapply(graph@nodeset, function(x) {
			nodename = NA
			if (length(x$inedges) == 0) {
				nodename = x$name
			}
			return(nodename)
		})
		
		inputnodes = inputnodes[!is.na(inputnodes)]
		names(inputnodes) = NULL
		
		return(inputnodes)
	}
)

#----GraphGRN:getNode----
setGeneric(
	name = 'getNode',
	def = function(graph, nodename) {
		standardGeneric('getNode')
	}
)

setGeneric(
	name = 'getNode<-',
	def = function(graph, nodename, value) {
		standardGeneric('getNode<-')
	}
)

setMethod(
	f = 'getNode',
	signature = c('GraphGRN', 'character'),
	definition = function(graph, nodename) {
		nodeObj = graph@nodeset[[nodename]]
		return(nodeObj)
	}
)

setReplaceMethod(
	f = 'getNode',
	signature = c('GraphGRN', 'character', 'Node'),
	definition = function(graph, nodename, value) {
		graph@nodeset[[nodename]] = value
		return(graph)
	}
)

#----GraphGRN: generateODE----
setGeneric(
	name = 'generateODE',
	def = function(graph) {
		standardGeneric('generateODE')
	}
)

setMethod(
	f = 'generateODE',
	signature = c('GraphGRN'),
	definition = function(graph) {
		fn = 'function(exprs, graph, externalInputs) {'
		
		#define the activation function
		fn = paste(fn, '\tfAct <- function(TF, EC50 = 0.5, n = 1.39) {', sep = '\n')
		fn = paste(fn, '\t\tB = (EC50 ^ n - 1) / (2 * EC50 ^ n - 1)', sep = '\n')
		fn = paste(fn, '\t\tK_n = (B - 1)', sep = '\n')
		fn = paste(fn, '\t\tact = B * TF ^ n / (K_n + TF ^ n)', sep = '\n')
		fn = paste(fn, '\t\t', sep = '\n')
		fn = paste(fn, '\t\treturn(act)', sep = '\n')
		fn = paste(fn, '\t}', sep = '\n')
		fn = paste(fn, '\t', sep = '\n')
		
		#function body
		fn = paste(fn, '\tparms = c(list(), externalInputs, exprs)', sep = '\n')
		fn = paste(fn, '\trates = exprs * 0', sep = '\n')
		fn = paste(fn, '\trates = with(parms, {', sep = '\n')
		#start with: create equations
		inputNodes = getInputNodes(graph)
		for (node in graph@nodeset) {
			if(node$name %in% inputNodes)
				next
			eqn = paste('\t\t', 'rates[\"', node$name, '\"] = ', generateRateEqn(node), sep = '')
			fn = paste(fn, eqn, sep = '\n')
		}
		
		#end with
		fn = paste(fn, '\t\treturn(rates)', sep = '\n')
		fn = paste(fn, '\t})', sep = '\n')
		fn = paste(fn, '\n\treturn(rates)', sep = '\n')
		#end function
		fn = paste(fn, '}', sep = '\n')
		
		return(eval(parse(text = fn)))
	}
)




