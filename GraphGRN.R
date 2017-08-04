#----Node----
setClass(
	Class = 'Node',
	slots = list(
		name = 'character',
		spmax = 'numeric',
		spdeg = 'numeric',
		inedges = 'list'
	)
)
setValidity('Node', validNode)

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
    mxp = 10
    
    cat(paste0('***', class(object), '***'), '\n')
    cat('Name:', object@name, '\n')
    cat('RNA max:', object@spmax, '\n')
    cat('RNA degradation rate:', object@spdeg, '\n')
    
    #print edge list
    edgelist = sapply(object@inedges, function(x) x$name)
    cat('Incoming edges:', maxPrint(edgelist, mxp), '\n')
  }
)

#----NodeRNA----
setClass(
	Class = 'NodeRNA',
	contains = 'Node',
	slots = list(
		tau = 'numeric'
	)
)
setValidity('NodeRNA', validNodeRNA)

setMethod(
	f = 'initialize',
	signature = 'NodeRNA',
	definition = initNodeRNA
)

setMethod(
	f = 'show',
	signature = 'NodeRNA',
	definition = function(object) {
	  callNextMethod()
		cat('Time const:', object@tau, '\n')
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
	signature = 'NodeRNA',
	definition = generateRateEqn
)

#----Edge----
setClass(
	Class = 'Edge',
	slots = list(
		from = 'list',
		to = 'Node',
		weight = 'numeric',
		name = 'character'
	)
)

setValidity('Edge', validEdge)

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
		mxp = 10
		
		cat(paste0('***', class(object), '***'), '\n')
		cat('Name:', object@name, '\n')
		cat('From:', maxPrint(sapply(object@from, function(x) x$name), mxp), '\n')
		cat('To:', object@to$name, '\n')
		cat('Weight:', object@weight, '\n')
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
	contains = 'Edge',
	slots = list(
	  EC50 = 'numeric',
	  n = 'numeric',
	  activation = 'logical'
	)
)

setValidity('EdgeOr', validEdgeOr)

setMethod(
  f = 'initialize',
  signature = 'EdgeOr',
  definition = initEdgeOr
)

setMethod(
	f = 'generateActivationEqn',
	signature = 'EdgeOr',
	definition = generateActivationEqnOr
)

setMethod(
  f = 'show',
  signature = 'EdgeOr',
  definition = function(object) {
    mxp = 10
    
    callNextMethod()
    cat('EC50:', maxPrint(object@EC50, mxp), '\n')
    cat('Hill constant (n):', maxPrint(object@n, mxp), '\n')
    cat('Activation:', maxPrint(object@activation, mxp), '\n')  }
)

#----EdgeAnd----
setClass(
	Class = 'EdgeAnd',
	contains = 'Edge',
	slots = list(
	  EC50 = 'numeric',
	  n = 'numeric',
	  activation = 'logical'
	)
)

setValidity('EdgeAnd', validEdgeAnd)

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

setMethod(
  f = 'show',
  signature = 'EdgeAnd',
  definition = function(object) {
    mxp = 10
    
    callNextMethod()
    cat('EC50:', maxPrint(object@EC50, mxp), '\n')
    cat('Hill constant (n):', maxPrint(object@n, mxp), '\n')
    cat('Activation:', maxPrint(object@activation, mxp), '\n')  }
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

setMethod(
	f = 'show',
	signature = 'GraphGRN',
	definition = function(object) {
		mxp = 10
		nodes = object@nodeset
		edges = object@edgeset
		cat('GraphGRN object with', length(nodes), 'nodes and',
			length(edges), 'edges','\n')
		nodenames = sapply(nodes[1:min(mxp+1, length(nodes))], 
						   function(x) x$name)
		edgenames = sapply(edges[1:min(mxp+1, length(edges))], 
						   function(x) x$name)
		if (length(nodes) == 0) {
		  cat('Nodes: (0)', '\n')
		} else{
		  cat('Nodes:', maxPrint(nodenames, mxp, length(nodes)), '\n')
		}
		if (length(edges) == 0) {
		  cat('Edges: (0)', '\n')
		} else{
		  cat('Edges:', maxPrint(edgenames, mxp, length(edges)), '\n')
		}
	}
)

#----GraphGRN:addNode----
setGeneric(
	name = 'addNode',
	def = function(graph, node, tau, spmax, spdeg, inedges) {
		standardGeneric('addNode')
	}
)

setMethod(
	f = 'addNode',
	signature = c('GraphGRN', 'NodeRNA', 'missing', 'missing', 'missing', 'missing'),
	definition = function(graph, node, tau, spmax, spdeg, inedges) {
		graph@nodeset = c(graph@nodeset, node)
		
		#node name is not empty
		if (node$name %in% '') {
			stop('Node name cannot be empty')
		}
		
		#named entry to graph structure
		names(graph@nodeset)[length(graph@nodeset)] = node$name
		validObject(graph)
		return(graph)
	}
)

setMethod(
	f = 'addNode',
	signature = c('GraphGRN', 'character', 'ANY', 'ANY', 'ANY', 'missing'),
	definition = function(graph, node, tau, spmax, spdeg, inedges) {
		#create default node
		nodeObj = new('NodeRNA', name = node)
		
		#modify default node with provided parameters
		if(!missing(tau) && !is.null(tau))
			edgeObj$tau = tau
		if(!missing(spmax) && !is.null(spmax))
			edgeObj$spmax = spmax
		if(!missing(spdeg) && !is.null(spdeg))
			edgeObj$spdeg = spdeg
		
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
		if(!missing(activation) && !is.null(activation))
			edgeObj$activation = activation
		if(!missing(weight) && !is.null(weight))
			edgeObj$weight = weight
		if(!missing(EC50) && !is.null(EC50))
			edgeObj$EC50 = EC50
		if(!missing(n) && !is.null(n))
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

#----GraphGRN:getInputNodes----
setGeneric(
	name = 'getInputNodes',
	def = function(graph) {
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

#----GraphGRN:nodenames----
setGeneric(
	name = 'nodenames',
	def = function(graph) {
		standardGeneric('nodenames')
	}
)

setMethod(
	f = 'nodenames',
	signature = c('GraphGRN'),
	definition = function(graph) {
		return(names(graph@nodeset))
	}
)

#----GraphGRN:edgenames----
setGeneric(
	name = 'edgenames',
	def = function(graph) {
		standardGeneric('edgenames')
	}
)

setMethod(
	f = 'edgenames',
	signature = c('GraphGRN'),
	definition = function(graph) {
		return(names(graph@edgeset))
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
		fn = 'function(exprs, externalInputs, graph) {'
		
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

#----GraphGRN: sample----


