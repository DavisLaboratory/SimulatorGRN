#----Node: functions----
validNode <- function(object) {
	#RNA maximum in range
	if (!is.na(object@spmax) & (object@spmax < 0 | object@spmax > 1)) {
		stop('RNA maximum expression has to be between 0 and 1')
	}
	
	#RNA degradation in range
	if (!is.na(object@spdeg) & (object@spdeg < 0 | object@spdeg > 1)) {
		stop('RNA degradation rate has to be between 0 and 1')
	}
	
	#Incoming edges
	if (!all(sapply(object@inedges, is, 'Edge'))) {
		stop('All inedges must be valid Edge objects')
	}
	
	#Incoming edges: check that all edges have to as this node
	testin = sapply(object@inedges, function(x) {
		return(identical(x$to$name, object$name))
	})
	
	if (!all(testin)) {
		stop('All interactions must have the \'to\' node as the current node')
	}
	
	#check that all regulators are unique
	regnames = sapply(object@inedges, function(x) {
		sapply(x$from, slot, 'name')
	})
	regnames = as.vector(regnames)
	if (length(unique(regnames)) < length(regnames)) {
		stop('All regulators must be unique')
	}
	
	return(TRUE)
}

initNode <- function(.Object, ..., name = '', spmax = 1, spdeg = 1, inedges = list()) {
	.Object@name = name
	.Object@spmax = spmax
	.Object@spdeg = spdeg
	.Object@inedges = inedges
	
	validObject(.Object)
	return(.Object)
}

#----NodeRNA: functions----
validNodeRNA <- function(object){
	#Time constant in range
	if (!is.na(object@tau) & (object@tau <= 0)) {
		stop('Time constant must be positive')
	}
	
	return(TRUE)
}

initNodeRNA <- function(.Object, ..., tau = 1) {
	.Object@tau = tau
	.Object = callNextMethod()
	
	validObject(.Object)
	return(.Object)
}

generateRateEqn <- function(object) {
	inedges = object$inedges
	#no rate equations for input nodes
	if (length(inedges) == 0) {
		return('')
	}
	
	#generate activation functions for each interaction
	actEqns = sapply(inedges, generateActivationEqn)
	act = paste(actEqns, collapse = ' + ')
	
	#subtract the combinations
	nEq = length(actEqns)
	if (nEq > 1) {
		combEqns = character(0)
		for (i in 2:nEq) {
			combs = combn(actEqns, i)
			combEqn = apply(combs, 2, paste, collapse = ' * ')
			
			sgn = ' - '
			if (i %% 2 != 0){
			  sgn = ' + '
			}
			
			act = paste(c(act, combEqn), collapse = sgn)
		}
	}
	
	#generate rate equation
	rateEqn = paste('(', act, ')', sep = '')
	rateEqn = paste(rateEqn, object@spmax, sep = ' * ')
	degradationEqn = paste(object@spdeg, object@name, sep = ' * ')
	rateEqn = paste(rateEqn, degradationEqn, sep = ' - ')
	return(rateEqn)
}

#----Edge: functions----
validEdge <- function(object) {
	#weight in range
	if (any(is.na(object@weight)) | sum(object@weight < 0 | object@weight > 1) > 0) {
		stop('Interaction weight has to be between 0 and 1')
	}
  
  #from are all of class Node
  if (!all(sapply(object@from, is, 'Node'))) {
    stop('All from nodes must be of class \'Node\'')
  }
}

validActivationParams <- function(object) {
	#EC50 in range
	if (any(is.na(object@EC50)) | sum(object@EC50 < 0 | object@EC50 > 1) > 0) {
		stop('EC50 has to be between 0 and 1')
	}
	
	#Hill constant in range
	if (any(is.na(object@n)) | sum(object@n == 1) > 0) {
		stop('Hill constant (n) cannot be 1')
	}
	
	return(TRUE)
}

#----EdgeOr: functions----
validEdgeOr <- function(object) {
  validActivationParams(object)
  
	#Weight length is 1
	if (length(object@weight) != 1) {
		stop('Only 1 parameter for the weight should be provided for OR interactions')
	}
	
	#EC50 length is 1
	if (length(object@EC50) != 1) {
		stop('Only 1 parameter for the EC50 should be provided for OR interactions')
	}
	
	#Hill constant length is 1
	if (length(object@n) != 1) {
		stop('Only 1 parameter for the Hill constant(n) should be provided for OR interactions')
	}
	
	#Number of regulators is 1
	if (length(object@from) != 1) {
		stop('Only 1 input node should be provided for OR interactions')
	}
	
	return(TRUE)
}

initEdgeOr <- function(.Object, ..., from, to, weight = 1, EC50 = 0.5, n = 1.39, activation = T) {
  .Object@from = from
  .Object@to = to
  .Object@weight = weight
  .Object@EC50 = EC50
  .Object@n = n
  .Object@activation = activation
  
  #generate name
  name = sapply(from, function(x) x$name)
  name = paste(name, collapse = '')
  name = paste(name, to$name, sep = '->')
  .Object@name = name
  
  validObject(.Object)
  return(.Object)
}

generateActivationEqnOr <- function(object) {
	e = object
	#generate activation eqn
	act = paste('fAct(', e$from[[1]]$name, ',', e$EC50, ',', e$n, ')', sep =	'')
	if (!e$activation) {
		act = paste('(1-', act, ')', sep = '')
	}
	act = paste(e$weight, act, sep = ' * ')
	
	return(act)
}

#----EdgeAnd: functions----
validEdgeAnd <- function(object) {
  validActivationParams(object)
	numint = length(object@from) #number of interactors
	
	#number of inputs
	if(numint == 1){
		stop('AND interaction requires at least 2 regulators')
	}
	
	#Weight length is 1
	if (length(object@weight) != 1) {
		stop('Only 1 parameter for the weight should be provided for AND interactions')
	}
	
	#EC50 length matches number of interactors
	if (length(object@EC50) != numint) {
		stop('Missing EC50 parameters for the AND interaction')
	}
	
	#Hill constant length matches number of interactors
	if (length(object@n) != numint) {
		stop('Missing Hill constant(n) parameters for the AND interaction')
	}
	
	#Hill constant length matches number of interactors
	if (length(object@activation) != numint) {
		stop('Missing activation/repression status for the AND interaction')
	}
	
	return(TRUE)
}

initEdgeAnd <- function(.Object, ..., from, to, weight = 1, EC50 = c(), n = c(), activation = c()) {
	numint = length(from)
	if (length(EC50) == 0)
		EC50 = rep(0.5, numint)
	if (length(n) == 0)
		n = rep(1.39, numint)
	if (length(activation) == 0)
		activation = rep(T, numint)
	
	.Object@from = from
	.Object@to = to
	.Object@weight = weight
	.Object@EC50 = EC50
	.Object@n = n
	.Object@activation = activation
	
	#generate name
	name = sapply(from, function(x) x$name)
	name = paste(name, collapse = '')
	name = paste(name, to$name, sep = '->')
	.Object@name = name
	
	validObject(.Object)
	return(.Object)
}

generateActivationEqnAnd <- function(object) {
	e = object
	
	#generate AND activation eqn
	act = character(length(e@from))
	for (i in 1:length(act)) {
		act[i] = paste('fAct(', e@from[[i]]$name, ',', e@EC50[[i]], ',', e@n[[i]], ')', sep =	'')
		if (!e$activation[[i]]) {
			act[i] = paste('(1-', act[i], ')', sep = '')
		}
	}
	
	act = paste(c(e$weight, act), collapse = ' * ')
	
	return(act)
}

#----GraphGRN: functions----
validGraphGRN <- function(object) {
	#nodeset are all of class Node
	if (!all(sapply(object@nodeset, is, 'Node'))) {
		stop('All nodes must be of class \'Node\'')
	}
	
	#check names of nodeset
	if (!all(sapply(object@nodeset, function(x) x$name) == names(object@nodeset))){
		stop('Invalid graph generated. Use the \'addNode\' method to add a node to the graph.')
	}
	
	#nodeset are all of class Node
	if (!all(sapply(object@edgeset, is, 'Edge'))) {
		stop('All nodes must be of class \'Edge\'')
	}
	
	#check names of nodeset
	if (!all(sapply(object@edgeset, function(x) x$name) == names(object@edgeset))){
		stop('Invalid graph generated. Use the \'addEdge\' method to add an edge to the graph.')
	}
	
	return(TRUE)
}

initGraphGRN <- function(.Object, ..., nodeset = list(), edgeset = list()) {
	.Object@nodeset = nodeset
	.Object@edgeset = edgeset
	
	validObject(.Object)
	
	return(.Object)
}

addEdgeHelper<-function(graph, edge){
	graph@edgeset = c(graph@edgeset, edge)
	
	#update node inedges information
	to$inedges = c(to$inedges, edge)
	graph@nodeset[[to$name]] = to
	validObject(graph)
	
	#named entry to graph structure
	names(graph@edgeset)[length(graph@edgeset)] = edge$name
	validObject(graph)
	
	return(graph)
}

getEdgeClass<-function(edgetype){
	if(edgetype %in% 'or')
		edgeclass = 'EdgeOr'
	else if(edgetype %in% 'and')
		edgeclass = 'EdgeAnd'
	else
		stop('Unrecognised edgetype')
	
	return(edgeclass)
}

#----GraphGRN: Conversion functions----
dfToGraphGRN <- function(edges, nodes) {
  if(missing(nodes) || is.null(nodes)){
    nodes = data.frame('node' = unique(c(edges[ , 1], edges[ , 3])), stringsAsFactors = F)
  }
  
  #names for the main columns should be consistent
  colnames(edges)[1:3] = c('from', 'activation', 'to')
  colnames(nodes)[1] = c('node')
  
  #missing interaction type
  if (!'type' %in% colnames(edges)) {
    edges['type'] = 'or'
  }
  
  #create graph object
  grn = new('GraphGRN')
  #add nodes
  for (i in 1:nrow(nodes)){
    n = nodes[i, , drop = F]
    grn = addNode(grn, n$node, n$tau, n$max, n$deg)
  }
  
  #add edges
  for (i in 1:nrow(edges)){
    e = edges[i, , drop = F]
    grn = addEdge(grn, e$from, e$to, e$type, e$activation, e$weight, e$EC50, e$n)
  }
  
  return(grn)
}

#----All classes----
maxPrint <- function(charvec, maxprint = 10, len = NULL) {
	if(is.null(len)){
		len = length(charvec)
	}
	txt = paste0('(', len, ')')
	
	if (length(charvec) > maxprint) {
		charvec = c(charvec[1:maxprint],  '...')
	}
	txt = paste(txt, paste(charvec, collapse = ', '))
	return(txt)
}

