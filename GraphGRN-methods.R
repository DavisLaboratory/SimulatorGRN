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
	
	#Outgoing edges
	if (!all(sapply(object@outedges, is, 'Edge'))) {
		stop('All outedges must be valid Edge objects')
	}
	
	#Incoming edges: check that all edges have to as this node
	testin = sapply(object@inedges, function(x) {
		return(identical(x$to$name, object$name))
	})
	
	if (!all(testin)) {
		stop('All interactions must have the \'to\' node as the current node')
	}
	
	#Outgoing edges: check that all edges have this node as part of from
	testout = sapply(object@outedges, function(x) {
		return(object$name %in% sapply(x$from, function (n) n$name))
	})
	
	if (!all(testout)) {
		stop('All interactions must have the current node as part of \'from\'')
	}
	
	#Incoming edges: check that all regulators are unique
	regnames = sapply(object@inedges, function(x) {
		sapply(x$from, slot, 'name')
	})
	regnames = as.vector(regnames)
	if (length(unique(regnames)) < length(regnames)) {
		stop('All regulators must be unique')
	}
	
	#outgoing edges: check that all targets are unique
	tgtnames = sapply(object@outedges, function(x) x$to$name)
	tgtnames = as.vector(tgtnames)
	if (length(unique(tgtnames)) < length(tgtnames)) {
		stop('All targets must be unique')
	}
	
	return(TRUE)
}

initNode <- function(.Object, ..., name = '', spmax = 1, spdeg = 1, inedges = list(), outedges = list()) {
	.Object@name = name
	.Object@spmax = spmax
	.Object@spdeg = spdeg
	.Object@inedges = inedges
	.Object@outedges = outedges
	
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

#----GraphGRN: core functions----
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

removeMissing <- function(graph) {
  #remove nodes with no interactions
  intnodes = sapply(graph@edgeset, function (e) sapply(e$from, function(n) n$name))
  intnodes = c(intnodes, sapply(graph@edgeset, function (e) e$to$name))
  intnodes = unique(intnodes)
  
  #warning and remove
  nonintnodes = setdiff(names(graph@nodeset), intnodes)
  if (length(nonintnodes) != 0) {
    msg = paste0('Nodes without interactions removed: ', paste(nonintnodes, collapse = ', '))
    warning(msg)
    
    graph@nodeset = graph@nodeset[intnodes]
  }
  return(graph)
}

#----GraphGRN: specific functions----
getODEFunc <- function(graph) {
  graph = removeMissing(graph)
  fn = 'function(exprs, externalInputs) {'
  
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

getEdgeClass <- function(edgetype) {
  if (edgetype %in% 'or')
    edgeclass = 'EdgeOr'
  else if (edgetype %in% 'and')
    edgeclass = 'EdgeAnd'
  else
    stop('Unrecognised edgetype')
  
  return(edgeclass)
}

orToAnd <- function(graph, from, to, weight) {
  #ensure or edges exist and are compatible to merge (i.e. same to node)
  if (length(to) > 1)
    stop('Multiple to nodes provided, expected 1')
  
  if (length(from) < 2)
    stop('Need 2 or more edges to merge into an AND edge')
  
  edgenames = paste(from, to, sep = '->')
  if (any(!edgenames %in% names(graph@edgeset))) {
    stop('Edges not found: ', paste(setdiff(edgenames, names(graph@edgeset)), collapse = ', '))
  }
  
  oredges = graph@edgeset[edgenames]
  #remove OR edges from graph
  for (e in oredges) {
    graph = removeEdge(graph, sapply(e$from, function (x) x$name), e$to$name)
  }
  
  #add AND edge
  graph = addEdge(
    graph,
    edgetype = 'and',
    from = from,
    to = to,
    activation = as.logical(sapply(oredges, function(e) e$activation)),
    weight = weight,
    EC50 = as.numeric(sapply(oredges, function(e) e$EC50)),
    n = as.numeric(sapply(oredges, function(e) e$n))
  )
  
  return(graph)
}

andToOr <- function(graph, from, to) {
  #ensure or edges exist and are compatible to merge (i.e. same to node)
  if (length(to) > 1)
    stop('Multiple to nodes provided, expected 1')

  if (length(from) < 2)
    stop('Need 2 or more edges to define an AND edge')

  andedge = getEdge(graph, from, to)
  if (is.null(andedge)) {
    stop('Edge not found')
  }

  #remove AND edge
  graph = graph = removeEdge(graph, sapply(andedge$from, function (x)
    x$name), andedge$to$name)

  #add OR edges
  for (i in 1:length(andedge$from)){
    graph = addEdge(
      graph,
      edgetype = 'or',
      from = andedge$from[[i]]$name,
      to = andedge$to$name,
      activation = andedge$activation[i],
      weight = andedge$weight,
      EC50 = andedge$EC50[i],
      n = andedge$n[i]
    )
  }

  return(graph)
}

subsetGraph <- function(graph, snodes) {
  #remove non-existing nodes from the graph
  for (n in setdiff(names(graph@nodeset), snodes)) {
    graph = removeNode(graph, n)
  }
  return(graph)
}

sampleSubNetwork <- function(graph, size, k, seed) {
  #get the adjacency matrix for the graph
  A = getAM(graph)
  
  #calculate total number of edges
  m = sum(diag(A)) + sum(A[upper.tri(A)])
  
  #calculate theoretical edge numbers between nodes
  degrees = rowSums(A)
  P = (degrees %*% t(degrees)) / (2 * m)
  B = A - P
  
  #start with empty network
  s = rep(-1, ncol(A)) # 1 or -1 for inc. or excl. resp.
  
  #set seed and start adding neighbours
  set.seed(seed)
  s[sample(which(s < 0), 1)] = 1
  
  for (i in 2:size) {
    #find neighbours
    neighbours = which(colSums(A[s > 0, , drop = F]) > 0 & s < 0)
    if (length(neighbours) == 0) {
      s[sample(which(s < 0), 1)] = 1
      next
    }
    
    #generate subnetworks
    subs = s %*% t(rep(1, length(s)))
    diag(subs) = 1
    subs = subs[ , neighbours]
    #calculate modulatiry
    Q = diag((t(subs) %*% B %*% subs) / (4 * m))
    cand = neighbours[Q >= quantile(Q, 1-k)] #candidates for addition
    s[cand[sample.int(length(cand),1)]] = 1
  }
  
  #return names of nodes in sampled subgraph
  return(colnames(A)[s > 0])
}

getAM <- function(graph, directed = F) {
  nodes = graph@nodeset
  edges = graph@edgeset
  
  A = matrix(rep(0, length(nodes) ^ 2), nrow = length(nodes))
  colnames(A) = rownames(A) = names(nodes)
  
  #generate adjacency matrix: true edge numbers between nodes
  for (e in edges) {
    from = sapply(e$from, function(x) x$name)
    to = e$to$name
    A[from, to] = 1
  }
  
  #if undirected
  if (!directed) {
    A = apply(A, 2, as.logical)
    A = A | t(A)
    A = apply(A, 2, as.numeric)
    rownames(A) = colnames(A)
  }
  
  return(A)
}

#----GraphGRN: Conversion functions----
dfToGraphGRN <- function(edges, nodes, propand = 0.3, loops = F, seed = sample.int(1E6, 1)) {
  if (missing(nodes) || is.null(nodes)) {
    nodes = data.frame('node' = unique(c(edges[ , 1], edges[ , 3])), stringsAsFactors = F)
  }
  
  #remove loops if required
  if (!loops) {
    edges = edges[edges[, 1] != edges[, 3], ]
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
  for (i in 1:nrow(nodes)) {
    n = nodes[i, , drop = F]
    grn = addNode(grn, n$node, n$tau, n$max, n$deg)
  }
  
  #add edges
  for (i in 1:nrow(edges)) {
    e = edges[i, , drop = F]
    grn = addEdge(grn, e$from, e$to, e$type, e$activation, e$weight, e$EC50, e$n)
  }
  
  #convert propand proportion of or's to and's
  #identify nodes with more than 2 inputs
  totaledges = length(grn@edgeset)
  andedges = round(totaledges * propand / 2)
  edgeset = grn@edgeset
  nodesin = sapply(nodeset, function(x) sum(sapply(x$inedges, is, 'EdgeOr')))
  
  #sample and convert to and edges
  # set.seed(seed)
  # for (i in 1:andedges) {
  #   candtgts = names(nodesin)[nodesin > 1]
  #   toNode = sample(candtgts, 1)
  #   
  #   #sample 2 OR edges to combine
  #   inedges = getNode(grn, toNode)$inedges
  #   inedges = inedges[sapply(inedges, is, 'EdgeOr')]
  #   fromNodes = sapply(sample(inedges, 2), function (x) x$from[[1]]$name)
  #   
  #   #convert OR to AND edge
  #   grn = mergeOr(grn, fromNodes, toNode, 1)
  #   nodesin[toNode] = nodesin[toNode] - 2
  # }
  
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

