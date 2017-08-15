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

  return(TRUE)
}

initNode <- function(.Object, ..., name = '', spmax = 1, spdeg = 1, inedges = character(), outedges = character()) {
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

generateRateEqn <- function(node, graph) {
	inedges = node$inedges
	#no rate equations for input nodes
	if (length(inedges) == 0) {
		return('')
	}
	
	#generate activation functions for each interaction
	actEqns = sapply(graph@edgeset[inedges], generateActivationEqn)
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
	rateEqn = paste(rateEqn, node@spmax, sep = ' * ')
	degradationEqn = paste(node@spdeg, node@name, sep = ' * ')
	rateEqn = paste(rateEqn, degradationEqn, sep = ' - ')
	return(rateEqn)
}

#----Edge: functions----
validEdge <- function(object) {
	#weight in range
	if (any(is.na(object@weight)) | sum(object@weight < 0 | object@weight > 1) > 0) {
		stop('Interaction weight has to be between 0 and 1')
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
  name = paste(sort(from), collapse = '')
  name = paste(name, to, sep = '->')
  .Object@name = name
  
  validObject(.Object)
  return(.Object)
}

generateActivationEqnOr <- function(object) {
	e = object
	#generate activation eqn
	act = paste('fAct(', e$from, ',', e$EC50, ',', e$n, ')', sep =	'')
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
	
	#Activation length matches number of interactors
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
	name = paste(sort(from), collapse = '')
	name = paste(name, to, sep = '->')
	.Object@name = name
	
	validObject(.Object)
	return(.Object)
}

generateActivationEqnAnd <- function(object) {
	e = object
	
	#generate AND activation eqn
	act = character(length(e@from))
	for (i in 1:length(act)) {
		act[i] = paste('fAct(', e@from[i], ',', e@EC50[i], ',', e@n[i], ')', sep =	'')
		if (!e$activation[i]) {
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
	
	#edgeset are all of class Edge
	if (!all(sapply(object@edgeset, is, 'Edge'))) {
		stop('All nodes must be of class \'Edge\'')
	}
	
	#check names of edgeset
	if (!all(sapply(object@edgeset, function(x) x$name) == names(object@edgeset))){
		stop('Invalid graph generated. Use the \'addEdge\' method to add an edge to the graph.')
	}
  
  #edge checks
  nnames = names(object@nodeset)
  for (e in object@edgeset) {
    if (!all(e$from %in% nnames)) {
      stop('Some/all from nodes do not exist in edge: ', e$name)
    }
    if (!e$to %in% nnames) {
      stop('To node does not exist in edge: ', e$name)
    }
  }
  
  #node checks
  #check node duplication
  if (length(unique(nnames)) != length(nnames)) {
    stop('All nodes in the graph must be unique')
  }
  
  enames = names(object@edgeset)
  for (n in object@nodeset) {
    #check whether the incoming and outgoing edges exist
    if (!all(n$inedges %in% enames)) {
      stop('Some/all inbound edges do not exist in node: ', n$name)
    }
    if (!all(n$outedges %in% enames)) {
      stop('Some/all outbound edges do not exist in node: ', n$name)
    }
    
    #Incoming edges: check that all edges have to as this node
    infroms = unlist(sapply(object@edgeset[n$inedges], function(e) e$from))
    intos = sapply(object@edgeset[n$inedges], function(e) e$to)
    
    if (!all(intos %in% n$name)) {
      stop('All interactions must have the \'to\' node as current node for node: ', n$name)
    }
    
    #Incoming edges: check that all regulators are unique
    if (length(unique(infroms)) != length(infroms)) {
      stop('All regulators must be unique for node: ', n$name)
    }
    
    #Outgoing edges: check that all edges have from as this node
    outfroms = lapply(object@edgeset[n$outedges], function(e) e$from)
    outtos = sapply(object@edgeset[n$outedges], function(e) e$to)
    
    if (!all(sapply(outfroms, function (flist) any(flist %in% n$name)))) {
      stop('All interactions must have the current node as part of \'from\' for node: ', n$name)
    }
    
    #Outgoing edges: check that all targets are unique
    if (length(unique(outtos)) != length(outtos)) {
      stop('All targets must be unique for node: ', n$name)
    }
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
  intnodes = unlist(sapply(graph@edgeset, function (e) e$from))
  intnodes = c(intnodes, sapply(graph@edgeset, function (e) e$to))
  intnodes = unique(intnodes)
  
  #warning and remove
  nonintnodes = setdiff(names(graph@nodeset), intnodes)
  if (length(nonintnodes) != 0) {
    msg = paste0('Nodes without interactions removed: ', paste(nonintnodes, collapse = ', '))
    warning(msg)
    
    graph = removeNode(graph, nonintnodes)
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
    eqn = paste('\t\t', 'rates[\"', node$name, '\"] = ', generateRateEqn(node, graph), sep = '')
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

getEdgeClass <- function(edgetype) {
  if (edgetype %in% 'or')
    edgeclass = 'EdgeOr'
  else if (edgetype %in% 'and')
    edgeclass = 'EdgeAnd'
  else
    stop('Unrecognised edgetype')
  
  return(edgeclass)
}

rmnode <- function(graph, nodename) {
  #Check nodename exists
  if (length(nodename) > 1)
    stop('Multiple nodes provided, expected 1')
  
  if (!nodename %in% names(graph@nodeset)) {
    stop('Node not found')
  }
  
  #retrieve node
  node = getNode(graph, nodename)
  
  #remove from incoming interactions
  edges = graph@edgeset[node$inedges]
  for (e in edges) {
    graph = removeEdge(graph, e$from, e$to) #remove edge
  }
  
  #remove from inedges of targets
  edges = graph@edgeset[node$outedges]
  for (e in edges) {
    from = e$from
    graph = removeEdge(graph, from, e$to) #remove edge
    
    if (is(e, 'EdgeAnd')) {
      newfrom = !from %in% node$name
      
      if (sum(newfrom) > 1) {
        edgetype = 'and'
      } else{
        edgetype = 'or'
      }
      
      #create a new edge without the current node
      graph = addEdge(
        graph,
        edgetype = edgetype,
        from = e$from[newfrom],
        to = e$to,
        activation = e$activation[newfrom],
        weight = e$weight,
        EC50 = e$EC50[newfrom],
        n = e$n[newfrom]
      )
    }
  }
  
  #remove node from nodeset
  graph@nodeset = graph@nodeset[!names(graph@nodeset) %in% node$name]
  
  return(graph)
}

orToAnd <- function(graph, from, to, weight) {
  #ensure or edges exist
  if (length(to) > 1)
    stop('Multiple to nodes provided, expected 1')
  
  if (length(from) < 2)
    stop('Need 2 or more edges to merge into an AND edge')
  
  edgenames = paste(sort(from), to, sep = '->')
  if (any(!edgenames %in% names(graph@edgeset))) {
    stop('Edges not found: ', paste(setdiff(edgenames, names(graph@edgeset)), collapse = ', '))
  }
  
  oredges = graph@edgeset[edgenames]
  #remove OR edges from graph
  for (e in oredges) {
    graph = removeEdge(graph, e$from, e$to)
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
  #ensure or edges exist
  if (length(to) > 1)
    stop('Multiple to nodes provided, expected 1')

  if (length(from) < 2)
    stop('Need 2 or more edges to define an AND edge')

  andedge = getEdge(graph, from, to)
  if (is.null(andedge)) {
    stop('Edge not found')
  }

  #remove AND edge
  graph = graph = removeEdge(graph, andedge$from, andedge$to)

  #add OR edges
  for (i in 1:length(andedge$from)){
    graph = addEdge(
      graph,
      edgetype = 'or',
      from = andedge$from[i],
      to = andedge$to,
      activation = andedge$activation[i],
      weight = andedge$weight,
      EC50 = andedge$EC50[i],
      n = andedge$n[i]
    )
  }

  return(graph)
}

subsetGraph <- function(graph, snodes) {
  #remove nodes from the graph
  for (n in setdiff(names(graph@nodeset), snodes)) {
    graph = removeNode(graph, n)
  }
  return(graph)
}

sampleSubNetwork <- function(graph, size, minregs, k, seed) {
  #identify regulators
  regs = rowSums(getAM(graph))
  regs = names(regs)[regs > 0]
  
  #get the adjacency matrix for the graph
  A = getAM(graph, directed = F)
  hdegree = rowSums(A)
  hdegree = names(hdegree)[hdegree > 1]
  
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
    
    #sample the minimum number of regulators required
    if (minregs > 0) {
      newn = intersect(names(neighbours), regs)
      newnhdeg = intersect(names(neighbours), hdegree)
      
      if (length(newn) != 0) {
        neighbours = neighbours[newn]
        minregs = minregs - 1
      } else if (length(newnhdeg) != 0) {
        #favour higher degree neighbours to improve chances of hitting a reg
        neighbours = neighbours[newnhdeg]
      }
    }
    
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

getAMC <- function(graph, directed = T) {
  nodes = graph@nodeset
  edges = graph@edgeset
  
  A = matrix(rep(0, length(nodes) ^ 2), nrow = length(nodes))
  colnames(A) = rownames(A) = names(nodes)
  
  #generate adjacency matrix: true edge numbers between nodes
  for (e in edges) {
    A[e$from, e$to] = 1
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
df2GraphGRN <- function(edges, nodes, propand = 0.3, loops = F, seed = sample.int(1E6, 1)) {
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
  edgeset = grn@edgeset
  totaledges = length(edgeset)
  nodesin = sapply(grn@nodeset, function(x) {
    sum(as.numeric(sapply(edgeset[x$inedges], is, 'EdgeOr')))
  })
  
  #sample and convert to and edges
  set.seed(seed)
  andsize = c(2, 3, 4)
  andprobs = c(1, 0, 0)
  pAnd = 0
  nfroms = 2
  while (pAnd < propand) {
    candtgts = names(nodesin)[nodesin >= nfroms]
    
    if (length(candtgts) == 0) {
      if (nfroms == 2) {
        msg = paste0('Only ', round(pAnd, digits = 2), '/', propand,
                     ' AND edges could be created')
        warning(msg)
        break
      } else{
        nfroms = nfroms - 1
        #change probalilities
        andprobs[andsize > nfroms] = 0
        andprobs = andprobs / sum(andprobs)
        next
      }
    }
    toNode = sample(candtgts, 1)

    #sample nfrom(2) OR edges to combine
    inedges = edgeset[getNode(grn, toNode)$inedges]
    inedges = inedges[sapply(inedges, is, 'EdgeOr')]
    fromNodes = sapply(sample(inedges, nfroms), function (x) x$from)

    #convert OR to AND edge
    grn = mergeOr(grn, fromNodes, toNode, 1)
    nodesin[toNode] = nodesin[toNode] - nfroms
    
    #calculate proportions
    pAnd = sum(sapply(grn@edgeset, is, 'EdgeAnd')) / length(grn@edgeset)
    nfroms = sample(andsize, 1, prob = andprobs)
  }
  
  return(grn)
}

GraphGRN2df <- function(graph) {
  nodes = graph@nodeset
  edges = graph@edgeset
  
  #create node df
  nodedf = data.frame('name' = names(nodes), stringsAsFactors = F)
  nodedf$rnamax = sapply(nodes, slot, 'spmax')
  nodedf$rnadeg = sapply(nodes, slot, 'spdeg')
  nodedf$tau = sapply(nodes, slot, 'tau')
  nodedf$type = 'or'
  
  #create edge df
  #convert oredges
  oredges = edges[sapply(edges, is, 'EdgeOr')]
  edgedf = data.frame('from' = sapply(oredges, slot, 'from'), stringsAsFactors = F)
  edgedf$type = sapply(oredges, slot, 'activation')
  edgedf$to = sapply(oredges, slot, 'to')
  edgedf$weight = sapply(oredges, slot, 'weight')
  edgedf$EC50 = sapply(oredges, slot, 'EC50')
  edgedf$n = sapply(oredges, slot, 'n')
  
  #convert andedges
  andedges = edges[sapply(edges, is, 'EdgeAnd')]
  andedgem = c()
  andnodem = c()
  for (e in andedges) {
    es = c()
    newnode = paste(sort(e$from), collapse = '')
    andnodem = c(andnodem, newnode)
    
    #from to intermediate
    es = cbind(e$from, e$activation, newnode, NA, e$EC50, e$n)
    
    #intermediate
    es = rbind(es, c(newnode, T, e$to, e$weight, NA, NA))
    andedgem = rbind(andedgem, es)
  }
  
  colnames(andedgem) = colnames(edgedf)
  edgedf = rbind(edgedf, andedgem)
  
  andnodem = cbind(andnodem, NA, NA, NA, 'and')
  colnames(andnodem) = colnames(nodedf)
  nodedf = rbind(nodedf, andnodem)
  
  #convert types
  for (i in 2:4){
    nodedf[ , i] = as.numeric(nodedf[ , i])
  }
  
  for (i in 4:6){
    edgedf[ , i] = as.numeric(edgedf[ , i])
  }
  
  #create list of results
  rownames(nodedf) = NULL
  rownames(edgedf) = NULL
  dflist = list('nodes' = nodedf, 'edges' = edgedf)
  
  return(dflist)
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

