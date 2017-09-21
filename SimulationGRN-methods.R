#----Simulation: functions----
validSimulationGRN <- function(object) {
	#graph is a GraphGRN object
	if (!is(object@graph, 'GraphGRN')) {
		stop('graph must be a GraphGRN object')
	}
	
	#local noise
	if (object@expnoise<0) {
		stop('Ecperimental noise standard deviation must be greater than 0')
	}
	
	#global noise
	if (object@bionoise<0) {
		stop('Biological noise standard deviation must be greater than 0')
	}
	
	return(TRUE)
}

initSimulationGRN <- function(.Object, ..., graph, expnoise = 0, bionoise = 0, seed = sample.int(1e6,1), inputModels = list(), propBimodal = 0) {
	.Object@graph = graph
	.Object@expnoise = expnoise
	.Object@bionoise = bionoise
	.Object@seed = seed
	.Object@inputModels = inputModels
	
	if (length(inputModels) == 0) {
	  .Object = generateInputModels(.Object, propBimodal)
	}
	
	validObject(.Object)
	return(.Object)
}

solveSteadyState <- function(object, externalInputs) {
  #external inputs
  if (is.null(names(externalInputs)) |
      !all(names(externalInputs) %in% getInputNodes(object@graph))) {
    stop('Invalid external inputs vector, named vector expected for ALL input nodes')
  }
  
  #set random seed
  set.seed(object@seed)
  
  #solve ODE
	ode = generateODE(object@graph)
	ext = externalInputs
	graph = object@graph
	nodes = setdiff(nodenames(graph), names(ext))
	exprs = rnorm(length(nodes), mean = 0.5, sd = 0.3 / 3) #3sd = 0.3 range
	exprs[exprs < 0] = 0
	exprs[exprs > 1] = 1
	names(exprs) = nodes
	
	soln = nleqslv(exprs, ode, jac = NULL, ext)
	
	#check if convergence is reached or not
	if(soln$termcd != 1) {
	  warning('Solution not achieved. use \'diagnostics(simulation)\' to get details')
	}
	return(soln)
}

createInputModels <- function(simulation, propBimodal) {
  set.seed(simulation@seed)
  
  #create input models
  innodes = getInputNodes(simulation@graph)
  inmodels = list()
  
  for (n in innodes) {
    parms = list()
    mxs = sample(c(1, 2), 1, prob = c(1 - propBimodal, propBimodal))
    
    if (mxs == 2) {
      parms = c(parms, 'prop' = runif(1, 0.2, 0.8))
      parms$prop = c(parms$prop, 1 - parms$prop)
    } else {
      parms$prop = 1
    }
    
    parms$mean = runif(mxs, 0.1, 0.9)
    maxsd = pmin(parms$mean, 1 - parms$mean) / 3
    parms$sd = sapply(maxsd, function(x) runif(1, 0.01, x))
    inmodels = c(inmodels, list(parms))
  }
  
  names(inmodels) = innodes
  simulation@inputModels = inmodels
  
  return(simulation)
}

generateInputData <- function(simulation, numsamples) {
  set.seed(simulation@seed)
  
  innodes = getInputNodes(simulation@graph)
  externalInputs = matrix(-1,nrow = numsamples, ncol = length(innodes))
  colnames(externalInputs) = innodes
  
  #create input models
  if (length(simulation@inputModels) == 0) {
    simulation = generateInputModels(simulation)
  }
  
  #simulate external inputs
  inmodels = simulation@inputModels
  for (n in innodes) {
    m = inmodels[[n]]
    mix = sample(1:length(m$prop), numsamples, prob = m$prop, replace = T)
    
    outbounds = 1
    while (sum(outbounds) > 0){
      outbounds = externalInputs[ , n] < 0 | externalInputs[ , n] > 1
      externalInputs[outbounds & mix == 1, n] = rnorm(sum(outbounds & mix == 1), m$mean[1], m$sd[1])
      if (length(m$prop) > 1) {
        externalInputs[outbounds & mix == 2, n] = rnorm(sum(outbounds & mix == 2), m$mean[2], m$sd[2])
      }
    }
  }
  
  colnames(externalInputs) = innodes
  return(externalInputs)
}

simDataset <- function(simulation, numsamples, externalInputs) {
  #generate input matrix
  innodes = getInputNodes(simulation@graph)
  if (!missing(externalInputs) && !is.null(externalInputs)) {
    if (nrow(externalInputs) != numsamples |
        length(setdiff(innodes, colnames(externalInputs))) != 0) {
          stop('Invalid externalInputs matrix provided')
    }
    externalInputs = externalInputs[, innodes]
  } else{
    externalInputs = generateInputData(simulation, numsamples)
  }
  
  #set random seed
  set.seed(simulation@seed)

  #solve ODE
  graph = simulation@graph
  ode = generateODE(graph)
  
  #generate LN noise for simulation
  lnnoise = exp(rnorm(numsamples * length(nodenames(graph)), 0, simulation@bionoise))
  lnnoise = matrix(lnnoise, nrow = numsamples, byrow = T)
  colnames(lnnoise) = nodenames(graph)
  
  #initialize solutions
  nodes = setdiff(nodenames(graph), colnames(externalInputs))
  exprs = rnorm(length(nodes) * numsamples, mean = 0.5, sd = 0.2 / 3) #3sd = 0.3 range
  exprs[exprs < 0] = 0
  exprs[exprs > 1] = 1
  exprs = matrix(exprs, nrow = numsamples)
  colnames(exprs) = nodes

  #solve ODEs for different inputs
  res = foreach(i = 1:numsamples, .packages = c('nleqslv'), .combine = cbind) %dopar% {
    soln = nleqslv(exprs[i, ], ode, externalInputs = externalInputs[i, ], lnnoise = lnnoise[i, ])
    return(c(soln$x, soln$termcd))
  }
  
  res = rbind(res)
  
  termcd = res[nrow(res),]
  emat = res[-(nrow(res)), , drop = F]
  emat = rbind(emat, t(externalInputs))
  colnames(emat) = paste0('sample_', 1:numsamples)

  #check for errors
  if (!all(termcd == 1)) {
    nc = termcd != 1
    msg = 'Simulations for the following samples did not converge:'
    sampleids = paste(colnames(emat)[nc], ' (', termcd[nc], ')', sep = '')
    msg = paste(c(msg, sampleids), collapse = '\n\t')
    msg = paste(msg, 'format: sampleid (termination condition)', sep = '\n\n\t')
    warning(msg)

    emat = emat[, !nc]
  }

  return(emat)
}

#netbenchmark strategy
addNoiseC <- function(simulation, simdata){
  genesds = apply(simdata, 1, sd)
  expnoisesds = runif(length(genesds), 0.8 * simulation@expnoise, 1.2 * simulation@expnoise)
  expnoisesds = expnoisesds * genesds
  
  bionoisesds = runif(length(genesds), 0.8 * simulation@bionoise, 1.2 * simulation@bionoise)
  bionoisesds = bionoisesds * mean(genesds)
  
  #generate noise matrices
  set.seed(simulation@seed)
  noisematL = c()
  noisematG = c()
  for (i in 1:length(genesds)) {
    noisematL = rbind(noisematL, rnorm(ncol(simdata), 0, expnoisesds[i]))
    noisematG = rbind(noisematG, rnorm(ncol(simdata), 0, bionoisesds[i]))
  }
  
  rownames(noisematL) = rownames(noisematG) = rownames(simdata)
  colnames(noisematL) = colnames(noisematG) = colnames(simdata)
  
  #add generated noise to data
  noisydata = simdata
  noisydata = log(noisydata + exp(noisematG))
  noisydata = noisydata + noisematL
  
  #ensure range of data is 0-1
  noisydata[noisydata > 1] = 1
  noisydata[noisydata < 0] = 0
  
  return(noisydata)
}

generateSensMat <- function(simulation, pertb, inputs = NULL, pertbNodes = NULL) {
  set.seed(simulation@seed)
  graph = simulation@graph
  
  if (is.null(inputs)) {
    inputs = runif(length(getInputNodes(graph)), pertb + 1E-4, 1)
    names(inputs) = getInputNodes(graph)
  }else if (!all(getInputNodes(graph) %in% names(inputs))) {
    stop('Missing Inputs')
  }
  
  if (is.null(pertbNodes)) {
    pertbNodes = nodenames(graph)
  } else{
    pertbNodes = intersect(pertbNodes, nodenames(graph))
  }
  
  #generate ODE functions
  odes = list()
  graph = simulation@graph
  ode = generateODE(graph)
  
  for (n in pertbNodes) {
    if (!n %in% names(inputs)) {
      getNode(graph, n)$spmax = getNode(graph, n)$spmax - pertb
      odes = c(odes, getODEFunc(graph))
      getNode(graph, n)$spmax = getNode(graph, n)$spmax + pertb
    }
  }
  names(odes) = setdiff(pertbNodes, names(inputs))
  
  #generate LN noise for simulation - 0 noise
  lnnoise = exp(rep(0, length(nodenames(graph))))
  names(lnnoise) = nodenames(graph)
  
  #initialize solutions
  nodes = setdiff(nodenames(graph), names(inputs))
  exprs = rnorm(length(nodes) * (1 + length(pertbNodes)), mean = 0.5, sd = 0.3 / 3) #3sd = 0.3 range
  exprs[exprs < 0] = 0
  exprs[exprs > 1] = 1
  exprs = matrix(exprs, nrow = (1 + length(pertbNodes)))
  colnames(exprs) = nodes
  
  #original solution with no perturbations
  soln0 = nleqslv(exprs[(1 + length(pertbNodes)), ], ode, externalInputs = inputs, lnnoise = lnnoise)$x
  soln0 = c(inputs, soln0)
  
  #solve ODEs for different perturbations
  res = foreach(i = 1:length(pertbNodes), .packages = c('nleqslv'), .combine = rbind) %dopar% {
    n = pertbNodes[i]
    tempinputs = inputs
    if (n %in% names(inputs)){
      odefn = ode
      tempinputs[n] = max(tempinputs[n] - pertb, 1E-4)
    } else{
      odefn = odes[[n]]
    }
    
    soln = nleqslv(exprs[i, ], odefn, externalInputs = tempinputs, lnnoise = lnnoise)$x
    soln = c(tempinputs, soln)
    return(soln)
  }
  res = rbind(res)
  
  sensmat = c()
  for (i in 1:length(pertbNodes)) {
    n = pertbNodes[i]
    pertbamt = pertb
    if (n %in% names(inputs) & inputs[n] == 1E-4) {
      pertbamt = inputs[n]
    }
    
    #calculate sensitivity
    sensmat = rbind(sensmat, (res[i, ] - soln0)/-pertbamt * getNode(graph, n)$spmax/soln0)
  }
  
  rownames(sensmat) = pertbNodes
  return(sensmat)
}



