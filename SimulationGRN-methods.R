#----Simulation: functions----
validSimulationGRN <- function(object) {
	#graph is a GraphGRN object
	if (!is(object@graph, 'GraphGRN')) {
		stop('graph must be a GraphGRN object')
	}
	
	#local noise
	if (object@noiseL<0 | object@noiseL>1) {
		stop('Local noise ratio must be between 0 and 1')
	}
	
	#global noise
	if (object@noiseG<0 | object@noiseG>1) {
		stop('Global noise ratio must be between 0 and 1')
	}
	
	return(TRUE)
}

initSimulationGRN <- function(.Object, ..., graph, noiseL = 0, noiseG = 0, seed = sample.int(1e6,1), inputModels = list()) {
	.Object@graph = graph
	.Object@noiseL = noiseL
	.Object@noiseG = noiseG
	.Object@seed = seed
	.Object@inputModels = inputModels
	
	if (length(inputModels) == 0) {
	  .Object = generateInputModels(.Object)
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

createInputModels <- function(simulation) {
  set.seed(simulation@seed)
  
  #create input models
  innodes = getInputNodes(simulation@graph)
  inmodels = list()
  
  for (n in innodes) {
    parms = list()
    mxs = sample(c(1, 2), 1, prob = c(0.7, 0.3))
    
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

simDataset <- function(simulation, numsamples, externalInputs) {
  set.seed(simulation@seed)
  
  #generate input matrix
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
    mix = sample(1:length(m$prop), prob = m$prop, replace = T)
    
    outbounds = 1
    while (sum(outbounds) > 0){
      outbounds = externalInputs[ , n] < 0 | externalInputs[ , n] > 1
      externalInputs[outbounds & mix == 1, n] = rnorm(sum(outbounds & mix == 1), m$mean[1], m$sd[1])
      if (length(m$prop) > 1) {
        externalInputs[outbounds & mix == 2, n] = rnorm(sum(outbounds & mix == 2), m$mean[2], m$sd[2])
      }
    }
  }
  
  #solve ODEs for different inputs
  emat = c()
  termcd = c()
  for (i in 1:numsamples) {
    soln = solveSteadyState(simulation, externalInputs[i, ])
    emat = cbind(emat, soln$x)
    termcd = c(termcd, soln$termcd)
  }
  emat = rbind(emat, t(externalInputs))
  colnames(emat) = paste0('sample', 1:numsamples)
  
  #check for errors
  if (!all(termcd == 1)) {
    nc = termcd != 1
    msg = 'Simulations for the following samples did not converge:'
    sampleids = paste(colnames(emat)[nc], ' (', termcd[nc], ')', sep = '')
    msg = paste(c(msg, sampleids), collapse = '\n\t')
    msg = paste(msg, 'format: sampleid (termination condition)', sep = '\n\n\t')
    warning(msg)
  }
  
  return(emat)
}



