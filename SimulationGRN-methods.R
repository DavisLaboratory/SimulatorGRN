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
	
	#external inputs
	if (is.null(names(object@externalInputs)) |
		!all(names(object@externalInputs) %in% getInputNodes(object@graph))) {
		stop('Invalid external inputs vector, named vector expected')
	}
	
	return(TRUE)
}

initSimulationGRN <- function(.Object, ..., graph, externalInputs, noiseL = 0, noiseG = 0, seed = sample.int(1e12,1)) {
	if(missing(externalInputs)){
		inputNodes = getInputNodes(graph)
		externalInputs = numeric(length(inputNodes)) + 0.5
		names(externalInputs) = inputNodes
	}
	
	.Object@graph = graph
	.Object@noiseL = noiseL
	.Object@noiseG = noiseG
	.Object@seed = seed
	.Object@externalInputs = externalInputs
	.Object@solution = solveSteadyState(.Object)
	
	validObject(.Object)
	return(.Object)
}

solveSteadyState <- function(object) {
	ode = generateODE(object@graph)
	ext = object@externalInputs
	graph = object@graph
	nodes = setdiff(nodenames(graph), names(ext))
	exprs = runif(length(nodes))
	names(exprs) = nodes
	
	soln = nleqslv(exprs, ode, jac = NULL, ext, graph)
	return(soln)
}

