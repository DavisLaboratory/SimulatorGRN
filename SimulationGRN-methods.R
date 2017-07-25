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

initSimulationGRN <- function(.Object, ..., graph, noiseL = 0, noiseG = 0) {
	.Object@graph = graph
	.Object@noiseL = noiseL
	.Object@noiseG = noiseG
	
	validObject(.Object)
	return(.Object)
}