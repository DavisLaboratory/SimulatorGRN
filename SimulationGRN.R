#----SimulationGRN----
setClass(
	Class = 'SimulationGRN',
	slots = list(
		graph = 'GraphGRN',
		externalInputs = 'numeric',
		noiseL = 'numeric',
		noiseG = 'numeric',
		seed = 'numeric',
		solution = 'numeric'
	)
)

setValidity('SimulationGRN', validSimulationGRN)

setMethod(
	f = 'initialize',
	signature = 'SimulationGRN',
	definition = initSimulationGRN
)

setMethod(
	f = '$',
	signature = 'SimulationGRN',
	definition = function(x, name) {
		return(slot(x, name))
	}
)

setMethod(
	f = '$<-',
	signature = 'SimulationGRN',
	definition = function(x, name, value) {
		if(name %in% 'solution'){
			stop('Steady state solution cannot be modified')
		}
		
		slot(x, name)<-value
		if(name %in% 'externalInputs'){
			x@solution = solveSteadyState(x)
		}
		
		validObject(x)
		return(x)
	}
)



