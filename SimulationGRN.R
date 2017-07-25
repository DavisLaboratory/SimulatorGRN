#----SimulationGRN----
setClass(
	Class = 'SimulationGRN',
	slots = list(
		graph = 'GraphGRN',
		noiseL = 'numeric',
		noiseG = 'numeric'
	)
)

setValidity(validSimulationGRN)

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
		slot(x, name)<-value
		validObject(x)
		return(x)
	}
)