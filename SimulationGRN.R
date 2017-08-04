#----SimulationGRN----
setClass(
	Class = 'SimulationGRN',
	slots = list(
		graph = 'GraphGRN',
		externalInputs = 'numeric',
		noiseL = 'numeric',
		noiseG = 'numeric',
		seed = 'numeric',
		solution = 'list'
	)
)

setValidity('SimulationGRN', validSimulationGRN)

setMethod(
	f = 'initialize',
	signature = 'SimulationGRN',
	definition = initSimulationGRN
)

setMethod(
	f = 'show',
	signature = 'SimulationGRN',
	definition = function(object) {
		mxp = 10
		
		cat('Graph:', length(object@graph@nodeset), 'nodes,', length(object@graph@edgeset), 'edges', '\n')
		cat('Local noise ratio:', object@noiseL, '\n')
		cat('Global noise ratio:', object@noiseG, '\n')
		cat('Randomization seed:', object@seed, '\n')
		
		#external inputs and solution print
		ext = round(object@externalInputs, digits = 2)
		soln = round(object$solution, digits = 2)
		if (length(ext) > mxp){
			ext = ext[1:mxp]
			ext = c(ext, '...' = paste0('...'))
		}
		if (length(soln) > mxp){
			soln = soln[1:mxp]
			soln = c(soln, '...' = '...')
		}
		
		cat(paste0('External Inputs (', length(object@externalInputs),' nodes):'), '\n')
		print(ext, quote = F)
		cat(paste0('Solution (', length(object$solution),' nodes):'), '\n')
		print(soln, quote = F)
	}
)

setMethod(
	f = '$',
	signature = 'SimulationGRN',
	definition = function(x, name) {
	  value = slot(x, name)
	  
	  if (name %in% 'solution'){
	    soln = x@solution
	    soln = c(soln$x, x@externalInputs)
	    soln = soln[order(names(soln))]
	    value = soln
	  }
	  
		return(value)
	}
)

setMethod(
	f = '$<-',
	signature = 'SimulationGRN',
	definition = function(x, name, value) {
		if(name %in% 'solution'){
			stop('Steady state solution cannot be modified directly')
		}
		
		slot(x, name)<-value
		if(name %in% 'externalInputs'){
			x@solution = solveSteadyState(x)
		}
		
		validObject(x)
		return(x)
	}
)



