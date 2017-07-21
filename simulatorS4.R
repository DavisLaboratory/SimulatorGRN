source('simulatorS4-methods.R')

#----Node----
setClass(
	Class = 'Node',
	slots = list(
		name = 'character',
		tau = 'numeric',
		rnamax = 'numeric',
		rnadeg = 'numeric',
		inedges = 'list'
	),
	prototype = list(
		tau = 1,
		rnamax = 1,
		rnadeg = 1,
		inedges = list()
	)
)

setValidity('Node',validNode)

setMethod(
	f = 'initialize',
	signature = 'Node',
	definition = initNode
)

setMethod(
	f = '$',
	signature = 'Node',
	definition = function(x, name) {
		return(slot(x, name))
	}
)

setMethod(
	f = '$<-',
	signature = 'Node',
	definition = function(x, name, value) {
		print(value)
		slot(x, name)<-value
		validObject(x)
		return(x)
	}
)

setMethod(
	f = 'show',
	signature = 'Node',
	definition = function(object) {
		# str=''
		# str = paste(str, paste('Name: ', object@name, sep = ''), '\n')
		# str = paste(str, paste('Time const: ', object@tau, sep = ''), '\n')
		# str = paste(str, paste('RNA max: ', object@rnamax, sep = ''), '\n')
		# str = paste(str, paste('RNA degradation rate: ', object@rnadeg, sep = ''), '\n')
		# print(str)
	}
)

setGeneric(
	name = 'generateEqn',
	def = function(object){
		standardGeneric('generateEqn')
	}
)

setMethod(
	f = 'generateEqn',
	signature = 'Node',
	definition = generateRateEqn
)

#----Edge----
setClass(
	Class = 'Edge',
	slots = list(
		from = 'list',
		to = 'Node',
		weight = 'numeric',
		EC50 = 'numeric',
		n = 'numeric',
		activation = 'logical'
	),
	prototype = list(
		weight = 1,
		EC50 = 0.5,
		n = 1.39,
		activation = T
	)
)

setValidity('Edge',validEdge)

setMethod(
	f = 'initialize',
	signature = 'Edge',
	definition = initEdge
)

setMethod(
	f = '$',
	signature = 'Edge',
	definition = function(x, name) {
		return(slot(x, name))
	}
)

setMethod(
	f = '$<-',
	signature = 'Edge',
	definition = function(x, name, value) {
		slot(x, name)<-value
		validObject(x)
		return(x)
	}
)

setMethod(
	f = 'show',
	signature = 'Edge',
	definition = function(object) {
		return('')
	}
)

setGeneric(
	name = 'generateActivationEqn',
	def = function(object){
		standardGeneric('generateActivationEqn')
	}
)

#----EdgeOr----
setClass(
	Class = 'EdgeOr',
	contains = 'Edge'
)

setValidity('EdgeOr',validEdgeOr)

setMethod(
	f = 'generateActivationEqn',
	signature = 'EdgeOr',
	definition = generateActivationEqnOr
)

#----EdgeAnd----
setClass(
	Class = 'EdgeAnd',
	contains = 'Edge'
)

setValidity('EdgeAnd',validEdgeAnd)

setMethod(
	f = 'initialize',
	signature = 'EdgeAnd',
	definition = initEdgeAnd
)

setMethod(
	f = 'generateActivationEqn',
	signature = 'EdgeAnd',
	definition = generateActivationEqnAnd
)
