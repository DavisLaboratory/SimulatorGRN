#----Node: functions----
validNode <- function(object) {
	#RNA maximum in range
	if (!is.na(object@rnamax) & (object@rnamax < 0 | object@rnamax > 1)) {
		stop('RNA maximum expression has to be between 0 and 1')
	}
	
	#RNA degradation in range
	if (!is.na(object@rnadeg) & (object@rnadeg < 0 | object@rnadeg > 1)) {
		stop('RNA degradation rate has to be between 0 and 1')
	}
	
	#Time constant in range
	if (!is.na(object@tau) & (object@tau <= 0)) {
		stop('Time constant must be positive')
	}
	
	#Incoming edges
	if (!all(sapply(object@inedges, is, 'Edge'))) {
		stop('All inedges must be valid Edge objects')
	}
	
	#Incoming edges: check that all edges have to as this node
	testin = sapply(object@inedges, function(x) {
		return(identical(x$to$name, object$name))
	})
	
	if (!all(testin)) {
		stop('All interactions must have the \'to\' node as the current node')
	}
	
	#check that all regulators are unique
	regnames = sapply(object@inedges, function(x) {
		sapply(x$from, slot, 'name')
	})
	if (length(unique(regnames)) < length(regnames)) {
		stop('All regulators must be unique')
	}
	
	return(TRUE)
}

initNode <- function(.Object, ..., name, rnamax = 1, rnadeg = 1, tau = 1, inedges = list()) {
	.Object@name = name
	.Object@rnamax = rnamax
	.Object@rnadeg = rnadeg
	.Object@tau = tau
	.Object@inedges = inedges
	
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
			combEqns = c(combEqns, apply(combs, 2, paste, collapse = ' * '))
		}
		
		act = paste(c(act, combEqns), collapse = ' - ')
	}
	
	#generate rate equation
	rateEqn = paste('(', act, ')', sep = '')
	rateEqn = paste(rateEqn, object@rnamax, sep = ' * ')
	degradationEqn = paste(object@rnadeg, object@name, sep = ' * ')
	rateEqn = paste(rateEqn, degradationEqn, sep = ' - ')
	return(rateEqn)
}

#----Edge: functions----
validEdge <- function(object) {
	#RNA maximum in range
	if (any(is.na(object@weight)) | sum(object@weight < 0 | object@weight > 1) > 0) {
		stop('Interaction weight has to be between 0 and 1')
	}
	
	#RNA degradation in range
	if (any(is.na(object@EC50)) | sum(object@EC50 < 0 | object@EC50 > 1) > 0) {
		stop('EC50 has to be between 0 and 1')
	}
	
	#Time constant in range
	if (any(is.na(object@n)) | sum(object@n == 1) > 0) {
		stop('Hill constant (n) cannot be 1')
	}
	
	#from are all of class Node
	if (!all(sapply(object@from, is, 'Node'))) {
		stop('All from nodes must be of class \'Node\'')
	}
	
	return(TRUE)
}

initEdge <- function(.Object, ..., from, to, weight = 1, EC50 = 0.5, n = 1.39, activation = T) {
	.Object@from = from
	.Object@to = to
	.Object@weight = weight
	.Object@EC50 = EC50
	.Object@n = n
	.Object@activation = activation
	
	validObject(.Object)
	return(.Object)
}

#----EdgeOr: functions----
validEdgeOr <- function(object) {
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
	numint = length(object@from) #number of interactors
	
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



