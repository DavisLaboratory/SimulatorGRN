source('GraphGRN-methods.R')
source('GraphGRN.R')
source('SimulationGRN-methods.R')
source('SimulationGRN.R')
library(nleqslv)
library(distr)

#----case 1----
n1=new('Node', name='A')
n2=new('Node', name='B')
n3=new('Node', name='C')
n4=new('Node', name='D')
e1=new('EdgeOr', from=c(n2), to=n1, activation=F)
e2=new('EdgeOr', from=c(n3), to=n1)
e3=new('EdgeOr', from=c(n4), to=n1)
e4=new('EdgeOr', from=c(n3), to=n2)
and1=new('EdgeAnd', from=c(n1,n2,n3), to=n4)

n1$inedges=list(e1,e2,e3)
n2$inedges=list(e4)
n4$inedges=list(and1)

print(generateRateEqn(n1))
print(generateRateEqn(n2))
print(generateRateEqn(n3))
print(generateRateEqn(n4))

#----case 2----
grn=new('GraphGRN')
grn=addNode(grn, 'A')
grn=addNode(grn, 'B')
grn=addNode(grn, 'C')
grn=addNode(grn, 'D')
grn=addEdge(grn,'B','A',activation=F)
grn=addEdge(grn,'C','A')
grn=addEdge(grn,'D','A')
grn=addEdge(grn,'C','B')
grn=addEdge(grn,c('A','B','C'),'D',edgetype='and')

print(generateRateEqn(grn@nodeset[['A']]))
print(generateRateEqn(grn@nodeset[['B']]))
print(generateRateEqn(grn@nodeset[['C']]))
print(generateRateEqn(grn@nodeset[['D']]))

#----case 3----
grn2=new('GraphGRN')
for (i in LETTERS[1:12]){
	grn2=addNode(grn2,i)
}

grn2=addEdge(grn2,'A','L',activation=T)
grn2=addEdge(grn2,'A','D',activation=T)
grn2=addEdge(grn2,'D','G',activation=T)
grn2=addEdge(grn2,'B','E',activation=T)
grn2=addEdge(grn2,'E','I',activation=T)
grn2=addEdge(grn2,'C','F',activation=T)
grn2=addEdge(grn2,'F','K',activation=T)
grn2=addEdge(grn2,c('D','E'),'H',edgetype='and')
grn2=addEdge(grn2,c('E','F'),'J',edgetype='and',activation=c(T,F))
ode=generateODE(grn2)
exprs=numeric(9)
names(exprs)=LETTERS[4:12]
ext=c(0.6,0.6,0.6)
names(ext)=LETTERS[1:3]
nleqslv(exprs,ode,jac=NULL,grn2,ext,jacobian=T)

#----case 4----
sim1=new('SimulationGRN', graph = grn2, externalInputs = c('A' = 0.5, 'B' = 0.7, 'C' = 0.5))
mix = UnivarMixingDistribution(Norm(0.3, 0.2/3), Norm(0.7, 0.2/3), mixCoeff = c(0.5, 0.5))
rmix = r(mix)

samp = 500
ext = data.frame('A' = rmix(samp), 'B' = rnorm(samp, 0.7, 0.2/3), 'C' = rmix(samp))
ext2 = data.frame('A' = rnorm(samp, 0.5, 0.2/3), 'B' = rnorm(samp, 0.5, 0.2/3), 'C' =rmix(samp))
ma = foreach(i=1:samp, .combine = rbind, .packages = c('nleqslv')) %do% {
	e = as.numeric(ext[i, ])
	names(e) = colnames(ext)
	sim1$externalInputs = e
	return(sim1$solution)
}

hist(ma, breaks = 100)
ma = as.data.frame(ma)
p1 = ggplot(ma, aes(E, J, colour = F)) + geom_point() + scale_color_distiller(palette = 'YlOrRd', direction = 1)
p2 = ggplot(ma, aes(E, H, colour = D)) + geom_point() + scale_color_distiller(palette = 'YlOrRd', direction = 1)
multiplot(p1, p2, cols = 2)



