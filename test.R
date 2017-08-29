source('GraphGRN-methods.R')
source('GraphGRN.R')
source('SimulationGRN-methods.R')
source('SimulationGRN.R')
library(nleqslv)
library(distr)
library(foreach)
library(ggplot2)
library(grid)
library(parallel)
library(doParallel)
library(foreach)
library(igraph)

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
nleqslv(exprs,ode,jac=NULL,ext,jacobian=T)

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

#----case 5----
df = read.csv('sourceNets/Yeast_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

grnEColi = df2GraphGRN(edges, loops = F, propand = 0.1, seed = 34234)
simEColi =new('SimulationGRN', graph = grnEColi, seed = 439591)
t = Sys.time()
emat = simulateDataset(simEColi, 500)
t = Sys.time() - t
print(t)

#plot graph
grnSmall = sampleGraph(grnEColi, 20, minregs = 10, seed = 173153)
dfl = GraphGRN2df(grnSmall)
dfl$edges = dfl$edges[, c(1, 3, 2, 4:6)]
g = graph_from_data_frame(d = dfl$edges, directed = T, vertices = dfl$nodes)
plot(g, vertex.size = 10, mark.group = dfl$nodes$name[dfl$nodes$type %in% 'and'])

simSmall =new('SimulationGRN', graph = grnSmall, seed = 528245)
t = Sys.time()
emat = simulateDataset(simSmall, 500)
t = Sys.time() - t

ma = as.data.frame(t(emat))
ggplot(ma, aes(nlpD_rpoS, aldB, colour = yhdG_fis)) + geom_point() + scale_color_distiller(palette = 'YlOrRd', direction = 1)
ggplot(ma, aes(yhdG_fis, aldB, colour = nlpD_rpoS)) + geom_point() + scale_color_distiller(palette = 'YlOrRd', direction = 1)
par(mfrow = c(1,2))
hist(emat['nlpD_rpoS',])
hist(emat['yhdG_fis',])


