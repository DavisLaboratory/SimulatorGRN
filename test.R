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
library(stringr)

#----case 2----
grn=new('GraphGRN')
grn=addNodeRNA(grn, 'A')
grn=addNodeRNA(grn, 'B')
grn=addNodeRNA(grn, 'C')
grn=addNodeRNA(grn, 'D')
grn=addEdgeReg(grn,'B','A',activation=F)
grn=addEdgeReg(grn,'C','A')
grn=addEdgeReg(grn,'D','A')
grn=addEdgeReg(grn,'C','B')
grn=addEdgeReg(grn,'A','D')
grn=addEdgeReg(grn,'B','D')
grn=addEdgeReg(grn,'C','D')
getNode(grn, 'D')$logiceqn = 'A & B & C'

print(generateRateEqn(grn@nodeset[['A']], grn))
print(generateRateEqn(grn@nodeset[['B']], grn))
print(generateRateEqn(grn@nodeset[['C']], grn))
print(generateRateEqn(grn@nodeset[['D']], grn))

#----case 3----
grn2=new('GraphGRN')
for (i in LETTERS[1:12]){
	grn2=addNodeRNA(grn2,i)
}

grn2=addEdgeReg(grn2,'A','L',activation=T)
grn2=addEdgeReg(grn2,'A','D',activation=T)
grn2=addEdgeReg(grn2,'D','G',activation=T)
grn2=addEdgeReg(grn2,'B','E',activation=T)
grn2=addEdgeReg(grn2,'E','I',activation=T)
grn2=addEdgeReg(grn2,'C','F',activation=T)
grn2=addEdgeReg(grn2,'F','K',activation=T)
grn2=addEdgeReg(grn2,'D','H',activation=T)
grn2=addEdgeReg(grn2,'E','H',activation=T)
grn2=addEdgeReg(grn2,'E','J',activation=T)
grn2=addEdgeReg(grn2,'F','J',activation=F)

ode=generateODE(grn2)
exprs=numeric(9)
names(exprs)=LETTERS[4:12]
ext=c(0.6,0.6,0.6)
names(ext)=LETTERS[1:3]
nleqslv(exprs,ode,jac=NULL,ext,jacobian=T)

#----case 4----
sim1=new('SimulationGRN', graph = grn2, seed = 360)
samp = 500
ma = simulateDataset(sim1, samp)

hist(ma, breaks = 100)
ma = as.data.frame(t(ma))
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


