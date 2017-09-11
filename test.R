source('GraphGRN-methods.R')
source('GraphGRN.R')
source('SimulationGRN-methods.R')
source('SimulationGRN.R')
library(nleqslv)
library(distr)
library(ggplot2)
library(grid)
library(parallel)
library(doParallel)
library(foreach)
library(igraph)
library(stringr)
library(MASS)
library(mclust)

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
df = read.csv('sourceNets/EColi_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

#read in network and create a sampled network
grnEColi = df2GraphGRN(edges, loops = F, propor = 0.1, seed = 36)
grnSmall = sampleGraph(grnEColi, 20, minregs = 10, seed = 36)

#create simulators
simEColi =new('SimulationGRN', graph = grnEColi, seed = 36, propBimodal = 0)
simSmall =new('SimulationGRN', graph = grnEColi, seed = 36, propBimodal = 0)

#Total simulation params
nsamp = 100
simseed = 36
set.seed(simseed)

#modify input models
prop = runif(1, 0.2, 0.8)
prop = c(prop, 1 - prop)
mu = c(runif(1, 0.1, 0.5), runif(1, 0.5, 0.9))
maxsd = pmin(mu, 1 - mu) / 3
sdev = sapply(maxsd, function(x) runif(1, 0.01, x))

simEColi$inputModels$purR = list('prop' = prop,
                    'mean' = mu,
                    'sd' = sdev)
datamat = simulateDataset(simEColi, nsamp)


