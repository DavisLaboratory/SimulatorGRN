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
grn=new('GraphGRN')
grn=addNodeRNA(grn, 'A')
grn=addNodeRNA(grn, 'B')
grn=addNodeRNA(grn, 'C')
grn=addEdgeReg(grn,'A','C')
grn=addEdgeReg(grn,'B','C')
sim = new('SimulationGRN', graph = grn, seed = 36)

sim$inputModels = list('A' = list('prop' = c(0.5, 0.5),
                                  'mean' = c(0.5, 0.6),
                                  'sd' = c(0.05, 0.05)),
                       'B' = list('prop' = 1,
                                  'mean' = 0.5,
                                  'sd' = 0.1))
sim$noiseL = 0.3

dm = simulateDataset(sim, 100)
dm = addNoise(sim, dm)
df = as.data.frame(t(dm))
df$cond = Mclust(df$A, verbose = F)$classification
ggplot(df, aes(B, C)) +geom_point() + facet_wrap(~cond)

set.seed(36)
dm = log(dm * exp(rnorm(3, 8, 2)) %*% t(rep(1, 100)))
df = as.data.frame(t(dm))
df$cond = Mclust(df$A, verbose = F)$classification
ggplot(df, aes(B, C)) +geom_point() + facet_wrap(~cond)


