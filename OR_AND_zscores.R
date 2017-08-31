source('GraphGRN-methods.R')
source('GraphGRN.R')
source('SimulationGRN-methods.R')
source('SimulationGRN.R')
source('figures/analysis_methods.R')
library(nleqslv)
library(foreach)
library(ggplot2)
library(grid)
library(parallel)
library(doParallel)
library(mclust)

nsamp = 100
ec50range = c(0.4, 0.6)
nrange = c(1.01, 1.7)
murange = c(0.2, 0.8)
sdrange = c(0.025, 0.3)
iter = 100
andScores = numeric(iter)
orScores = numeric(iter)

#2 activators ABC->OR, DEF->AND
set.seed = 36
seeds = sample.int(1E6, iter)

t = Sys.time()
foreach(i=1:iter, .combine = c) %do% {
  ec50 = runif(4, ec50range[1], ec50range[2])
  n = runif(4, nrange[1], nrange[2])
  
  g = new('GraphGRN')
  for (nd in LETTERS[1:6]) {
    g = addNode(g, nd)
  }
  g = addEdge(g, 'A', 'C', EC50 = ec50[1], n = n[1], activation = T)
  g = addEdge(g, 'B', 'C', EC50 = ec50[2], n = n[2], activation = T)
  g = addEdge(g, c('D', 'E'), 'F', EC50 = ec50[3:4], n = n[3:4], activation = c(T, T), edgetype = 'and')
  sim = new('SimulationGRN', graph = g, seed = seeds[i])
  
  #set B & E to bimodal distributions
  mus = runif(3, murange[1], murange[2])
  sds = runif(3, sdrange[1], sdrange[2])
  inputM = list('A' = list('prop' = 1, 'mean' = mus[1], 'sd' = sds[1]),
                   'B' = list('prop' = 2, 'mean' = mus[2:3], 'sd' = sds[2:3]),
                   'D' = list('prop' = 1, 'mean' = mus[1], 'sd' = sds[1]),
                   'E' = list('prop' = 2, 'mean' = mus[2:3], 'sd' = sds[2:3])
  )
  
  #simulate data
  sim$inputModels = inputM
  dmat = as.data.frame(t(simulateDataset(sim, nsamp)))
  classB = Mclust(dmat$B, G = 2, verbose = F)$classification
  classE = Mclust(dmat$E, G = 2, verbose = F)$classification
  scoresB = z.score(t(dmat), classB)
  scoresE = z.score(t(dmat), classE)
  orScores[i] = scoresB['C', 'A']
  andScores[i] = scoresE['F', 'D']
}
t = difftime(Sys.time(), t, units = 'sec')

