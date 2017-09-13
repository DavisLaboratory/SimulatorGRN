source('GraphGRN-methods.R')
source('GraphGRN.R')
source('SimulationGRN-methods.R')
source('SimulationGRN.R')
source('figures/analysis_methods.R')
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
library(EBcoexpress)

#----read network----
df = read.csv('sourceNets/EColi_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

#read in network and create simulator object
grnEColi = df2GraphGRN(edges, loops = F, propor = 0.1, seed = 36)
simEColi =new('SimulationGRN', graph = grnEColi, seed = 36, propBimodal = 0, noiseL = 0.1, noiseG = 0.1)

#----generate data----
#Total simulation params
nsamp = 100
simseed = 36
set.seed(simseed)

#modify input models
modinput = 'purR'
prop = runif(1, 0.2, 0.8)
prop = c(prop, 1 - prop)
mu = c(runif(1, 0.1, 0.5), runif(1, 0.5, 0.9))
maxsd = pmin(mu, 1 - mu) / 3
sdev = sapply(maxsd, function(x) runif(1, 0.01, x))

simEColi$inputModels[[modinput]] = list('prop' = prop,
                                 'mean' = mu,
                                 'sd' = sdev)
#generate data and add noise
datamat = simulateDataset(simEColi, nsamp)
noisydatamat = addNoise(simEColi, datamat)

#----plot networks----
plotInference <- function(graph, modulator, scores, topn = 100) {
  scoredf = cbind(expand.grid('from' = rownames(datamat), 'to' = rownames(datamat)),
                  'score' = as.numeric(scores))
  
  dflist = GraphGRN2df(graph)
  dflist$edges = merge(dflist$edges, scoredf, all = T, by = c('from', 'to'))
  dflist$edges[is.na(dflist$edges$type), 'type'] = 'inferred'
  
  net = graph_from_data_frame(d = dflist$edges, directed = T, vertices = dflist$nodes)
  V(net)$size = log(degree(net - E(net)[E(net)$type %in% 'inferred'], mode = 'out') + 3) * 2
  V(net)$label = NA
  V(net)$color = 1
  V(net)$color[V(net)$name %in% modulator] = 2
  E(net)$arrow.size = .05
  E(net)$color = 'gray80'
  E(net)$color[!is.na(E(net)$activation) & E(net)$activation] = 'darkgreen'
  E(net)$color[!is.na(E(net)$activation) & !E(net)$activation] = 'red'
  E(net)$width = 1.5
  
  goldnet = net - E(net)[E(net)$type %in% 'inferred']
  infnet = net - E(net)[rank(-E(net)$score) > topn]
  set.seed(360)
  l = layout_with_fr(goldnet)
  
  plot(goldnet, layout = l,
       mark.group = list(c(modulator, adjacent_vertices(goldnet, modulator)[[1]]$name)),
       main = 'Gold standard network')
  plot(infnet, layout = l,
       mark.group = list(c(modulator, adjacent_vertices(goldnet, modulator)[[1]]$name)),
       main = 'Inferred network')
}

#generate a random sensitivity matrix
sensmat = sensitivityAnalysis(simEColi)

#make inference
classf = Mclust(noisydatamat[modinput,], verbose = F)$classification
zsc = z.score(noisydatamat, classf)
zsc[is.infinite(zsc)] = 0
# ftgi = ftgi.score(noisydatamat, classf)
# ebsc = ebcoexpress.score(noisydatamat, classf)
# caisc = cai.score(noisydatamat, classf)
# magicsc = magic.score(noisydatamat, classf)
# dicersc = dicer.score(noisydatamat, classf)
# diffcoexsc = diffcoex.score(noisydatamat, classf)
# mindysc = mindy.score(noisydatamat, classf)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
plotInference(grnEColi, modinput, zsc)
layout(matrix(c(1), 1, 1, byrow = TRUE))

















#plot inference results
df = as.data.frame(t(noisydatamat))
df$cond = Mclust(df[,modinput], verbose = F)$classification
ggplot(df, aes(metR, glyA, colour = purR)) + geom_point() + facet_wrap(~cond)

expr1=noisydatamat[,classf==1]
expr2=noisydatamat[,classf==2]
r1=cor(t(expr1), method = 'spearman')
r2=cor(t(expr2), method = 'spearman')
plotdf = data.frame('r1' = as.numeric(r1), 'r2' = as.numeric(r2), 'z' = as.numeric(abs(zsc)))
plotdf = cbind(plotdf, expand.grid(rownames(datamat), rownames(datamat)))
ggplot(plotdf[plotdf$z > 3.5,], aes(r1, r2, colour = z)) + geom_point()

