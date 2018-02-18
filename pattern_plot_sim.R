setwd('/wehisan/home/allstaff/b/bhuva.d/R_projects/cGRNsimulator')
source('/wehisan/home/allstaff/b/bhuva.d/R_projects/scripts/multiplot.R')
source('Classes.R')
source('GraphGRN_core.R')
source('GraphGRN.R')
source('SimulationGRN_core.R')
source('SimulationGRN.R')
source('analysis_methods.R')
library(nleqslv)
library(distr)
library(ggplot2)
library(grid)
library(parallel)
library(doParallel)
library(doSNOW)
library(foreach)
library(igraph)
library(stringr)
library(MASS)
library(mclust)
library(EBcoexpress)
library(diggit)
library(plyr)
library(ggsci)
library(pROC)
library(reshape2)
library(infotheo)
library(plotly)
library(COSINE)
library(fastLiquidAssociation)

originseed = 60897

#----read network----
org = 'EColi'
net_srcfile = file.path('sourceNets', paste0(org, '_full.sif'))
net_rdfile = file.path('sourceNets', paste0('grn', org, '_', originseed, '.RData'))
if (!file.exists(net_rdfile)){
  df = read.csv(net_srcfile, sep = '\t', header = F, stringsAsFactors = F)
  edges = df
  maplogical = c('ac' = T, 're' = F, 'du' = F)
  edges[,2] = maplogical[edges[,2]]
  
  #read in network and randomize params
  grnEColi = df2GraphGRN(edges, loops = F, propor = 0, seed = originseed)
  save(grnEColi, edges, file = net_rdfile)
} else {
  load(net_rdfile)
}
grnEColi = randomizeParams(grnEColi, 'linear-like', seed = originseed)

ncores = 10#detectCores()
set.seed(originseed)
nsamp = 500

#simulation specific parameters
expnoise = 0.05
bionoise = 0.05
mu = c(0.1, 0.5)
sdev = c(0.03, 0.15)
prop = c(0.5, 0.5)
netSize = 80
minTFs = 10

#----simulate data----
originseed = 60897
grnSmall = sampleGraph(grnEColi, netSize, minTFs, seed = originseed)
grnSmall = randomizeParams(grnSmall, 'linear-like', originseed)
simSmall = new('SimulationGRN', graph = grnSmall, seed = originseed, propBimodal = 0.3, expnoise = expnoise, bionoise = bionoise)

# #select a modulator
# iSmall = GraphGRN2igraph(grnSmall)
# modulatedNodes = degree(iSmall, mode = 'in')
# modulatedNodes = names(modulatedNodes)[modulatedNodes>=2]
# moddists = distances(iSmall)
# moddists[is.infinite(moddists)] = NA
# modinput = rowSums(moddists[,modulatedNodes]!=0, na.rm = T)
# modinput = names(modinput)[modinput != 0]
# modinput = intersect(getInputNodes(grnSmall), modinput)
# set.seed(originseed)
# bimodal = sample(modinput, 2)
# # bimodal = c('himA', 'fnr')
# 
# for (b in bimodal){
#   simSmall$inputModels[[b]] = list('prop' = prop,
#                                    'mean' = mu,
#                                    'sd' = sdev)
# }

#natural modulator
bimodal = unlist(lapply(simSmall$inputModels, function(x) length(x$prop)))
bimodal = names(bimodal)[bimodal == 2]

cl = makeCluster(ncores, outfile = '')
clusterExport(cl, c('addBioNoise'))
registerDoSNOW(cl)
datamat = simulateDataset(simSmall, nsamp)
stopCluster(cl)

#----generate a truth adjacency matrix----
inputs = sapply(simSmall$inputModels, function(x) x$mean[1])
inputs[bimodal] = 0.5
names(inputs) = getInputNodes(grnSmall)
cl = makeCluster(ncores, outfile = '')
clusterExport(cl, c('addBioNoise'))
registerDoSNOW(cl)
sensmat1 = sensitivityAnalysis(simSmall, 0.25, inputs, nodenames(grnSmall))
stopCluster(cl)

triplets = getGoldStandard(simSmall, threshold = 0.7, assocnet = T, sensmat = sensmat1)

#----determine direct correlations----


#----fit models and store parameters----
fitModels <- function(x, y){
  #rows are Y's and cols are X's in y=a+bx
  xbar = colMeans(x)
  ybar = colMeans(y)
  sxx = colSums((x-matrix(rep(1,nrow(x)))%*%xbar)^2)
  sxy = t(t(x-matrix(rep(1,nrow(x)))%*%xbar)%*%(y-matrix(rep(1, nrow(y)))%*%ybar))
  b = sxy/matrix(rep(sxx,ncol(y)),ncol = ncol(x), byrow=T)
  a = matrix(rep(ybar, ncol(x)), ncol = ncol(x)) - b*matrix(rep(xbar, ncol(y)), ncol = ncol(x), byrow = T)
  
  return(cbind(expand.grid('y' = rownames(b), 'x' = colnames(b)),
               'm' = as.numeric(b),
               'c' = as.numeric(a),
               'p' = as.numeric(cor(y, x))))
}

#fit linear models for known conditional interactions
cl = makeCluster(ncores, outfile = '')
registerDoParallel(cl)
fitsKnown <- foreach(i = bimodal, .combine = rbind, .packages = c('mclust')) %dopar% {
  cond = Mclust(datamat[i, ], G = 2)$classification
  res1 = fitModels(t(datamat[, cond == 1, drop = F]), t(datamat[, cond == 1, drop = F]))
  res2 = fitModels(t(datamat[, cond == 2, drop = F]), t(datamat[, cond == 2, drop = F]))
  colnames(res1)[3:5] = paste0(colnames(res1)[3:5], '.1')
  colnames(res2)[3:5] = paste0(colnames(res2)[3:5], '.2')
  res = cbind('cond' = i,res1, res2[, 3:5])
}
stopCluster(cl)

#fit linear models for random conditional interactions
cl = makeCluster(ncores, outfile = '')
registerDoParallel(cl)
fitsRandom <- foreach(i = 1:10, .combine = rbind, .packages = c('mclust')) %dopar% {
  res = c()
  for (j in bimodal) {
    bimmodel = Mclust(datamat[j, ], G = 1:2)
    cond = bimmodel$classification
    cond = sample(cond, length(cond))
    res1 = fitModels(t(datamat[, cond == 1, drop = F]), t(datamat[, cond == 1, drop = F]))
    res2 = fitModels(t(datamat[, cond == 2, drop = F]), t(datamat[, cond == 2, drop = F]))
    colnames(res1)[3:5] = paste0(colnames(res1)[3:5], '.1')
    colnames(res2)[3:5] = paste0(colnames(res2)[3:5], '.2')
    res = rbind(res, cbind('cond' = j, res1, res2[, 3:5]))
  }
  return(res)
}
stopCluster(cl)
# fitsRandom = fitsRandom[!(fitsRandom$y %in% bimtargets & fitsRandom$x %in% bimtargets), ]
fitsRandom = fitsRandom[fitsRandom$y != fitsRandom$x, ]
# fitsRandom = fitsRandom[!fitsRandom$y %in% condcoex & !fitsRandom$x %in% condcoex, ]

#restructure data for plotting
#for known interactions
fts1 = fitsKnown[, 1:6]
fts2 = fitsKnown[, c(1:3, 7:9)]
colnames(fts1)[4:6] = colnames(fts2)[4:6] = c('m', 'c', 'p')
fts = melt(fts1, id.vars = c('cond', 'x', 'y'))
fts$condition2 = melt(fts2, id.vars = c('cond', 'x', 'y'))$value
colnames(fts)[5] = 'condition1'
ftsKnown = fts
remove(fts1, fts2, fts)
ftsKnown = merge(triplets, ftsKnown)

#for random interactions
fitsRandom = fitsRandom[sample.int(nrow(fitsRandom), min(0.5E6, nrow(fitsRandom))), ] # down sample
fts1 = fitsRandom[, 1:6]
fts2 = fitsRandom[, c(1:3, 7:9)]
colnames(fts1)[4:6] = colnames(fts2)[4:6] = c('m', 'c', 'p')
fts = melt(fts1, id.vars = c('cond', 'x', 'y'))
fts$condition2 = melt(fts2, id.vars = c('cond', 'x', 'y'))$value
colnames(fts)[5] = 'condition1'
ftsRandom = fts
remove(fts1, fts2, fts)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plotdf = ftsKnown
plotdf = plotdf[as.character(plotdf$x) != as.character(plotdf$y),]

#compute sample numbers in condition 1 and then the z-score
sampcount <- foreach(i = bimodal, .combine = rbind) %do% {
  numsamp1 = sum(Mclust(datamat[i, ], G = 2)$classification == 1)
  return(data.frame('cond' = i, 'count' = numsamp1))
}
plotdf = merge(plotdf, sampcount)
plotdf$zs = 0
ps = plotdf[plotdf$variable %in% 'p',]
zs = (atanh(ps$condition1) - atanh(ps$condition2)) / sqrt(1/(ps$count - 3) + 1/(nsamp - ps$count - 3))
plotdf$zs[plotdf$variable %in% 'm'] = plotdf$zs[plotdf$variable %in% 'c'] = plotdf$zs[plotdf$variable %in% 'p'] = zs
bgzs = (atanh(fitsRandom$p.1) - atanh(fitsRandom$p.2)) / sqrt(1/(50 - 3) + 1/(nsamp - 50 - 3))

#text for interactive plots
plotdf$text = apply(plotdf[,1:5], 1, function (x) paste(paste(colnames(plotdf)[1:5], as.vector(x), sep = ': '), collapse = '\n'))

#plot correlation coefficients
textSize = 1.2
rng = c(-1, 1) * ceiling(max(abs(plotdf$zs)))
pp = ggplot(plotdf[plotdf$variable %in% 'p', ]) +
  geom_point(aes(condition1, condition2, color = zs, shape = cond, text = text), size = 2) +
  scale_colour_distiller(palette = 'PRGn', direction = 1, limits = rng) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_abline() +
  labs(subtitle = 'Pearson rho') + xlab('Condition 1') + ylab('Condition 2') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(textSize)),
    axis.text.x = element_text(angle = 0, size = rel(textSize)),
    axis.text.y = element_text(angle = 0, size = rel(textSize)),
    strip.background = element_rect(colour = "gray80", fill = "gray80"),
    strip.text = element_text(size = rel(textSize)),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(face = "italic"),
    plot.title = element_text(
      face = "bold",
      size = rel(textSize),
      hjust = 0.5
    ),
    plot.subtitle = element_text(size = rel(textSize))
  )
pp = pp + geom_density2d(data = ftsRandom[ftsRandom$variable %in% 'p', ], aes(condition1, condition2))

pm = ggplot(plotdf[plotdf$variable %in% 'm', ]) +
  geom_point(aes(condition1, condition2, color = zs, shape = cond, text = text), size = 2) +
  scale_colour_distiller(palette = 'PRGn', direction = 1, limits = rng) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_abline() +
  labs(subtitle = 'slope of lm') + xlab('Condition 1') + ylab('Condition 2') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(textSize)),
    axis.text.x = element_text(angle = 0, size = rel(textSize)),
    axis.text.y = element_text(angle = 0, size = rel(textSize)),
    strip.background = element_rect(colour = "gray80", fill = "gray80"),
    strip.text = element_text(size = rel(textSize)),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(face = "italic"),
    plot.title = element_text(
      face = "bold",
      size = rel(textSize),
      hjust = 0.5
    ),
    plot.subtitle = element_text(size = rel(textSize))
  ) +
  xlim(-2,2) +
  ylim(-2,2)
pm = pm + geom_density2d(data = ftsRandom[ftsRandom$variable %in% 'm', ], aes(condition1, condition2))

pc = ggplot(plotdf[plotdf$variable %in% 'c', ]) +
  geom_point(aes(condition1, condition2, color = zs, shape = cond, text = text), size = 2) +
  scale_colour_distiller(palette = 'PRGn', direction = 1, limits = rng) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_abline() +
  labs(subtitle = 'y-intercept of lm') + xlab('Condition 1') + ylab('Condition 2') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(textSize)),
    axis.text.x = element_text(angle = 0, size = rel(textSize)),
    axis.text.y = element_text(angle = 0, size = rel(textSize)),
    strip.background = element_rect(colour = "gray80", fill = "gray80"),
    strip.text = element_text(size = rel(textSize)),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(face = "italic"),
    plot.title = element_text(
      face = "bold",
      size = rel(textSize),
      hjust = 0.5
    ),
    plot.subtitle = element_text(size = rel(textSize))
  ) +
  xlim(-2,2) +
  ylim(-2,2)
pc = pc + geom_density2d(data = ftsRandom[ftsRandom$variable %in% 'c', ], aes(condition1, condition2))

ptitle = ggplot() + ggtitle('Conditional interaction characteristics in simulations') + 
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = rel(textSize + 0.5),
      hjust = 0.5
    )
  )

gridsz = 10
layout = c(rep(c(1, rep(2, gridsz)), gridsz),
           rep(c(1, rep(3, gridsz)), gridsz),
           rep(c(1, rep(4, gridsz)), gridsz))
layout = matrix(layout, ncol = 3 * gridsz)

# png(filename = 'figures/sim_correlation_profiles_EColifull_bimodal0.1_colour.png', width = 12, height = 6, units = 'in', res = 300)
multiplot(ptitle, pm, pc, pp, layout = layout)
# dev.off()


