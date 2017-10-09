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

#----read network----
df = read.csv('sourceNets/EColi_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

#read in network and create simulator object
grnEColi = df2GraphGRN(edges, loops = F, propor = 0.1, seed = 36)
# grnEColi = sampleGraph(grnEColi, 100, 15, seed = 26363)
grnEColi = randomizeParams(grnEColi, 'linear-like', seed = 36)
simEColi =new('SimulationGRN', graph = grnEColi, seed = 36, propBimodal = 0, expnoise = 0, bionoise = 0.1)

#----generate data----
#Total simulation params
nsamp = 100
simseed = 438756928
set.seed(simseed)

#modify input models
modinput = 'crp'
prop = runif(1, 0.2, 0.8)
prop = c(prop, 1 - prop)
mu = c(runif(1, 0.1, 0.5), runif(1, 0.5, 0.9))
maxsd = pmin(mu, 1 - mu) / 3
sdev = sapply(maxsd, function(x) runif(1, 0.01, x))

simEColi$inputModels[[modinput]] = list('prop' = prop,
                                 'mean' = mu,
                                 'sd' = sdev)
#generate data and add noise
ncores = 10
cl = makeCluster(ncores, outfile = '')
clusterExport(cl, c('addBioNoise'))
registerDoSNOW(cl)
datamat = simulateDataset(simEColi, nsamp)
stopCluster(cl)
noisydatamat = addNoise(simEColi, datamat)

#----plot networks----
plotInference <- function(graph, modulator, scores, topn = 100, truth = NULL) {
  scoredf = cbind(expand.grid('from' = rownames(datamat), 'to' = rownames(datamat)),
                  'score' = as.numeric(scores))
  scoredf['truth'] = NA
  if (!is.null(truth)) {
    truth = truth[rownames(zsc), colnames(zsc)]
    scoredf$truth = as.numeric(truth)
  }
  
  dflist = GraphGRN2df(graph)
  dflist$edges = merge(dflist$edges, scoredf, all = T, by = c('from', 'to'))
  dflist$edges[is.na(dflist$edges$type), 'type'] = 'inferred'
  
  net = graph_from_data_frame(d = dflist$edges, directed = T, vertices = dflist$nodes)
  V(net)$size = log(degree(net - E(net)[E(net)$type %in% 'inferred'], mode = 'out') + 3) * 2
  V(net)$label = NA
  V(net)$color = 1
  V(net)$color[V(net)$name %in% modulator] = 2
  E(net)$arrow.size = .075
  E(net)$color = 'gray36'
  E(net)$color[!is.na(E(net)$activation) & E(net)$activation] = 'darkgreen'
  E(net)$color[!is.na(E(net)$activation) & !E(net)$activation] = 'red'
  E(net)$color[!is.na(E(net)$truth) & E(net)$truth == 1] = 'blue'
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

#generate a sensitivity matrix
iEcoli = GraphGRN2igraph(grnEColi)
inputs = sapply(simEColi$inputModels, function(x) x$mean[1])
pertbNodes = adjacent_vertices(iEcoli, modinput)[[1]]$name
# pertbNodes = unique(c(pertbNodes, adjacent_vertices(iEcoli, pertbNodes)[[1]]$name))
pertbNodes = nodenames(grnEColi)

cl = makeCluster(ncores, outfile = '')
clusterExport(cl, c('addBioNoise'))
registerDoSNOW(cl)
inputs[modinput] = mu[2]
sensmat1 = sensitivityAnalysis(simEColi, 0.25, inputs, pertbNodes)
sensmat2 = sensitivityAnalysis(simEColi, mu[2] - mu[1], inputs, modinput)
sensmat2 = rep(1,nrow(sensmat1)) %*% sensmat2
stopCluster(cl)

sensthresh = 0.1
s4 = sensmat1 * sensmat2
s4 = s4[!rownames(s4) %in% modinput, ]
diffpairs = expand.grid('TF' = rownames(s4), 'Target' = colnames(s4), stringsAsFactors = F)
diffpairs$sensitivity = as.numeric(s4)
diffpairs = diffpairs[abs(diffpairs$sensitivity) > sensthresh, ]
diffpairs = diffpairs[order(diffpairs$TF), ]
diffpairs = diffpairs[diffpairs$TF != diffpairs$Target,]
rownames(diffpairs) = NULL

#exhaust all pairs and include all targets of TFs
truthmat = matrix(0, nrow = length(nodenames(grnEColi)), ncol = length(nodenames(grnEColi)))
rownames(truthmat) = colnames(truthmat) = nodenames(grnEColi)
for(i in 1:nrow(diffpairs)) {
  tf = diffpairs[i, 'TF']
  tgt = diffpairs[i, 'Target']
  tf = colnames(sensmat1)[abs(sensmat1[tf, ]) > sensthresh * 2]
  truthmat[tf, tgt] = 1
  truthmat[tgt, tf] = 1
}
diag(truthmat) = 0

#make inference
classf = Mclust(noisydatamat[modinput,], verbose = F)$classification
zsc = z.score.sp(noisydatamat, classf)
zsc[is.infinite(zsc)] = 0

cl = makeCluster(ncores)
registerDoParallel(cl)
ftgisc = ftgi.score(noisydatamat, classf)
stopCluster(cl)

ebsc = ebcoexpress.score(noisydatamat, classf)
caisc = cai.score(noisydatamat, classf, thresh.method = 'score')
magicsc = magic.score(noisydatamat, classf)
dicersc = dicer.score(noisydatamat, classf)
diffcoexsc = diffcoex.score(noisydatamat, classf)
mindysc = mindy.score(noisydatamat, classf, ncores = ncores)

truthmat = truthmat[rownames(zsc), colnames(zsc)]
#store output of runs
topn = 100
fname = paste0('figures/LN_infPerfTOP', topn, '_', modinput, '_', simseed, '.pdf')
pdf(file = fname, width = 12, onefile = T)

par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
plotInference(grnEColi, modinput, truthmat, sum(truthmat))
mtext("Truth", outer = T, cex = 1.5)

# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, abs(zsc), topn, truthmat)
# mtext("Z-score", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, ebsc, topn, truthmat)
# mtext("EBCo-express", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, mindysc, topn)
# mtext("MINDy", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, magicsc, topn)
# mtext("MAGIC", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, caisc, topn)
# mtext("Cai", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, dicersc, topn)
# mtext("DICER", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, diffcoexsc, topn)
# mtext("DiffCoEx", outer = T, cex = 1.5)
# 
# par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
# plotInference(grnEColi, modinput, ftgisc, topn)
# mtext("FTGI", outer = T, cex = 1.5)
par(mfrow=c(1,1), oma = c(0,0,0,0))

#ROC analysis
zroc = roc(as.numeric(truthmat),abs(as.numeric(zsc)))
ebroc = roc(as.numeric(truthmat),abs(as.numeric(ebsc)))
mindyroc = roc(as.numeric(truthmat),abs(as.numeric(mindysc)))
magicroc = roc(as.numeric(truthmat),abs(as.numeric(magicsc)))
cairoc = roc(as.numeric(truthmat),abs(as.numeric(caisc)))
dicerroc = roc(as.numeric(truthmat),abs(as.numeric(dicersc)))
diffcoexroc = roc(as.numeric(truthmat),abs(as.numeric(diffcoexsc)))
ftgiroc = roc(as.numeric(truthmat),abs(as.numeric(ftgisc)))

rocdata = data.frame('Method' = 'z-score', 'Specificity' = zroc$specificities, 'Sensitivity' = zroc$sensitivities)
rocdata = rbind(rocdata, data.frame('Method' = 'EBCo-express', 'Specificity' = ebroc$specificities, 'Sensitivity' = ebroc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'MINDy', 'Specificity' = mindyroc$specificities, 'Sensitivity' = mindyroc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'MAGIC', 'Specificity' = magicroc$specificities, 'Sensitivity' = magicroc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'Cai', 'Specificity' = cairoc$specificities, 'Sensitivity' = cairoc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'DICER', 'Specificity' = dicerroc$specificities, 'Sensitivity' = dicerroc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'DiffCoEx', 'Specificity' = diffcoexroc$specificities, 'Sensitivity' = diffcoexroc$sensitivities))
rocdata = rbind(rocdata, data.frame('Method' = 'FTGI', 'Specificity' = ftgiroc$specificities, 'Sensitivity' = ftgiroc$sensitivities))
ggplot(rocdata, aes(1 - Specificity, Sensitivity, colour = Method)) + geom_line() +
  scale_colour_brewer(palette = 'Paired') +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  ggtitle('Inference methods ROC curves') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1)),
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(size = rel(1.2)),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(face="italic", size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    plot.title = element_text(face = "bold", size = rel(1.5), hjust = 0.5)
  )

dev.off()

aucs = c(zroc$auc, ebroc$auc, mindyroc$auc, magicroc$auc, cairoc$auc, dicerroc$auc, diffcoexroc$auc, ftgiroc$auc)
names(aucs) = c('z-score', 'EBCo-express', 'MINDy', 'MAGIC', 'Cai', 'DICER', 'DiffCoEx', 'FTGI')







# 
# simEColi$noiseG = 0.1
# simEColi$noiseL = 0
# noisydatamat = addNoise(simEColi, datamat)
# 
# #plot inference results
# df = as.data.frame(t(noisydatamat))
# df$cond = Mclust(df[,modinput], verbose = F)$classification
# ggplot(df, aes(cpxAR, skp_lpxDA_fabZ, colour = rpoE_rseABC)) + geom_point() + facet_wrap(~cond)
# ggplot(df, aes(rpoH, cytR, colour = rpoE_rseABC)) + geom_point() + facet_wrap(~cond)
# ggplot(df, aes(rpoE_rseABC, mglBAC, colour = purR)) + geom_point() + facet_wrap(~cond)
# 
# expr1=noisydatamat[,classf==1]
# expr2=noisydatamat[,classf==2]
# r1=cor(t(expr1), method = 'spearman')
# r2=cor(t(expr2), method = 'spearman')
# plotdf = data.frame('r1' = as.numeric(r1), 'r2' = as.numeric(r2), 'z' = as.numeric(abs(zsc)))
# plotdf = cbind(plotdf, expand.grid(rownames(datamat), rownames(datamat)))
# ggplot(plotdf[plotdf$z > 3.5,], aes(r1, r2, colour = z)) + geom_point()

