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
library(reshape2)

#----read network----
df = read.csv('sourceNets/EColi_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

#read in network and randomize params
grnEColi = df2GraphGRN(edges, loops = F, propor = 0.1, seed = 36)
grnEColi = randomizeParams(grnEColi, 'linear-like', seed = 36)

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
  V(net)$size = log(igraph::degree(net - E(net)[E(net)$type %in% 'inferred'], mode = 'out') + 3) * 2
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

relent <- function(truth, prediction){
  D = truth
  M = prediction
  term1 = suppressWarnings(D*log(D/M))
  term2 = suppressWarnings((1 - D)*log((1 - D)/(1 - M)))
  term1[is.nan(term1)] = 0
  term2[is.nan(term2)] = 0
  H = sum(term1+term2)
  return(H)
}

#----create small graph----
originseed = 360
set.seed(originseed)
netSize = 50
minTFs = 10
expnoise = 0
bionoise = 0.1
simseeds = sample.int(1E7, 25)
ncores = 10
allaucs = c()
allrelents = c()

topn = 100
fname = paste0('figures/Nruns', '_', originseed, '_', netSize, '_', minTFs,'_N', length(simseeds), '.pdf')
pdf(file = fname, width = 12, onefile = T)

for (i in 1:length(simseeds)) {
  message(i)
  grnSmall = sampleGraph(grnEColi, netSize, minTFs, seed = simseeds[i])
  simSmall =new('SimulationGRN', graph = grnSmall, seed = simseeds[i], propBimodal = 0, expnoise = expnoise, bionoise = bionoise)
  
  #select modulator
  iSmall = GraphGRN2igraph(grnSmall)
  modulatedNodes = degree(iSmall, mode = 'in')
  modulatedNodes = names(modulatedNodes)[modulatedNodes>=2]
  
  modinput = rowSums(distances(iSmall)[,modulatedNodes]!=0)
  modinput = names(modinput)[modinput != 0]
  modinput = intersect(getInputNodes(grnSmall), modinput)
  modinput = sample(modinput, 1)
  
  #----generate data----
  #Total simulation params
  nsamp = 100
  set.seed(simseeds[i])
  
  #modify input models
  prop = runif(1, 0.3, 0.7)
  prop = c(prop, 1 - prop)
  mu = c(runif(1, 0.1, 0.5), runif(1, 0.5, 0.9))
  maxsd = pmin(mu, 1 - mu) / 3
  sdev = sapply(maxsd, function(x) runif(1, 0.01, x))
  
  simSmall$inputModels[[modinput]] = list('prop' = prop,
                                          'mean' = mu,
                                          'sd' = sdev)
  #generate data and add noise
  cl = makeCluster(ncores, outfile = '')
  clusterExport(cl, c('addBioNoise'))
  registerDoSNOW(cl)
  datamat = simulateDataset(simSmall, nsamp)
  stopCluster(cl)
  noisydatamat = addNoise(simSmall, datamat)
  
  classf = Mclust(noisydatamat[modinput,], verbose = F)$classification
  if (length(unique(classf)) != 2) {
    warning('Non-bimodal condition encountered')
    next
  }
  
  #generate a sensitivity matrix
  iEcoli = GraphGRN2igraph(grnSmall)
  inputs = sapply(simSmall$inputModels, function(x) x$mean[1])
  pertbNodes = adjacent_vertices(iEcoli, modinput)[[1]]$name
  pertbNodes = nodenames(grnSmall)
  
  cl = makeCluster(ncores, outfile = '')
  clusterExport(cl, c('addBioNoise'))
  registerDoSNOW(cl)
  inputs[modinput] = mu[2]
  sensmat1 = sensitivityAnalysis(simSmall, 0.25, inputs, pertbNodes)
  sensmat2 = sensitivityAnalysis(simSmall, mu[2] - mu[1], inputs, modinput)
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
  
  if (nrow(diffpairs) == 0) {
    warning('No differential signal present')
    next
  }
  
  #exhaust all pairs and include all targets of TFs
  truthmat = matrix(0, nrow = length(nodenames(grnSmall)), ncol = length(nodenames(grnSmall)))
  rownames(truthmat) = colnames(truthmat) = nodenames(grnSmall)
  for(i in 1:nrow(diffpairs)) {
    tf = diffpairs[i, 'TF']
    tgt = diffpairs[i, 'Target']
    tf = colnames(sensmat1)[abs(sensmat1[tf, ]) > sensthresh * 2]
    truthmat[tf, tgt] = 1
    truthmat[tgt, tf] = 1
  }
  diag(truthmat) = 0
  
  #make inference
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
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  plotInference(grnSmall, modinput, truthmat, sum(truthmat))
  mtext(paste0("Truth: ", modinput), outer = T, cex = 1.5)
  par(mfrow=c(1,1), oma = c(0,0,0,0))
  
  #ROC analysis
  zroc = roc(as.numeric(truthmat),abs(as.numeric(zsc)))
  ebroc = roc(as.numeric(truthmat),as.numeric(ebsc))
  mindyroc = roc(as.numeric(truthmat),abs(as.numeric(mindysc)))
  magicroc = roc(as.numeric(truthmat),abs(as.numeric(magicsc)))
  cairoc = roc(as.numeric(truthmat),abs(as.numeric(caisc)))
  dicerroc = roc(as.numeric(truthmat),abs(as.numeric(dicersc)))
  diffcoexroc = roc(as.numeric(truthmat),as.numeric(diffcoexsc))
  ftgiroc = roc(as.numeric(truthmat),as.numeric(ftgisc))
  
  rocdata = data.frame('Method' = 'z-score', 'Specificity' = zroc$specificities, 'Sensitivity' = zroc$sensitivities)
  rocdata = rbind(rocdata, data.frame('Method' = 'EBCo-express', 'Specificity' = ebroc$specificities, 'Sensitivity' = ebroc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'MINDy', 'Specificity' = mindyroc$specificities, 'Sensitivity' = mindyroc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'MAGIC', 'Specificity' = magicroc$specificities, 'Sensitivity' = magicroc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'Cai', 'Specificity' = cairoc$specificities, 'Sensitivity' = cairoc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'DICER', 'Specificity' = dicerroc$specificities, 'Sensitivity' = dicerroc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'DiffCoEx', 'Specificity' = diffcoexroc$specificities, 'Sensitivity' = diffcoexroc$sensitivities))
  rocdata = rbind(rocdata, data.frame('Method' = 'FTGI', 'Specificity' = ftgiroc$specificities, 'Sensitivity' = ftgiroc$sensitivities))
  p = ggplot(rocdata, aes(1 - Specificity, Sensitivity, colour = Method)) + geom_line() +
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
  
  print(p)
  aucs = c(zroc$auc, ebroc$auc, mindyroc$auc, magicroc$auc, cairoc$auc, dicerroc$auc, diffcoexroc$auc, ftgiroc$auc)
  names(aucs) = c('z-score', 'EBCo-express', 'MINDy', 'MAGIC', 'Cai', 'DICER', 'DiffCoEx', 'FTGI')
  allaucs = rbind(allaucs, aucs)
  
  #Relative entropy analysis
  zrelent = relent(as.numeric(truthmat),abs(as.numeric(zsc)))
  ebrelent = relent(as.numeric(truthmat),as.numeric(ebsc))
  mindyrelent = relent(as.numeric(truthmat),abs(as.numeric(mindysc)))
  magicrelent = relent(as.numeric(truthmat),abs(as.numeric(magicsc)))
  cairelent = relent(as.numeric(truthmat),abs(as.numeric(caisc)))
  dicerrelent = relent(as.numeric(truthmat),abs(as.numeric(dicersc)))
  diffcoexrelent = relent(as.numeric(truthmat),as.numeric(diffcoexsc))
  ftgirelent = relent(as.numeric(truthmat),as.numeric(ftgisc))
  relents = c(zrelent, ebrelent, mindyrelent, magicrelent, cairelent, dicerrelent, diffcoexrelent, ftgirelent)
  names(relents) = c('z-score', 'EBCo-express', 'MINDy', 'MAGIC', 'Cai', 'DICER', 'DiffCoEx', 'FTGI')
  
  allrelents = rbind(allrelents, relents)
}

#plot boxplot of all aucs
allaucs = melt(allaucs)[, 2:3]
colnames(allaucs) = c('Method', 'AUROC')
p = ggplot(allaucs, aes(Method, AUROC, colour = Method)) + geom_boxplot() +
  scale_colour_brewer(palette = 'Paired') +
  ggtitle('Inference methods AUROC values') +
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
print(p)

#plot boxplot of all relents
allrelents = melt(allrelents)[, 2:3]
colnames(allrelents) = c('Method', 'Relative Entropy')
p = ggplot(allrelents, aes(Method, `Relative Entropy`, colour = Method)) + geom_boxplot() +
  scale_colour_brewer(palette = 'Paired') +
  ggtitle('Inference methods AUROC values') +
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
print(p)

dev.off()

