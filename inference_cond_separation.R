setwd('/wehisan/home/allstaff/b/bhuva.d/R_projects/cGRNsimulator')
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
df = read.csv('sourceNets/EColi_full.sif', sep = '\t', header = F, stringsAsFactors = F)
edges = df
maplogical = c('ac' = T, 're' = F, 'du' = F)
edges[,2] = maplogical[edges[,2]]

#read in network and randomize params
grnEColi = df2GraphGRN(edges, loops = F, propor = 0.1, seed = originseed)
grnEColi = randomizeParams(grnEColi, 'linear-like', seed = originseed)

#----plot networks----
plotInference <- function(graph, modulator, scores, topn = 100, truth = NULL) {
  #convert scores matrix to dataframe
  scoredf = cbind(expand.grid('from' = rownames(scores), 'to' = rownames(scores),
                              stringsAsFactors = T), 'score' = as.numeric(scores))
  scoredf['truth'] = NA
  if (!is.null(truth)) {
    truth = truth[rownames(scores), colnames(scores)]
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
  infnet = net - E(net)[rank(-E(net)$score, ties.method = 'min') > topn]
  set.seed(360)
  l = layout_with_fr(goldnet)
  
  plot(goldnet, layout = l,
       mark.group = list(c(modulator, adjacent_vertices(goldnet, modulator)[[1]]$name)),
       main = 'Gold standard network')
  plot(infnet, layout = l,
       mark.group = list(c(modulator, adjacent_vertices(goldnet, modulator)[[1]]$name)),
       main = 'Inferred network')
}

scVectorize <- function(scoremat, modin) {
  modpos = which(rownames(scoremat) %in% modin)
  scoremat = scoremat[-modpos, -modpos]
  scoremat = pmax(scoremat, t(scoremat))
  scoremat = scoremat[upper.tri(scoremat)]
  
  return(scoremat)
}

#----Run multiple simulations----
numsims = 300
numpertb = 1
nsamp = 100
netSize = 80
minTFs = 10
pdfname = paste0('simdata/sep_INF_', originseed, '_', netSize, '_', minTFs,'_N', numsims, '_ecoli.pdf')
rdname = paste0('simdata/separation_INF_', originseed, '_', netSize, '_Ptb', numpertb, '_N', numsims, '_ecoli.RData')

ncores = detectCores()
set.seed(originseed)
allaucs = c()
allsingscores = c()
modnames = c()
simseeds = sample.int(1E7, numsims)

pdf(file = pdfname, width = 12, onefile = T)

#simulation specific parameters
expnoise = seq(0, 0, length.out = numpertb)
bionoise = seq(0.0/3, 0.0/3, length.out = numpertb)
mus = expand.grid('mu1' = 0.3,
                  'mu2' = seq(0.7, 0.7, length.out = numpertb))
sdevs = expand.grid('sd1' = 0.3,
                    'sd2' = seq(0.3, 0.3, length.out = numpertb))/3
xs = seq(1, 1, length.out = numpertb)
props = data.frame('p1' = xs/(1 + xs),
                   'p2' = 1 - xs/(1 + xs))

for (moditr in 1:numsims) {
  set.seed(simseeds[moditr])
  message(paste0('Modulator iter: ', moditr))
  
  #sample a smaller network OR use full
  grnSmall = sampleGraph(grnEColi, netSize, minTFs, seed = simseeds[moditr])
  grnSmall = randomizeParams(grnSmall, 'linear-like', simseeds[moditr])
  
  #select a modulator
  iSmall = GraphGRN2igraph(grnSmall)
  modulatedNodes = degree(iSmall, mode = 'in')
  modulatedNodes = names(modulatedNodes)[modulatedNodes>=2]
  moddists = distances(iSmall)
  moddists[is.infinite(moddists)] = NA
  modinput = rowSums(moddists[,modulatedNodes]!=0, na.rm = T)
  modinput = names(modinput)[modinput != 0]
  modinput = intersect(getInputNodes(grnSmall), modinput)
  modinput = sample(modinput, 1)
  
  for (itr in 1:numpertb) {
    message(paste0('Iter: ', itr))
    simSmall =new('SimulationGRN', graph = grnSmall, seed = simseeds[moditr], propBimodal = 0, expnoise = expnoise[itr], bionoise = bionoise[itr])
    
    #----generate data----
    #Total simulation params
    set.seed(simseeds[moditr])
    
    #modify input models
    prop = as.numeric(props[itr,])
    mu = as.numeric(mus[itr,])
    sdev = as.numeric(sdevs[itr,])
    
    simSmall$inputModels[[modinput]] = list('prop' = prop,
                                            'mean' = mu,
                                            'sd' = sdev)
    #generate data and add noise
    cl = makeCluster(ncores, outfile = '')
    clusterExport(cl, c('addBioNoise'))
    registerDoSNOW(cl)
    datamat = simulateDataset(simSmall, nsamp)
    stopCluster(cl)
    
    classf = attr(datamat, 'classf')[modinput, ]
    if (length(unique(classf)) != 2) {
      warning('Non-bimodal condition encountered')
      next
    }
    
    #----generate a truth adjacency matrix----
    iEcoli = GraphGRN2igraph(grnSmall)
    inputs = sapply(simSmall$inputModels, function(x) x$mean[1])
    #pertbNodes = adjacent_vertices(iEcoli, modinput)[[1]]$name
    pertbNodes = nodenames(grnSmall)
    
    cl = makeCluster(ncores, outfile = '')
    clusterExport(cl, c('addBioNoise'))
    registerDoSNOW(cl)
    inputs[modinput] = mu[2]
    sensmat1 = sensitivityAnalysis(simSmall, 0.25, inputs, pertbNodes)
    sensmat2 = sensitivityAnalysis(simSmall, mu[2] - mu[1], inputs, modinput)
    sensmat2 = rep(1,nrow(sensmat1)) %*% sensmat2
    stopCluster(cl)
    
    sensthresh = 0.05
    sensmat1[abs(sensmat1) < sensthresh] = 0
    sensmat2[abs(sensmat2) < sensthresh] = 0
    sensmat1 = sensmat1[colnames(sensmat2), colnames(sensmat2)]
    sensmat1[sensmat2[1, ] != 0, ] = 0
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
      truthmat = matrix(0, nrow = length(nodenames(grnSmall)), ncol = length(nodenames(grnSmall)))
      rownames(truthmat) = colnames(truthmat) = nodenames(grnSmall)
      next
    }
    
    #exhaust all pairs and include all targets of TFs
    truthmat = matrix(0, nrow = length(nodenames(grnSmall)), ncol = length(nodenames(grnSmall)))
    rownames(truthmat) = colnames(truthmat) = nodenames(grnSmall)
    for(j in 1:nrow(diffpairs)) {
      tf = diffpairs[j, 'TF']
      tgt = diffpairs[j, 'Target']
      tf = colnames(sensmat1)[abs(sensmat1[tf, ]) > sensthresh * 2]
      truthmat[tf, tgt] = 1
      truthmat[tgt, tf] = 1
    }
    diag(truthmat) = 0
    
    #----make inference----
    scorelist = c()
    zsc = z.score(datamat, classf)
    zsc[is.infinite(zsc)] = 0
    zspsc = z.score.sp(datamat, classf)
    zspsc[is.infinite(zspsc)] = 0
    cl = makeCluster(ncores)
    registerDoParallel(cl)
    ftgisc = ftgi.score(datamat, classf)
    stopCluster(cl)
    ebsc = ebcoexpress.score(datamat, classf)
    caisc = cai.score(datamat, classf, thresh.method = 'score')
    magicsc = magic.score(datamat, classf)
    dicersc = dicer.score(datamat, classf)
    diffcoexsc = diffcoex.score(datamat, classf)
    ecfsc = ecf.score(datamat, classf)
    entsc = ent.score(datamat, classf)
    # lasc = la.score(datamat, datamat[modinput,])
    # mindysc = mindy.score(datamat, classf, ncores = ncores)
    
    truthmat = truthmat[rownames(zsc), colnames(zsc)]#reorder names in truth matrix
    scorelist = list('z-score' = scVectorize(abs(zsc), modinput),
                     'z-score-s' = scVectorize(abs(zspsc), modinput),
                     'FTGI' = scVectorize(abs(ftgisc), modinput),
                     'EBcoexpress' = scVectorize(ebsc, modinput),
                     'Cai' = scVectorize(abs(caisc), modinput),
                     'MAGIC' = scVectorize(abs(magicsc), modinput),
                     'DICER' = scVectorize(abs(dicersc), modinput),
                     'DiffCoEx' = scVectorize(diffcoexsc, modinput),
                     'ECF' = scVectorize(ecfsc, modinput),
                     'ENT' = scVectorize(entsc, modinput)
                     # 'LA*' = scVectorize(abs(lasc), modinput)
                     # 'MINDy' = as.numeric(mindysc)
    )
    
    #network properties
    indens = edge_density(graph_from_adjacency_matrix(truthmat))
    hubscore = hub_score(iSmall)$vector[modinput]
    ec = eigen_centrality(iSmall)$vector[modinput]
    dens = edge_density(iSmall)
    tfsds = mean(apply(datamat[diffpairs$TF, , drop = F], 1, sd))
    
    #ROC evaluation
    roclist = lapply(scorelist, function(x) roc(scVectorize(truthmat, modinput), x, direction = '<'))
    aucs = lapply(roclist, function(x) x$auc)
    allaucs = rbind(
      allaucs,
      cbind(
        unlist(aucs),
        'MuDiff' = mu[2] - mu[1],
        'VarRatio' = sdev[2] / sdev[1],
        'PropRatio' = prop[2] / prop[1],
        'ExpNoise' = expnoise[itr],
        'BioNoise' = bionoise[itr],
        'DensityInf' = indens,
        'HubScore' = hubscore,
        'EigenCent' = ec,
        'Density' = dens,
        'AvgSensitivity' = mean(abs(diffpairs$sensitivity)),
        'AvgTFsd' = tfsds
      )
    )
    
    #mean rank of true positives
    singscore = lapply(scorelist, function(x) mean(rank(x)[as.logical(scVectorize(truthmat, modinput))]/length(x)))
    allsingscores = rbind(
      allsingscores,
      cbind(
        unlist(singscore),
        'MuDiff' = mu[2] - mu[1],
        'VarRatio' = sdev[2] / sdev[1],
        'PropRatio' = prop[2] / prop[1],
        'ExpNoise' = expnoise[itr],
        'BioNoise' = bionoise[itr],
        'DensityInf' = indens,
        'HubScore' = hubscore,
        'EigenCent' = ec,
        'Density' = dens,
        'AvgSensitivity' = mean(abs(diffpairs$sensitivity)),
        'AvgTFsd' = tfsds
      )
    )
    
    modnames = c(modnames, paste(modinput, moditr, 'sep' = '_'))
  }
  
  #store output of runs
  message(modinput)
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  plotInference(grnSmall, modinput, truthmat, sum(truthmat))
  mtext(paste0("Truth: ", modinput), outer = T, cex = 1.5)
  par(mfrow=c(1,1), oma = c(0,0,0,0))
  
  #regular saving
  if (moditr %% 50 == 0) {
    oldaucs = allaucs
    oldsingscores = allsingscores
    allaucs = as.data.frame(allaucs, 'Modulator' = modnames)
    allsingscores = as.data.frame(allsingscores, 'Modulator' = modnames)
    
    save(allaucs, allsingscores, originseed, numpertb, numsims,
         expnoise, bionoise, mus, sdevs, props, modnames, file = rdname)
    allaucs = oldaucs
    allsingscores = oldsingscores
  }
}

dev.off()

#add modulator info to measurement dfs
allaucs = as.data.frame(allaucs, 'Modulator' = modnames)
allsingscores = as.data.frame(allsingscores, 'Modulator' = modnames)

save(allaucs, allsingscores, originseed, numpertb, numsims,
     expnoise, bionoise, mus, sdevs, props, modnames, file = rdname)

# #----plot results----
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   #data_sum <- rename(data_sum, c("mean" = varname))
#   colnames(data_sum)[colnames(data_sum) %in% 'mean'] = varname
#   return(data_sum)
# }
# 
# allaucs = data.frame(allaucs, rownames(allaucs), 'AUROC')
# allMIs = data.frame(allMIs, rownames(allMIs), 'MI')
# allrelents = data.frame(allrelents, rownames(allrelents), 'RelativeEntropy')
# allsingscores = data.frame(allsingscores, rownames(allsingscores), 'singscore')
# colnames(allaucs) = colnames(allMIs) = colnames(allrelents) =
#   colnames(allsingscores) = c('Measurement', 'Distance', 'Difference', 'VarRatio', 'PropRatio', 'Method', 'Measure')
# 
# measures = rbind(allaucs, allMIs, allrelents, allsingscores)
# measures = data_summary(measures, varname="Measurement",
#                         groupnames=c("Method", 'PropRatio', "Measure"))
# 
# p1 = ggplot(measures, aes(PropRatio, Measurement, colour = Method)) +
#   geom_line(position=position_dodge(0.05)) +
#   geom_point(position=position_dodge(0.05)) +
#   geom_errorbar(aes(ymin=Measurement-sd, ymax=Measurement+sd), width=.2,
#                 position=position_dodge(0.05))+
#   facet_wrap(~Measure, scales = 'free_y') +
#   scale_color_brewer(palette = 'Paired')
# 
# p2 = ggplot(measures[measures$Measure %in% 'AUROC',], aes(PropRatio, Measurement, colour = Method)) +
#   geom_line() +
#   geom_point() +
#   scale_color_brewer(palette = 'Paired')
# 
# p3 = ggplot(measures[measures$Measure %in% 'AUROC',], aes(Method, Measurement, colour = Method)) +
#   geom_boxplot() +
#   scale_color_brewer(palette = 'Paired')