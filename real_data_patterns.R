#patterns in TCGA
#----set path and other parameters for the analysis----
datapath = '/wehisan/home/allstaff/b/bhuva.d/data'
source('/wehisan/home/allstaff/b/bhuva.d/conditional_interactions/code/functions/generate_conditions.R')

#----load required libraries----
liblist = c('Biobase', 'mclust', 'doParallel', 'parallel', 'foreach', 'plyr',
            'iterators', 'ggplot2', 'grid', 'viridis', 'reshape2')
for (lib in liblist) {
  library(lib, character.only = T)
}
remove(lib, liblist)

#----load data----
dataFile = file.path(datapath, 'easy_load', 'exprmat_TCGA_seq.RData')
tfFile = file.path(datapath, 'databases', 'NURSA_nhrs.csv')
cofFile = file.path(datapath, 'databases', 'NURSA_cofactor_aliases.csv')
load(dataFile)
nursaNHRs = read.csv(tfFile, colClasses = rep('character', 2))
nursaCofs = read.csv(cofFile,
                     header = T,
                     colClasses = c('character'))$cof_aliases
nursaCofs = subset(nursaCofs, nursaCofs %in% rownames(exprmat))
remove(dataFile, cofFile, tfFile)

#----fit models and store parameters----
#x and y have variables along the columns and samples on the rows
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

cl = makeCluster(5, outfile = '')
registerDoParallel(cl)

#----generate bimodal genes list----
conditions = getBimodalGenes(exprmat,Hellinger.thresh = 0.28, WVRS.thresh = exp(-0.35))
modelinfo = conditions[,1:7]
conditions = conditions[,-(1:7)]
stopCluster(cl)

cl = makeCluster(5, outfile = '')
registerDoParallel(cl)

#----Known interactions----
targetlist = rownames(exprmat)[rank(-apply(exprmat, 1, var))<5000]
# tflist = nursaNHRs$Symbol[nursaNHRs$Symbol %in% rownames(exprmat)]
tflist = c('ESR1', 'AR', 'PGR', 'VDR')
conds = c(nursaCofs[nursaCofs %in% rownames(conditions)], 'FOXA1', 'GATA3')

itcoreg = iter(conditions[conds, , drop = F], by = 'row')
fitsKnown <- foreach(coreg = itcoreg, .combine = rbind) %dopar% {
  cond1 =  !is.na(coreg) & coreg == 0
  cond2 =  !is.na(coreg) & coreg == 1
  res1 = fitModels(t(exprmat[tflist, cond1, drop = F]), t(exprmat[targetlist, cond1, drop = F]))
  res2 = fitModels(t(exprmat[tflist, cond2, drop = F]), t(exprmat[targetlist, cond2, drop = F]))
  colnames(res1)[3:5] = paste0(colnames(res1)[3:5], '.1')
  colnames(res2)[3:5] = paste0(colnames(res2)[3:5], '.2')
  res = cbind('cond' = rownames(coreg),res1, res2[, 3:5])
  
  return(res)
}
save(fitsKnown, file = 'simdata/lmFits.RData')

#----Random interactions----
set.seed(360)
targetlist = rownames(exprmat)[rank(-apply(exprmat, 1, var))<5000]
tflist = sample(rownames(exprmat), 100)
conds = sample(rownames(conditions), length(conds))

itcoreg = iter(conditions[conds, , drop = F], by = 'row')
fitsRandom <- foreach(coreg = itcoreg, .combine = rbind) %dopar% {
  cond1 =  !is.na(coreg) & coreg == 0
  cond2 =  !is.na(coreg) & coreg == 1
  res1 = fitModels(t(exprmat[tflist, cond1, drop = F]), t(exprmat[targetlist, cond1, drop = F]))
  res2 = fitModels(t(exprmat[tflist, cond2, drop = F]), t(exprmat[targetlist, cond2, drop = F]))
  colnames(res1)[3:5] = paste0(colnames(res1)[3:5], '.1')
  colnames(res2)[3:5] = paste0(colnames(res2)[3:5], '.2')
  res = cbind('cond' = rownames(coreg),res1, res2[, 3:5])
  
  return(res)
}
save(fitsKnown, fitsRandom, file = 'simdata/lmFits.RData')

stopCluster(cl)

#----plots and statistics----
fitsRandom = fitsRandom[sample.int(nrow(fitsRandom), 0.5E6), ]
# fitsKnown = fitsKnown[sample.int(nrow(fitsKnown), 0.5E6), ]

fts1 = fitsRandom[, 1:6]
fts2 = fitsRandom[, c(1:3, 7:9)]
colnames(fts1)[4:6] = colnames(fts2)[4:6] = c('m', 'c', 'p')
fts = melt(fts1, id.vars = c('cond', 'x', 'y'))
fts$condition2 = melt(fts2, id.vars = c('cond', 'x', 'y'))$value
colnames(fts)[5] = 'condition1'
ftsRandom = fts

fts1 = fitsKnown[, 1:6]
fts2 = fitsKnown[, c(1:3, 7:9)]
colnames(fts1)[4:6] = colnames(fts2)[4:6] = c('m', 'c', 'p')
fts = melt(fts1, id.vars = c('cond', 'x', 'y'))
fts$condition2 = melt(fts2, id.vars = c('cond', 'x', 'y'))$value
colnames(fts)[5] = 'condition1'
ftsKnown = fts
remove(fts1, fts2, fts)

textSize = 1.5
cof = 'FOXA1'#HSDL1

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plotdf = ftsKnown[ftsKnown$cond %in% cof,]
plotdf = plotdf[as.character(plotdf$x) != as.character(plotdf$y),]
datadens = dlply(plotdf, 'variable', function(x) get_density(x$condition1, x$condition2))
plotdf$density = 0
plotdf$density[plotdf$variable %in% 'm'] = datadens$m
plotdf$density[plotdf$variable %in% 'c'] = datadens$c
plotdf$density[plotdf$variable %in% 'p'] = datadens$p

#sample points to plot: not plotting all points
probp = plotdf$density[plotdf$variable %in% 'p']
probp = 1/(probp + 1)
probp = probp/sum(probp)
probm = plotdf$density[plotdf$variable %in% 'm']
probm = 1/(probm + 1)
probm = probm/sum(probm)
probc = plotdf$density[plotdf$variable %in% 'c']
probc = 1/(probc + 1)
probc = probc/sum(probc)
nsamp=5000

#plot for coefficient of correlation
pp = ggplot(plotdf[plotdf$variable %in% 'p', ][sample(1:length(probp), nsamp, prob = probp),]) +
  geom_point(aes(condition1, condition2, color = density), size = 1) +
  scale_colour_distiller(palette = 'PRGn', direction = 1) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
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

#plot for slope
pm = ggplot(plotdf[plotdf$variable %in% 'm', ][sample(1:length(probm), nsamp, prob = probm),]) +
  geom_point(aes(condition1, condition2, color = density), size = 1) +
  scale_colour_distiller(palette = 'PRGn', direction = 1) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
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
  )
pm = pm + geom_density2d(data = ftsRandom[ftsRandom$variable %in% 'm', ], aes(condition1, condition2))

#plot for intercept
pc = ggplot(plotdf[plotdf$variable %in% 'c', ][sample(1:length(probc), nsamp, prob = probc),]) +
  geom_point(aes(condition1, condition2, color = density), size = 1) +
  scale_colour_distiller(palette = 'PRGn', direction = 1) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
  geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
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
  )
pc = pc + geom_density2d(data = ftsRandom[ftsRandom$variable %in% 'c', ], aes(condition1, condition2))

ptitle = ggplot() + ggtitle(paste0('Conditional interaction characteristics of gene: ', cof)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = rel(textSize+0.5),
      hjust = 0.5
    )
  )

gridsz = 10
layout = c(rep(c(1, rep(2, gridsz)), gridsz),
           rep(c(1, rep(3, gridsz)), gridsz),
           rep(c(1, rep(4, gridsz)), gridsz))
layout = matrix(layout, ncol = 3 * gridsz)

svg(filename = 'figures/fig 5 - real_patterns_FOXA1.svg', width = 12, height = 6)
multiplot(ptitle, pm, pc, pp, layout = layout)
dev.off()




# p = ggplot(plotdf, aes(condition1, condition2)) + stat_bin2d(aes(fill = (..count..)), bins = 200, geom="tile") +
#   scale_fill_distiller(palette = 'PRGn', direction = 1) +
#   facet_wrap(~variable, scales = 'free', labeller = label_both, ncol = 3) +
#   geom_hline(yintercept = 0, size = 0.5, linetype = 'dotdash') +
#   geom_vline(xintercept = 0, size = 0.5, linetype = 'dotdash') +
#   ggtitle(paste0('Coregulatory behaviour of: ', cof)) +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = rel(textSize)),
#     axis.text.x = element_text(angle = 0, size = rel(textSize)),
#     axis.text.y = element_text(angle = 0, size = rel(textSize)),
#     strip.background = element_rect(colour = "gray80", fill = "gray80"),
#     strip.text = element_text(size = rel(textSize)),
#     axis.line = element_line(colour = "black"),
#     axis.ticks = element_line(),
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.margin = margin(unit(0, "cm")),
#     legend.title = element_text(face = "italic"),
#     plot.title = element_text(
#       face = "bold",
#       size = rel(textSize),
#       hjust = 0.5
#     )
#   )
# 
# p = p + geom_density2d(data = ftsRandom, aes(condition1, condition2)) +
#   scale_colour_viridis(option = 'viridis') +
#   facet_wrap(~variable, scales = 'free', labeller = label_both, ncol = 3) +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = rel(textSize)),
#     axis.text.x = element_text(angle = 0, size = rel(textSize)),
#     axis.text.y = element_text(angle = 0, size = rel(textSize)),
#     strip.background = element_rect(colour = "gray80", fill = "gray80"),
#     strip.text = element_text(size = rel(textSize)),
#     axis.line = element_line(colour = "black"),
#     axis.ticks = element_line(),
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.margin = margin(unit(0, "cm")),
#     legend.title = element_text(face = "italic"),
#     plot.title = element_text(
#       face = "bold",
#       size = rel(textSize),
#       hjust = 0.5
#     )
#   )
# p







