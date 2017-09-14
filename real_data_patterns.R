#patterns in TCGA
#----set path and other parameters for the analysis----
datapath = '/wehisan/home/allstaff/b/bhuva.d/data'
source('/wehisan/home/allstaff/b/bhuva.d/conditional_interactions/code/functions/generate_conditions.R')

#----load required libraries----
liblist = c('Biobase', 'mclust', 'doParallel', 'parallel', 'foreach', 'plyr', 'iterators')
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

#generate bimodal genes list
conditions = getBimodalGenes(exprmat,Hellinger.thresh = 0.28, WVRS.thresh = exp(-0.35))
modelinfo = conditions[,1:7]
conditions = conditions[,-(1:7)]
stopCluster(cl)

targetlist = rownames(exprmat)
tflist = nursaNHRs$Symbol[nursaNHRs$Symbol %in% rownames(exprmat)]

itcoreg = iter(conditions[1, , drop = F], by = 'row')
fits <- foreach(coreg = itcoreg, .combine = rbind) %dopar% {
  cond1 =  !is.na(coreg) & coreg == 0
  cond2 =  !is.na(coreg) & coreg == 1
  res1 = fitModels(t(exprmat[tflist, cond1]), t(exprmat[targetlist, cond1]))
  res2 = fitModels(t(exprmat[tflist, cond2]), t(exprmat[targetlist, cond2]))
  colnames(res1)[3:5] = paste0(colnames(res1)[3:5], '.1')
  colnames(res2)[3:5] = paste0(colnames(res2)[3:5], '.2')
  res = cbind(res1, res2[, 3:5])
  
  return(res)
}




