source("renv/activate.R")
source("../helper_functions/load_renv.R")
source("FindModules.R")
      
allowWGCNAThreads(nThreads = snakemake@threads)

expr <- as.data.frame(fread(snakemake@input[[1]]))
sampleindexEnd <- snakemake@config[['sampleindexEnd']]
if(!is.numeric(sampleindexEnd)) {
  ## If not numeric, must be an R expression, i.e. ncol(expr)
  sampleindexEnd <- eval(parse(text=sampleindexEnd))
} 
samplegroups <- snakemake@config[['samplegroups']]
if(!is.null(samplegroups)) {
  ## If not NULL, must be name of .RDS file that contains a vector of length sampleindex
  samplegroups <- readRDS(samplegroups)
} 
subset <- snakemake@config[['subset']]
if(!is.null(subset)) {
  ## If not NULL, must be name of .RDS file that contains a logical vector of length nrow(expr)
  subset <- readRDS(subset)
} 
projectname <- snakemake@config[['projectname']]
geneinfo <- snakemake@config[['geneinfo']]
sampleindex <- snakemake@config[['sampleindexStart']]:sampleindexEnd
simMat <- snakemake@config[['simMat']]
saveSimMat <- snakemake@config[['saveSimMat']]
simType <- snakemake@config[['simType']]
overlapType <- snakemake@config[['overlapType']]
TOtype <- snakemake@config[['TOtype']]
TOdenom <- snakemake@config[['TOdenom']]
beta <- snakemake@config[['beta']]
MIestimator <- snakemake@config[['MIestimator']]
MIdisc <- snakemake@config[['MIdisc']]
signumType <- snakemake@config[['signumType']]
iterate <- snakemake@config[['iterate']]
signumvec <- snakemake@config[['signumvec']]
minsizevec <- snakemake@config[['minsizevec']]
merge.by <- snakemake@config[['merge.by']]
merge.param <- snakemake@config[['merge.param']]
export.merge.comp <- snakemake@config[['export.merge.comp']]
ZNCcut <-snakemake@config[['ZNCcut']]
calcSW <- snakemake@config[['calcSW']]
loadTree <- snakemake@config[['loadTree']]
writeKME <- snakemake@config[['writeKME']]
calcBigModStat <- snakemake@config[['calcBigModStat']]
writeModSnap <- snakemake@config[['writeModSnap']]

setwd("Snakemake/results")

FindModules(
  projectname,
  expr,
  geneinfo,
  sampleindex,
  samplegroups,
  subset,
  simMat,
  saveSimMat,
  simType,
  overlapType,
  TOtype,
  TOdenom,
  beta,
  MIestimator,
  MIdisc,
  signumType,
  iterate,
  signumvec,
  minsizevec,
  signum,
  minSize,
  merge.by,
  merge.param,
  export.merge.comp,
  ZNCcut,
  calcSW,
  loadTree,
  writeKME,
  calcBigModStat,
  writeModSnap
)