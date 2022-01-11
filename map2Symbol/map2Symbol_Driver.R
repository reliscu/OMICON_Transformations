source("../helper_functions/renv_functions/load_renv.R")
load_renv(lock.path = "../map2Symbol/renv.lock")

source("../map2Symbol/map2Symbol_Function.R")
map2Symbol(snakemake@input[['expr']], snakemake@input[['mappingTable']], snakemake@output[[1]])
