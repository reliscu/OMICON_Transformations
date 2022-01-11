source("../helper_functions/makeSessionInfo/makeSessionInfo_Function.R")

makeSessionInfo(analyst = snakemake@config['analyst'], output = snakemake@output[[1]])
