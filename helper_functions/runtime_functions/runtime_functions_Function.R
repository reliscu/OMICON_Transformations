.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")

library(data.table)

makeSessionInfo <- function(analyst, output) {
  
  date <- format(Sys.time(), '%B-%d-%Y')
  sif <- data.frame(Attribute = c("Analyst", "Date"), Value = c(as.character(analyst), date))
  write.csv(sif, file = output, row.names = F)
  
}

makeSessionInfo(analyst = snakemake@config['analyst'], output = snakemake@output[[1]])
