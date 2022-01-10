.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(data.table)

map2Symbol <- function(expr, mappingTable, output) {
  
  expr <- as.data.frame(fread(as.character(expr)))
  mappingTable <- read.csv(as.character(mappingTable))

  nrow1 <- nrow(mappingTable)
  
  mappingTable <- mappingTable %>%
    na_if("") %>%
    group_by(!!as.name(colnames(mappingTable)[1])) %>%
    slice(1) %>%
    as.data.frame()
  
  nrow2 <- nrow(mappingTable)
  
  if(nrow2<nrow1) {
    warning("Mapping table contains many-to-one (key:target) mappings. The gene symbol for these keys was arbitrarily selected")
  }
  
  mappedIds <- merge(expr[,1], mappingTable, by = 1)
  mappedIds[is.na(mappedIds[,2]),2] <- mappedIds[is.na(mappedIds[,2]),1]
  
  colnames(mappedIds)[2] <- "Gene"
  expr <- data.frame(mappedIds, expr[,2:ncol(expr)])
  fwrite(expr, file = as.character(output))
  
}

map2Symbol(snakemake@input['expr'], snakemake@input['mappingTable'], snakemake@output[1])












# colnames(expr)[2:ncol(expr)] <- expr[1,2:ncol(expr)]
# expr <- expr[-c(1), -c(1)]
# expr <- data.frame(expr[,1], as.matrix(expr[,2:ncol(expr)]))