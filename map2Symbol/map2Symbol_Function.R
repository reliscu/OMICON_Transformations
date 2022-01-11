library(dplyr)
library(data.table)

map2Symbol <- function(expr, mappingTable, output) {
  
  expr <- as.data.frame(fread(expr))
  mappingTable <- read.csv(mappingTable)

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
  mapped_expr <- data.frame(mappedIds, expr[,2:ncol(expr)])
  fwrite(mapped_expr, file = output)
  
}
