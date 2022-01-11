makeSessionInfo <- function(analyst, output) {

  date <- format(Sys.time(), '%B-%d-%Y')
  sif <- data.frame(Attribute = c("Analyst", "Date"), Value = c(as.character(analyst), date))
  write.csv(sif, file = output, row.names = F)
  
}
