load_renv <- function(lock.path) {
  
  Sys.setenv(RENV_PATHS_LIBRARY = "/home/rebecca/omicon/transformations/renv_lib")
  source("/home/rebecca/omicon/transformations/helper_functions/renv_functions/activate.R")
  options(renv.consent=T)
  renv::restore(confirm=F,lockfile=lock.path)
  
}
