.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")

library(yaml)

compareConfig <- function() {
  
  default <- read_yaml("DefaultParamConfig.yaml")
  instance <- read_yaml("Snakemake/config/ParamConfig.yaml")
  instance <- instance[is.element(names(instance), names(default))]
  default[match(names(instance), names(default))] <- instance
  write_yaml(default, file = "Snakemake/config/ParamConfig.yaml")
  
}