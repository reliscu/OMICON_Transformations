library(yaml)

updateConfig <- function() {
  
  default <- read_yaml("DefaultParamConfig.yaml")
  instance <- read_yaml("Snakemake/config/ParamConfig.yaml")
  instance <- instance[is.element(names(instance), names(default))]
  default[match(names(instance), names(default))] <- instance
  write_yaml(default, file = "Snakemake/config/ParamConfig.yaml")
  
}