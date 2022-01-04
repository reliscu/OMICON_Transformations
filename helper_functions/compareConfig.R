library(yaml)

compareConfig <- function() {
  
  default <- read_yaml("../../DefaultParamConfig.yaml")
  instance <- read_yaml("../config/ParamConfig.yaml")
  instance <- instance[is.element(names(instance), names(default))]
  default[match(names(instance), names(default))] <- instance
  write_yaml(default, file = "../config/ParamConfig.yaml")
  
}