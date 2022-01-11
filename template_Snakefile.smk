transformation = ""

workdir: "/home/rebecca/omicon/transformations/" + transformation
configfile: "Snakemake/config/ParamConfig.yaml"
include: "/home/rebecca/omicon/transformations/helper_functions/helper.smk" 
shell("Rscript /home/rebecca/omicon/transformations/helper_functions/updateConfig/updateConfig_Driver.R")

rule all:
   input:
      expand(["Snakemake/results/{projectname}_{transformation}/{projectname}_SessionInfo.csv", 
              "Snakemake/results/{projectname}_{transformation}/{projectname}_ParamConfig.yaml"], 
              projectname = config['projectname'], 
              transformation = transformation)
