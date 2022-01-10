rule session_info:
   output:
      sessionInfo = "Snakemake/results/{projectname}_SessionInfo.csv"
   script:
      "/home/rebecca/omicon/transformations/helper_functions/makeSessionInfo/makeSessionInfo.R"

rule param_config:
   output:
      paramConfig = "Snakemake/results/{projectname}_{transformation}/{projectname}_ParamConfig.yaml"
   shell:
      "cp Snakemake/config/ParamConfig.yaml Snakemake/results/{wildcards.projectname}_{transformation}/{wildcards.projectname}_ParamConfig.yaml"

