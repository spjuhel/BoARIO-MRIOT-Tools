import pandas as pd

# rule aggregate_eora26:
#     input:
#         mrio_pkl="mrio-files/pkls/eora26/eora26_{year}_full.pkl",
#         aggreg="aggregation-files/eora26/eora26_{aggregation}.csv"
#     params:
#         full_mrio_params = "mrio-files/params/eora26/eora26_full_params.csv"
#     output:
#         "mrio-files/pkls/eora26/eora26_{year}_{aggregation}.pkl"
#     conda:
#         "../envs/boario-tools-main.yml"
#     log:
#         "logs/aggregate_eora26/aggregate_eora26_{year}_{aggregation}.log"
#     script:
#         "../scripts/aggregate_eora26.py"


rule eora_sector_config:
    input:
        exio3_sectors_config="config/exiobase3_full_sectors.csv",
        aggregation_master="config/exiobase3_to_other_mrio_sectors.ods",
    output:
        "config/eora26_full_sectors.csv",
    conda:
        "../envs/boario-tools-main.yml"
    params:
        mrio_type="eora26",
    script:
        "../scripts/params_gen_from_exiobase3_full.py"


rule preparse_eora26:
    input:
        "autodownloads/eora26/Eora26_{year}_bp.zip",
    output:
        "mrio-files/pkls/eora26/eora26_{year}_full.pkl",
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/preparse_eora26/preparse_eora26_{year}.log",
    resources:
        mem_mb=6000,
    benchmark:
        "benchmarks/mrios/preparse_eora26_{year}.log"
    script:
        "../scripts/preparse_eora26.py"
