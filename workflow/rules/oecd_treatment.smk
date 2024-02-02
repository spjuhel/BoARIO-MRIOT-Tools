rule parse_icio2021:
    input:
        expand("{downloaded}/icio2021/ICIO2021_{{year}}.csv",downloaded=config['downloaded_mriot_dir']),
    output:
        expand("{parsed}/icio2021/icio2021_{{year}}_full.pkl",parsed=config["parsed_mriot_dir"]),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/parse_icio2021/parse_icio2021_{year}.log",
    script:
        "../scripts/parse_icio2021.py"


rule download_icio2021:
    output:
        expand("{downloaded}/icio2021/ICIO2021_{{year}}.csv",downloaded=config['downloaded_mriot_dir'])
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/download_icio2021/download_icio2021_{year}.log",
    script:
        "../scripts/download_icio2021.py"
