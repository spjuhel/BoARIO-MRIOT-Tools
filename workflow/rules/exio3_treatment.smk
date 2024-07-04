wildcard_constraints:
    system=r"ixi|pxp",


rule parse_exiobase3_test:
    input:
        expand(
            "{parsed}/exiobase3_ixi/exiobase3_ixi_1995_{base_aggreg}.pkl",
            parsed=config["parsed_mriot_dir"],
            base_aggreg=config["mriot_base_aggreg"]["exiobase3_ixi"],
        ),


rule parse_exiobase3:
    input:
        expand(
            "{downloaded}/exiobase3_{{system}}/IOT_{{year}}_{{system}}.zip",
            downloaded=config["downloaded_mriot_dir"],
        ),
    output:
        expand(
            "{parsed}/exiobase3_{{system}}/exiobase3_{{system}}_{{year}}_{base_aggreg}.pkl",
            parsed=config["parsed_mriot_dir"],
            base_aggreg=config["mriot_base_aggreg"]["exiobase3_ixi"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/parse_exiobase3/parse_exiobase3_{system}_{year}.log",
    resources:
        mem_mb_per_cpu=10000,
        disk_mb=3000,
        disk_mib=3000,
    benchmark:
        "benchmarks/parse_exiobase3_{year}_{system}.log"
    script:
        "../scripts/parse_exiobase3.py"


rule download_exiobase3_test:
    input:
        expand(
            "{downloaded}/exiobase3_{{system}}/IOT_1995_ixi.zip",
            downloaded=config["downloaded_mriot_dir"],
        ),


rule download_exiobase3:
    output:
        expand(
            "{downloaded}/exiobase3_{{system}}/IOT_{{year}}_{{system}}.zip",
            downloaded=config["downloaded_mriot_dir"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/download_exiobase3/download_exiobase3_{year}_{system}.log",
    script:
        "../scripts/download_exiobase3.py"
