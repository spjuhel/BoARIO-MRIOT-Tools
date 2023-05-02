rule aggregate_exiobase3:
    input:
        mrio_pkl="pkls/exiobase3/exiobase3_{year}_full.pkl",
        aggreg="aggregation-files/exiobase3/exiobase3_{aggregation}.csv",
    params:
        full_mrio_params="params/exiobase3/exiobase3_full_params.csv",
    output:
        "pkls/exiobase3/exiobase3_{year}_{aggregation}.pkl",
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/aggregate_exiobase3/aggregate_exiobase3_{year}_{aggregation}.log",
    benchmark:
        "benchmarks/aggregate_exiobase3_{year}_{aggregation}.log"
    script:
        "../scripts/aggregate_exiobase3.py"


rule preparse_exiobase3:
    input:
        "autodownloads/exiobase3/IOT_{year}_ixi.zip",
    output:
        "pkls/exiobase3/exiobase3_{year}_full.pkl",
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/preparse_exiobase3/preparse_exiobase3_{year}.log",
    resources:
        mem_mb=6000,
    benchmark:
        "benchmarks/preparse_exiobase3_{year}.log"
    script:
        "../scripts/preparse_exiobase3.py"


rule download_exiobase3_test:
    input:
        "autodownloads/exiobase3/IOT_1995_ixi.zip",


rule download_exiobase3:
    output:
        "autodownloads/exiobase3/IOT_{year}_ixi.zip",
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/download_exiobase3/download_exiobase3_{year}.log",
    script:
        "../scripts/download_exiobase3.py"
