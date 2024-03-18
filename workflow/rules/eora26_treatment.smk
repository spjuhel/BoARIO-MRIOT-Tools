rule parse_eora26:
    input:
        expand(
            "{downloaded}/eora26/Eora26_{{year}}_bp.zip",
            downloaded=config["downloaded_mriot_dir"],
        ),
    output:
        expand(
            "{parsed}/eora26/eora26_{{year}}_full.pkl",
            parsed=config["parsed_mriot_dir"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/parse_eora26/parse_eora26_{year}.log",
    resources:
        mem_mb_per_cpu=6000,
    benchmark:
        "benchmarks/mrios/parse_eora26_{year}.log"
    script:
        "../scripts/parse_eora26.py"
