include: "exio3_treatment.smk"
include: "euregio_treatment.smk"
include: "eora26_treatment.smk"
include: "oecd_treatment.smk"
include: "sectors_config.smk"


ALL_FULL_MRIOT = expand(
    "pkls/{mrio_type}/{mrio_type}_{year}_full.pkl",
    mrio_type=["euregio", "eora26", "oecd_v2021", "exiobase3_ixi"],
    year=[2000, 2010],
)


rule all_mriot:
    input:
        ALL_FULL_MRIOT,


rule aggregate_icio2021_common_test:
    input:
        expand(
            "{aggregated}/icio2021/icio2021_2003_{common_aggreg}.pkl",
            aggregated=config["aggregated_mriot_dir"],
            common_aggreg=config["common_aggreg"],
        ),


rule aggregate_eora26_common_test:
    input:
        expand(
            "{aggregated}/eora26/eora26_2010_{common_aggreg}.pkl",
            aggregated=config["aggregated_mriot_dir"],
            common_aggreg=config["common_aggreg"],
        ),


rule aggregate_euregio_common_test:
    input:
        expand(
            "{aggregated}/euregio/euregio_2000_{common_aggreg}.pkl",
            aggregated=config["aggregated_mriot_dir"],
            common_aggreg=config["common_aggreg"],
        ),


rule aggregate_exiobase3_common_test:
    input:
        expand(
            "{aggregated}/exiobase3_ixi/exiobase3_ixi_1995_{common_aggreg}.pkl",
            aggregated=config["aggregated_mriot_dir"],
            common_aggreg=config["common_aggreg"],
        ),


rule aggregate_exiobase3_test:
    input:
        expand(
            "{aggregated}/exiobase3_ixi/exiobase3_ixi_1995_74_sectors.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),


rule aggregate_euregio_test:
    input:
        expand(
            "{aggregated}/euregio/euregio_2000_{common_aggreg}.pkl",
            aggregated=config["aggregated_mriot_dir"],
            common_aggreg=config["common_aggreg"],
        ),


def get_full_mriot(wildcards):
    full_agg_name = config["mriot_base_aggreg"][wildcards.mriot_name]
    parsed = config["parsed_mriot_dir"]
    mriot_name = wildcards.mriot_name
    year = wildcards.year
    return f"{parsed}/{mriot_name}/{mriot_name}_{year}_{full_agg_name}.pkl"


rule aggregate_mriot:
    input:
        full_mriot_pkl=get_full_mriot,
    output:
        expand(
            "{aggregated}/{{mriot_name}}/{{mriot_name}}_{{year}}_{{aggregation}}.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/aggregate_{mriot_name}/aggregate_{mriot_name}_{year}_{aggregation}.log",
    benchmark:
        "benchmarks/aggregate_{mriot_name}_{year}_{aggregation}.log"
    script:
        "../scripts/aggregate_mriot.py"
