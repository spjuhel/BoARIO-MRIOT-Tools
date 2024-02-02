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
            "{aggregated}/icio2021/icio2021_2003_common_aggreg.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),


rule aggregate_eora26_common_test:
    input:
        expand(
            "{aggregated}/eora26/eora26_2016_common_aggreg.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),


rule aggregate_euregio_common_test:
    input:
        expand(
            "{aggregated}/euregio/euregio_2000_common_aggreg.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),


rule aggregate_exiobase3_common_test:
    input:
        expand(
            "{aggregated}/exiobase3_ixi/exiobase3_ixi_1995_common_aggreg.pkl",
            aggregated=config["aggregated_mriot_dir"],
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
            "{aggregated}/euregio/euregio_2000_common_aggreg.pkl",
            aggregated=config["aggregated_mriot_dir"],
        ),


def decide_aggreg_file(wildcards):
    return wildcards.aggregation == "common_aggreg"


rule aggregate_mriot:
    input:
        full_mriot_pkl=expand(
            "{parsed}/{{mriot_name}}/{{mriot_name}}_{{year}}_full.pkl",
            parsed=config["parsed_mriot_dir"],
        ),
        aggreg=branch(
            decide_aggreg_file,
            then=expand(
                "{aggregation_dir}/sectors_common_aggreg.ods",
                aggregation_dir=config["aggregation_csv_dir"],
            ),
            otherwise=expand(
                "{aggregation_dir}/{{mriot_name}}/{{mriot_name}}_{{aggregation}}.csv",
                aggregation_dir=config["aggregation_csv_dir"],
            ),
        ),
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
