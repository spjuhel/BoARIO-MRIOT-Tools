ruleorder: full_sector_config_from_exio3 > aggreg_sector_config

base_aggreg ={
    "exiobase3_ixi":"full_sectors",
    "euregio":"full_sectors",
    "eora26":"full_no_reexport_sectors"
}


rule all_mriots_full_sectors_config:
    input:
        sectors_config=expand(
            "{mriot_params_dir}/{mriot_name}_full_sectors.csv",
            mriot_params_dir=config["mriot_params_dir"],
            mriot_name=["eora26", "euregio", "icio2021", "exiobase3_ixi"],
        ),

rule full_sector_config_from_exio3:
    input:
        exio3_sectors_config=expand(
            "{mriot_params_dir}/exiobase3_ixi_full_sectors.csv",
            mriot_params_dir=config["mriot_params_dir"],
        ),
    params:
        alt_aggregation_master = None,
        mrio_dict = {
            "eora26_full_no_reexport_sectors": "Eora sectors",
            "euregio_full_sectors": "Euregio sectors",
            "icio2021_full_reworked_sectors": "ICIO2021_reworked sectors name",
        }
    log:
        "logs/{mriot_name_noexio_full_sectors}_config_from_exio3.log",
    wildcard_constraints:
        mriot_name_noexio_full_sectors = r"icio2021_full_reworked_sectors|eora26_full_no_reexport_sectors|euregio_full_sectors"
    output:
        sectors_config=expand(
            "{mriot_params_dir}/{{mriot_name_noexio_full_sectors}}.csv",
            mriot_params_dir=config["mriot_params_dir"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    script:
        "../scripts/params_gen_from_exiobase3_full.py"

def get_base_sectors_config(wildcards):
    return expand(
        "{mriot_params_dir}/{name}.csv",
        mriot_params_dir=config["mriot_params_dir"],
        name = wildcards.mriot_name+"_"+base_aggreg[wildcards.mriot_name]
    )

rule aggreg_sector_config:
    input:
        sectors_config=get_base_sectors_config
    output:
        agg_sectors_config=expand(
            "{mriot_params_dir}/{{mriot_name}}_{{sectors_aggregation}}.csv",
            mriot_params_dir=config["mriot_params_dir"],
        ),
    params:
        alt_aggregation_file = None,
        base_aggreg = base_aggreg
    log:
        "logs/{mriot_name}_{sectors_aggregation}_config.log",
    conda:
        "../envs/boario-tools-main.yml"
    script:
        "../scripts/sectors_params_aggreg.py"
