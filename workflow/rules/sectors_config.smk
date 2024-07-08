ruleorder: full_sector_config_from_exio3 > aggreg_sector_config

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
        "logs/full_configs_from_exio3.log",
    output:
        sectors_config=expand(
            "{mriot_params_dir}/{config}.csv",
            mriot_params_dir=config["mriot_params_dir"],
            config=["icio2021_full_reworked_sectors",
                     "eora26_full_no_reexport_sectors",
                     "euregio_full_sectors"]
        ),
    conda:
        "../envs/boario-tools-main.yml"
    script:
        "../scripts/params_gen_from_exiobase3_full.py"

def is_eora26(wildcards):
    return ("eora26" in wildcards.mriot_name)

def is_icio2021(wildcards):
    return ("icio2021" in wildcards.mriot_name)

rule aggreg_sector_config:
    input:
        sectors_config = branch(
                condition=is_eora26,
                then=f"{config['mriot_params_dir']}/eora26_full_no_reexport_sectors.csv",
                otherwise=branch(
                    condition=is_icio2021,
                    then=f"{config['mriot_params_dir']}/icio2021_full_reworked_sectors.csv",
                    otherwise=f"{config['mriot_params_dir']}"+"/{mriot_name}_full_sectors.csv"
                )
            )
    output:
        agg_sectors_config=expand(
            "{mriot_params_dir}/{{mriot_name}}_{{sectors_aggregation_nofull}}.csv",
            mriot_params_dir=config["mriot_params_dir"],
        ),
    params:
        alt_aggregation_file = None,
    log:
        "logs/{mriot_name}_{sectors_aggregation_nofull}_config.log",
    conda:
        "../envs/boario-tools-main.yml"
    script:
        "../scripts/sectors_params_aggreg.py"
