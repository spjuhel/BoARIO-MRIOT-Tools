wildcard_constraints:
    mriot_name=r"eora26|euregio|icio2021",


rule sector_config_from_exio3:
    input:
        exio3_sectors_config=expand("{mriot_params_dir}/exiobase3_ixi_full_sectors.csv",mriot_params_dir=config['mriot_params_dir']),
        aggregation_master=expand("{mriot_aggreg_dir}/exiobase3_ixi/exiobase3_ixi_to_other_mrio_sectors.ods",mriot_aggreg_dir=config['aggregation_csv_dir']),
    output:
        exio3_sectors_config=expand("{mriot_params_dir}/{{mriot_name}}_full_sectors.csv",mriot_params_dir=config['mriot_params_dir'])
    conda:
        "../envs/boario-tools-main.yml"
    params:
        mrio_type="{mriot_name}",
    script:
        "../scripts/params_gen_from_exiobase3_full.py"
