import os
import zipfile
import shutil

# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
# HTTPS = HTTPRemoteProvider()

OFFICE_BIN = shutil.which("libreoffice")
if OFFICE_BIN is not None:
    OFFICE_EXISTS = True
else:
    OFFICE_EXISTS = False

if shutil.which("unoserver") is not None:
    UNOSERVER_EXISTS = False
else:
    UNOSERVER_EXISTS = False


wildcard_constraints:
    aggregation=r"\d+_sectors|common_aggreg",


rule start_unoserver:
    params:
        office=OFFICE_BIN,
    log:
        "logs/unoserver.log",
    output:
        service("unoserver"),
    shell:
        "unoserver --executable {params.office} 2> {log};"


rule parse_euregio_test:
    input:
        expand(
            "{downloaded}/euregio/euregio_2000_full.pkl",
            downloaded=config["parsed_mriot_dir"],
        ),


rule parse_euregio:
    input:
        expand(
            "{downloaded}/euregio/euregio_{{year}}.csv",
            downloaded=config["downloaded_mriot_dir"],
        ),
    output:
        expand(
            "{parsed}/euregio/euregio_{{year}}_full.pkl",
            parsed=config["parsed_mriot_dir"],
        ),
    conda:
        "../envs/boario-tools-main.yml"
    log:
        "logs/parse_euregio/parse_euregio_{year}.log",
    script:
        "../scripts/parse_euregio.py"


rule create_euregio_csvs_test:
    input:
        expand(
            "{downloaded}/euregio/euregio_2000.csv",
            downloaded=config["downloaded_mriot_dir"],
        ),


rule create_euregio_csvs:
    input:
        expand(
            "{downloaded}/euregio/EURegionalIOtable_{{year}}.xlsx",
            downloaded=config["downloaded_mriot_dir"],
        ),
    output:
        expand(
            "{downloaded}/euregio/euregio_{{year}}.csv",
            downloaded=config["downloaded_mriot_dir"],
        ),
    shell:
        """
        xlsx2csv -s 3 {input} {output}
        """


rule create_euregio_xlsx_test:
    input:
        expand(
            "{downloaded}/euregio/EURegionalIOtable_2000.xlsx",
            downloaded=config["downloaded_mriot_dir"],
        ),


rule create_euregio_xlsx:
    input:
        inp_file=expand(
            "{downloaded}/euregio/EURegionalIOtable_{{year}}.ods",
            downloaded=config["downloaded_mriot_dir"],
        ),
    params:
        office_exists=OFFICE_EXISTS,
        uno_exists=UNOSERVER_EXISTS,
    resources:
        libre_office_instance=1,
    output:
        files=expand(
            "{downloaded}/euregio/EURegionalIOtable_{{year}}.xlsx",
            downloaded=config["downloaded_mriot_dir"],
        ),
    log:
        "logs/parse_euregio/convert_euregio_xlsx_{year}.log",
    script:
        "../scripts/euregio_convert_xlsx.py"

rule download_euregio:
    output:
        folder=directory(expand(
            "{downloaded}/euregio/", downloaded=config["downloaded_mriot_dir"]
        )),
        files=expand(
            "{downloaded}/euregio/{files}", downloaded=config["downloaded_mriot_dir"],
            files=["2000-2005-ODS.zip","2006-2010-ODS.zip"]
        )
    shell:
        """
        mkdir -p {output.folder}
        wget -O {output.files[0]} "https://dataportaal.pbl.nl/downloads/PBL_Euregio/PBL-EUREGIO-2000-2005-ODS.zip"
        wget -O {output.file[1]} "https://dataportaal.pbl.nl/downloads/PBL_Euregio/PBL-EUREGIO-2006-2010-ODS.zip"
        """


rule extract_euregio:
    input:
        folder=expand(
            "{downloaded}/euregio/", downloaded=config["downloaded_mriot_dir"]
        ),
        files=expand(
            "{downloaded}/euregio/{files}", downloaded=config["downloaded_mriot_dir"],
            files=["2000-2005-ODS.zip","2006-2010-ODS.zip"]
        )
    output:
        expand(
            "{downloaded}/euregio/EURegionalIOtable_{years}.ods",
            downloaded=config["downloaded_mriot_dir"],
            years=[2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010],
        ),
    log:
        "logs/download_euregio/download_euregio.log",
    shell:
        """
        unzip -o {input.files[0]} -d {input.folder}
        unzip -o {input.files[1]} -d {input.folder}
        """
