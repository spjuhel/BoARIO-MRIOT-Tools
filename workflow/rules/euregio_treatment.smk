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
    resources:
        mem_mb_per_cpu=6000,
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
        inp_file=expand(
            "{downloaded}/euregio/RegionalIOtable_{{year}}.xlsb",
            downloaded=config["downloaded_mriot_dir"],
        ),
    params:
        office_exists=OFFICE_EXISTS,
    resources:
        libre_office_instance=1,
    output:
        files=expand(
            "{downloaded}/euregio/euregio_{{year}}.csv",
            downloaded=config["downloaded_mriot_dir"],
        ),
    log:
        "logs/parse_euregio/convert_euregio_csvs_{year}.log",
    script:
        "../scripts/euregio_convert_csvs.py"

rule download_euregio:
    output:
        folder=directory(expand(
            "{downloaded}/euregio/", downloaded=config["downloaded_mriot_dir"]
        )),
        files=expand(
            "{downloaded}/euregio/{files}", downloaded=config["downloaded_mriot_dir"],
            files=["2000-2010-XLSB.zip"]
        )
    shell:
        """
        mkdir -p {output.folder}
        wget -O {output.files[0]} "https://dataportaal.pbl.nl/downloads/PBL_Euregio/PBL-EUREGIO-2000-2010-XLSB.zip"
        """


rule extract_euregio:
    input:
        folder=ancient(expand(
            "{downloaded}/euregio/", downloaded=config["downloaded_mriot_dir"]
        )),
        files=expand(
            "{downloaded}/euregio/{files}", downloaded=config["downloaded_mriot_dir"],
            files=["2000-2010-XLSB.zip"]
        )
    output:
        expand(
            "{downloaded}/euregio/RegionalIOtable_{years}.xlsb",
            downloaded=config["downloaded_mriot_dir"],
            years=[2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010],
        ),
    log:
        "logs/download_euregio/download_euregio.log",
    shell:
        """
        unzip -o {input.files[0]} -d {input.folder}
        """
