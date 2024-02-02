# Snakemake workflow: `BoARIO-MRIOT-Tools`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for downloading, parsing and aggregating Multi-Regional Input Output Tables (MRIOTs), using [`pymrio`](https://pymrio.readthedocs.io/en/latest/), as well as associated configuration files for the purpose of being used in simulations with the `BoARIO` model [package](https://github.com/spjuhel/BoARIO).

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=spjuhel%2FBoARIO-MRIOT-Tools).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

### File naming structure

Generated `pymrio.IOSystem` files have the following format:

`{mriot_name}/{mriot_name}_{year}_{aggregation}.pkl`

Where:
 - `mriot_name` can be any of `exiobase3_ixi`, `eora26`, `euregio`, `icio2021`
 - `year` can be any of the existing year for the corresponding MRIOT
 - `aggregation` can be `common_aggreg`, or any value for which a aggregation correspondence file exists (`config/aggregation_files/{mriot_name}_{aggretation}`.

### Recommendations

Here are some recommendations:

- Keep track of the MRIOT sources you use and don't forget to cite the associated project.
- Monitor your disk usage, MRIOT consume quite the space.
- Don't hesitate to raise an issue !

## EORA 26 download

Accessing Eora26 MRIOT requires an account, thus automatic download is less easy. The current way for the workflow to work is to download manually[^1]:

- Connect to your account on https://worldmrio.com/eora26/
- Download the MRIOT you want (`Eora26_2000_bp.zip` for instance) and put it inside the download folder (`downloaded/eora26/` per default).
- Also download `indice.zip` from the same page, into the same folder.

[^1]: `pymrio` now allows to use credentials for this, which will be added in a future version of the pipeline.
