# General Configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

## Output folders

By default, this workflow stores the files obtained, parsed and aggregated into three different folders (`downloaded/`, `parsed/` `aggregated/`) in the working directory. This behavior can be changed with the entries `downloaded_mriot_dir:`, `parsed_mriot_dir:` and `aggregated_mriot_dir:`.

Files corresponding to a specific MRIOT will be stored in a subfolder of the corresponding name. For instance, Exiobase 3 (in ixi system) for year 1995 aggregated to 74 sectors will be stored under `parsed/exiobase3_ixi/exiobase3_ixi_1995_74_sectors.pkl`.

## Sector aggregation files

The workflow looks for aggregation files in the `aggregation_csv_dir:` (`config/aggregation_files`) by default.

Some existing aggregation are included in the repository, and new ones can be added. If you think an aggregation file you designed would be a valuable contribution, don't hesitate to pull-request it !

The folder structure follows the same as for outputs, for instance, the file used to aggregate Exiobase 3 to 74 sectors is looked under `config/aggregation_files/exiobase3_ixi/exiobase3_ixi_74_sectors.csv`

Files should be in csv format, with a first column `sector` with all the original sectors of the MRIOT, and `group` and `name` columns with the new ID and names to map the original sectors too.

## Common aggregation

As an additional way to aggregate different MRIOTs into a common one, this workflow also embarks a `config/aggregation_files/sectors_common_aggreg.ods` file for this purpose.

The sheet `common_aggreg` of the file defines the sector to aggregate to via `group_id` and `sector_name` columns. The following sheets, `<mriot_name>_to_common_aggreg`, have `original sector`, `to_id` and a `new sector` columns auto computed from `to_id`.

Asking the workflow to build `aggregated/icio2021/icio2021_2003_common_aggreg.pkl`, for instance, will aggregate `icio2021_2003_full.pkl` to the sectors defined in the sheet `common_aggreg` based on the correspondence in `icio2021_to_common_aggreg`.

Disclaimer: This work is still under progress and still requires to be streamlined.

## Parameters aggregation

The workflow also allows to build sector parameters files for the `BoARIO` model.
We provide a `config/mriot_params/exiobase3_ixi_full_sectors.csv` "master" file that contains values we used with the model and this MRIOT, and from which similar files can be obtained for the other MRIOT, using `config/aggregation_files/exiobase3_ixi/exiobase3_ixi_to_other_mrio_sectors.ods`.
