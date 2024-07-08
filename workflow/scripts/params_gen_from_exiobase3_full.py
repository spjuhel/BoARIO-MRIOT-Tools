from importlib import resources
from pathlib import Path
import pandas as pd
import logging
import sys, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

sectors_df = pd.read_csv(
    snakemake.input.exio3_sectors_config[0], index_col=0, decimal="."
)
if snakemake.params["alt_aggregation_master"] is not None:
    aggregation_master_df = pd.read_excel(
        snakemake.params.alt_aggregation_master, sheet_name=0, index_col=0
    )
else:
    with resources.path(
        "boario_tools.data.aggregation_files.exiobase3_ixi",
        "exiobase3_ixi_to_other_mrio_sectors.ods",
    ) as agg_path:
        aggregation_master_df = pd.read_excel(
            agg_path, sheet_name=0, index_col=0
        )

for mrio, colname in snakemake.params["mrio_dict"].items():
    res = (
        aggregation_master_df.join(sectors_df)[
            [
                colname,
                "affected",
                "rebuilding_factor",
                "inventory_size",
                "productive_capital_to_va_ratio",
                "inventory_tau",
            ]
        ]
        .groupby(colname)
        .agg(
            {
                "affected": any,
                "rebuilding_factor": "sum",
                "inventory_size": "mean",
                "productive_capital_to_va_ratio": "mean",
                "inventory_tau": "mean",
            }
        )
        .fillna("Infinity")
    )
    res.to_csv(Path(snakemake.output[0]).parent / (mrio+".csv") )
