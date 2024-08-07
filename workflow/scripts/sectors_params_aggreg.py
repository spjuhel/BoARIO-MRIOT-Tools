from importlib import resources
import re
from boario_tools.mriot import find_sectors_agg
from boario_tools.regex_patterns import MRIOT_BASENAME_REGEX, MRIOT_AGGREG_SECTORS_REGEX
import pandas as pd

import sys, traceback
import logging
import pandas as pd

sys.path.append("../others")

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

base_aggreg = re.search(MRIOT_BASENAME_REGEX+"_"+MRIOT_AGGREG_SECTORS_REGEX,snakemake.input.sectors_config).groupdict()["mrio_aggreg_sectors"]
sectors_df = pd.read_csv(snakemake.input.sectors_config, index_col=0, decimal=".")
if snakemake.params["alt_aggregation_file"] is not None:
    aggregation_master_df = pd.read_excel(
        snakemake.params.alt_aggregation_file, sheet_name=0, index_col=0
    )
else:
    with resources.path("boario_tools.data", "aggregation_files") as agg_path:
        aggregation_master_df = find_sectors_agg(
            snakemake.wildcards.mriot_name,
            base_aggreg,
            snakemake.wildcards.sectors_aggregation_nofull,
            agg_files_path=agg_path,
        )

res = (
    aggregation_master_df.join(sectors_df)[
        [
            "new sector",
            "affected",
            "rebuilding_factor",
            "inventory_size",
            "productive_capital_to_va_ratio",
            "inventory_tau",
        ]
    ]
    .groupby("new sector")
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

res.to_csv(snakemake.output[0])
