import re
import sys, os
import logging, traceback
from typing import Optional, Union
import pathlib
from pathlib import Path
import json
import pickle as pkl

import numpy as np
import pandas as pd
import pymrio as pym

from boario_tools.mriot import aggreg, load_mrio
from boario_tools.regex_patterns import MRIOT_AGGREG_REGEX

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

logger.info(
    f"Starting mrio aggregation to {snakemake.wildcards.aggregation} for: {snakemake.input.full_mriot_pkl}"
)

mriot_file = Path(snakemake.input.full_mriot_pkl)
mriot = load_mrio(mriot_file.stem, snakemake.config["parsed_mriot_dir"])
sectors_aggregation, regions_aggregation = re.match(MRIOT_AGGREG_REGEX, snakemake.wildcards.aggregation).groups()
save_dir = Path(snakemake.output[0]).parent

aggreg(
    mriot=mriot,
    sectors_aggregation=sectors_aggregation,
    regions_aggregation=regions_aggregation,
    save_dir=save_dir,
)
