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

from boario_tools.mriot import aggreg

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

aggreg(
    mrio_path=snakemake.input.full_mriot_pkl[0],
    sector_aggregator_path=snakemake.input.aggreg[0],
    save_path=snakemake.output,
)
