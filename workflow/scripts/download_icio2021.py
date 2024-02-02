import logging
import os
import pathlib
import sys
import traceback

import pymrio as pym

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
    f"Starting ICIO2021 (OECD) download for the year: {snakemake.wildcards.year}"
)

output_dir = pathlib.Path(snakemake.output[0]).parent
output_dir.mkdir(exist_ok=True, parents=True)

exio_meta = pym.download_oecd(
    storage_folder=output_dir,
    years=[snakemake.wildcards.year],
)
logger.info("ICIO2021 download completed successfully.")
