from pathlib import Path
import sys
import logging, traceback

from boario_tools.mriot import euregio_convert_xlsx2csv

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

logger.info(f"Converting {snakemake.input.inp_file} to csv")
out_dir = Path(snakemake.output.files[0]).parent
euregio_convert_xlsx2csv(
    snakemake.input[0],
    out_dir,
    snakemake.params.office_exists,
)



logger.info("Conversion done !")
