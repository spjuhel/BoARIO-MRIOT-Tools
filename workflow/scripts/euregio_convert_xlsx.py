import sys, os
import logging, traceback

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


def convert(inpt, output, folder, office_exists, uno_exists):
    if not office_exists:
        raise FileNotFoundError(
            "Creating xlsx files require at least libreoffice which wasn't found (and optionally unoserver). You may wan't to convert EUREGIO files by yourself if you are unable to install libreoffice"
        )
    if uno_exists:
        logger.info(f"Using unoserver with 'unoconvert --port 2002 --convert-to xlsx {inpt} {output}'")
        os.system(f"unoconvert --port 2002 --convert-to xlsx {inpt} {output}")
    else:
        logger.info(f"Using libreoffice with 'libreoffice --convert-to xlsx --outdir {folder} {inpt}'")
        os.system(f"libreoffice --convert-to xlsx --outdir {folder} {inpt}")

logger.info(f"Converting {snakemake.input.inp_file} to xlsx")

convert(
    snakemake.input.inp_file,
    snakemake.output.files,
    snakemake.output.folder[0],
    snakemake.params.office_exists,
    snakemake.params.uno_exists,
)

logger.info("Conversion done !")
