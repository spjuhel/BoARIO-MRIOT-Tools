import sys, os
import logging, traceback
from pathlib import Path

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


def convert(inpt, output, office_exists):
    folder = Path(output).parent
    if not office_exists:
        raise FileNotFoundError(
            "Creating csvs files require libreoffice which wasn't found. You may wan't to convert EUREGIO files by yourself if you are unable to install libreoffice"
        )
    logger.info(f"Executing: libreoffice --convert-to 'csv:Text - txt - csv (StarCalc):44,34,0,1,,,,,,,,3' {folder} {inpt}")
    os.system(f"libreoffice --convert-to 'csv:Text - txt - csv (StarCalc):44,34,0,1,,,,,,,,3' {folder} {inpt}")
    filename = Path(inpt[0]).name
    new_filename = "euregio_" + filename.split('_')[1].split('.')[0].replace('-', '_') + ".csv"
    old_path = folder / filename.replace('.ods', '-{}.csv'.format(filename.split('_')[1].split('.')[0]))
    new_path = folder / new_filename
    logger.info(f"Executing: mv {old_path} {new_path}")
    os.rename(old_path, new_path)

logger.info(f"Converting {snakemake.input.inp_file} to csv")

convert(
    snakemake.input.inp_file,
    snakemake.output.files[0],
    snakemake.params.office_exists,
)



logger.info("Conversion done !")
