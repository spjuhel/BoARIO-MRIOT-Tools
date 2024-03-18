import sys, os
import logging, traceback
import pymrio as pym
import pathlib
import pickle as pkl

from boario_tools.mriot import lexico_reindex

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


def parse_oecd_v2018(mrio_zip: str, output: str, year:int):
    logger.info(
        "Make sure you use the same python environment as the one loading the pickle file (especial pymrio and pandas version !)"
    )
    try:
        logger.info("Your current environment is: {}".format(os.environ["CONDA_PREFIX"]))
    except KeyError:
        logger.info("Could not find CONDA_PREFIX, this is normal if you are not using conda.")

    mrio_path = pathlib.Path(mrio_zip)
    mrio_pym = pym.parse_oecd(path=mrio_path, year=year)
    logger.info("Removing unnecessary IOSystem attributes")
    attr = [
        "Z",
        "Y",
        "x",
        "A",
        "L",
        "unit",
        "population",
        "meta",
        "__non_agg_attributes__",
        "__coefficients__",
        "__basic__",
    ]
    tmp = list(mrio_pym.__dict__.keys())
    for at in tmp:
        if at not in attr:
            delattr(mrio_pym, at)
    assert isinstance(mrio_pym, pym.IOSystem)
    logger.info("Done")
    logger.info("Computing the missing IO components")
    mrio_pym.calc_all()
    logger.info("Done")
    logger.info("Re-indexing lexicographicaly")
    mrio_pym = lexico_reindex(mrio_pym)
    logger.info("Done")
    save_path = pathlib.Path(output)
    logger.info("Saving to {}".format(save_path.absolute()))
    save_path.parent.mkdir(parents=True, exist_ok=True)
    setattr(mrio_pym, "monetary_factor", 1000000)
    mrio_pym.meta.change_meta("name", "oecd-icio2021")
    with open(save_path, "wb") as f:
        pkl.dump(mrio_pym, f)


parse_oecd_v2018(snakemake.input[0], snakemake.output[0], snakemake.wildcards.year)
