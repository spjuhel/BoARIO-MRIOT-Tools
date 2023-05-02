include: "exio3_treatment.smk"
include: "euregio_treatment.smk"
include: "eora26_treatment.smk"
include: "oecd_treatment.smk"

ALL_FULL_MRIOT = (
    expand(
        "pkls/{mrio_type}/{mrio_type}_full_{year}.pkl",
        mrio_type=["euregio", "eora26", "oecd_v2021", "exiobase3"],
        year=[2000, 2010],
    )
)


rule all_mriot:
    input:
        ALL_FULL_MRIOT,
