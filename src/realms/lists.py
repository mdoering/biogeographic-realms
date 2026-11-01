"""Taxa whose records say nothing about where they occur naturally, as a short list that can be defended name by name.

Every entry was checked against COL release 315557 (see DECISIONS.md): the name is an accepted species there, and nothing wild is filed
under it. A species whose wild forms COL files under the same accepted name stays out, as dropping it would drop the wild records too.
"""

import logging

import polars as pl

log = logging.getLogger(__name__)

# name -> (kingdom, reason)
NO_WILD_RANGE: dict[str, tuple[str, str]] = {
    "Homo sapiens": ("Animalia", "humans: records are people, observed or in photographs, wherever they were"),
    "Bos taurus": ("Animalia", "domestic cattle; the wild aurochs is Bos primigenius, a species of its own in COL"),
    "Bos indicus": ("Animalia", "zebu, domestic cattle with no wild population"),
    "Bubalus bubalis": ("Animalia", "domestic water buffalo; the wild water buffalo is Bubalus arnee in COL"),
    "Capra hircus": ("Animalia", "domestic goat; the wild bezoar is Capra aegagrus in COL (one local synonym, the Chiltan markhor, is filed here)"),
    "Equus caballus": ("Animalia", "domestic horse, whose feral herds have no natural range; the wild horse is Equus ferus in COL"),
    "Equus asinus": ("Animalia", "domestic donkey; the wild ass is Equus africanus in COL"),
    "Canis familiaris": ("Animalia", "domestic dog; COL holds it apart from the wolf, Canis lupus. The dingo is filed here, a feral dog with no natural range"),
    "Felis catus": ("Animalia", "domestic cat; the wildcat is Felis silvestris in COL"),
    "Sus domesticus": ("Animalia", "domestic pig; the wild boar is Sus scrofa in COL"),
    "Camelus dromedarius": ("Animalia", "dromedary, domestic, and extinct in the wild"),
    "Camelus bactrianus": ("Animalia", "domestic Bactrian camel; the wild camel is Camelus ferus in COL"),
    "Lama glama": ("Animalia", "llama, domestic; the guanaco is Lama guanicoe in COL"),
    "Vicugna pacos": ("Animalia", "alpaca, domestic; the vicuna is Vicugna vicugna in COL"),
    "Cavia porcellus": ("Animalia", "domestic guinea pig; wild cavies are Cavia aperea in COL"),
    "Triticum aestivum": ("Plantae", "bread wheat, a cultigen of the hexaploid with no wild population"),
    "Brassica napus": ("Plantae", "rapeseed, a cultigen of two other species with no wild population"),
}

# Entries looked at and left out, with what COL files under the accepted species. `Canis lupus` and `Sus scrofa` are left out for the same
# reason, as the wolf and the wild boar are the wild species, but nobody proposed them as entries.
NO_WILD_RANGE_DROPPED: dict[str, str] = {
    "Ovis aries": "COL files the wild mouflon under it: Ovis orientalis is a synonym, as are Ovis orientalis ophion (the Cyprus mouflon) and Ovis aries musimon",
    "Gallus gallus domesticus": "COL has no such species: Gallus domesticus is a synonym of Gallus gallus, which holds the wild red junglefowl (subspecies bankiva, spadiceus, jabouillei, murghi)",
    "Zea mays": "COL files the wild teosintes under it as accepted subspecies parviglumis and huehuetenangensis",
    "Hordeum vulgare": "COL files wild barley under it, as accepted subspecies spontaneum and the synonym Hordeum spontaneum",
    "Glycine max": "COL files wild soybean under it, as accepted subspecies soja and the synonym Glycine soja",
    "Solanum tuberosum": "COL files wild potatoes under it as species-rank synonyms (Solanum boyacense, Solanum parmentieri, Solanum maglia among 186), and WCVP gives it a native range in the Andes",
    "Oryza sativa": "COL files the wild perennial rice under it, as the synonym Oryza perennis",
    "Sorghum bicolor": "COL files the wild sorghum under it, as accepted subspecies verticilliflorum",
    "Capsicum annuum": "COL files the wild chiltepin under it, as accepted variety glabriusculum",
}

ACCEPTED_NAMES = ("accepted", "provisionally accepted")


def resolve_no_wild_range(index: pl.DataFrame, entries: dict[str, tuple[str, str]] | None = None) -> tuple[frozenset[str], dict[str, str]]:
    """The COL ids of the species in `entries`, and what became of each entry that gave none.

    An entry counts only where its name is that of an accepted or provisionally accepted species of its kingdom: a name that COL
    reads as a synonym of another species (`Gallus domesticus`, of the junglefowl) would take that species' wild records with it, so
    it is dropped. COL holds some names under several ids (a duplicate of the same species), and all of them count. `index` is
    `natives.name_index`.
    """
    entries = NO_WILD_RANGE if entries is None else entries
    ids: set[str] = set()
    left_out: dict[str, str] = {}
    own = index.filter(pl.col("status").is_in(ACCEPTED_NAMES))
    for name, (kingdom, _) in entries.items():
        found = own.filter(pl.col("name") == name, pl.col("kingdom") == kingdom)["species_id"].to_list()
        if found:
            ids.update(found)
            continue
        as_synonym = index.filter(pl.col("name") == name, pl.col("kingdom") == kingdom)["species_id"].to_list()
        left_out[name] = f"a synonym of {sorted(as_synonym)}" if as_synonym else f"no accepted {kingdom} species of that name"
        log.warning("NO_WILD_RANGE: %s dropped, %s", name, left_out[name])
    return frozenset(ids), left_out
