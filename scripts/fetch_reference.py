"""Download the reference data into data/reference and record what was found, or missing, in a manifest.

Existing ChecklistBank exports only: asking for a new one is a write on production, so anything without an export is listed for a
person to decide on. Requests are sequential. The COL release comes from what GBIF production uses today.
"""

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from realms.fetch import IncompleteDownload, fetch
from realms.gbif import production_release
from realms.reference import ExportInfo, griis_datasets, latest_export

ROOT = Path("data/reference")
WCVP = 2232
WRIMS = 26853
GAZETTEERS = "https://raw.githubusercontent.com/CatalogueOfLife/col-gazetteers/main"
NATURAL_EARTH_LAND = "https://naciscdn.org/naturalearth/10m/physical/ne_10m_land.zip"
# Admin-1 places the sub-national GRIIS lists (Hawaii, the Azores) on the grid
NATURAL_EARTH_ADMIN1 = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_1_states_provinces.zip"
NATURAL_EARTH_LICENCE = "public domain (https://www.naturalearthdata.com/about/terms-of-use/: 'in the public domain ... No permission is needed to use Natural Earth')"
# WCVP's distributions are TDWG level 3 codes, whose polygons put them on the grid
# pinned to master's last commit (2021-07-15), so that a rerun cannot change the polygons
WGSRPD_LEVEL3 = "https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson"
# The repository has no licence file. Its level3 shapefile metadata (level3/level3.shp.xml) says Access_Constraints None and Use_Constraints
# 'None, but dataset is in first release', with ESRI's copyright notice for the administrative boundaries it derives from.
WGSRPD_LICENCE = (
    "none stated: github.com/tdwg/wgsrpd has no licence file, and level3.shp.xml gives Access_Constraints 'None' and Use_Constraints 'None, but dataset is "
    "in first release' with the notice 'Administrative boundary files are the intellectual property of ESRI and its licensors and are using herein with "
    "permission. (c) 1992-1997 ESRI, GMi. All Rights Reserved.' Data by Justin Moat, RBG Kew, for Brummitt 2001 (TDWG, Plant Taxonomic Database Standards No. 2)"
)
ISEA3H_GRID = "https://download.gbif.org/grids/ISEA3H/ISEA3H-{}.gpkg"
ISEA3H_RESOLUTIONS = (5, 6, 7)


def get_export(name: str, dataset_key: int, export: ExportInfo | None, manifest: dict, flavour: str = "extended") -> Path | None:
    if export is None:
        manifest[name] = {"dataset": dataset_key, "status": "missing"}
        print(f"MISSING  {name} (dataset {dataset_key}) has no suitable export", flush=True)
        return None
    path = fetch(export.url, ROOT / "exports" / flavour / f"{name}.zip", expected_size=export.size)
    manifest[name] = {"dataset": dataset_key, "export": export.key, "created": export.created, "size": export.size, "flavour": flavour, "path": str(path)}
    print(f"ok       {name}: {export.size / 1e6:.1f} MB from {export.created[:10]}", flush=True)
    return path


def main() -> None:
    manifest: dict = {}
    release = production_release()
    print(f"GBIF production uses {release}")
    # only extended exports carry environment, extinct and distributions; the classification-only one adds the taxGroup
    export = latest_export(release, classification=True, extended=True)
    groups = latest_export(release, classification=True)
    if export is None or groups is None:
        sys.exit(f"release {release} lacks an extended or a classification export: stopping")
    get_export("col", release, export, manifest)
    manifest["col"]["release"] = release
    get_export("col_groups", release, groups, manifest)

    get_export("wcvp", WCVP, latest_export(WCVP, extended=True), manifest)
    get_export("wrims", WRIMS, latest_export(WRIMS, extended=True), manifest)

    griis = griis_datasets()
    print(f"{len(griis)} GRIIS checklists", flush=True)
    for dataset in griis:
        # most only have a plain export, which is still enough: every species on a GRIIS list is alien in that country
        extended = latest_export(dataset.key, extended=True)
        chosen = extended or latest_export(dataset.key)
        get_export(f"griis/{dataset.alias}", dataset.key, chosen, manifest, "extended" if extended else "plain")

    with urllib.request.urlopen(f"{GAZETTEERS}/realm/labels.tsv", timeout=60) as response:
        realm_ids = [line.split("\t")[0] for line in response.read().decode().splitlines() if line.strip()]
    realm_ids = [i for i in realm_ids if i.lower() not in ("id", "code")]
    for gazetteer, ids in (("costello", [str(i) for i in range(1, 31)]), ("realm", realm_ids)):
        for area in ids:
            fetch(f"{GAZETTEERS}/{gazetteer}/features/{urllib.parse.quote(area)}.geojson", ROOT / "gazetteers" / gazetteer / f"{area}.geojson")
        manifest[f"gazetteer/{gazetteer}"] = {"areas": len(ids), "path": str(ROOT / "gazetteers" / gazetteer)}
        print(f"ok       gazetteer {gazetteer}: {len(ids)} areas", flush=True)

    land = fetch(NATURAL_EARTH_LAND, ROOT / "natural_earth" / "ne_10m_land.zip")
    manifest["natural_earth_land"] = {"path": str(land), "size": land.stat().st_size}

    admin1 = fetch(NATURAL_EARTH_ADMIN1, ROOT / "natural_earth" / "ne_10m_admin_1_states_provinces.zip")
    manifest["natural_earth_admin1"] = {"path": str(admin1), "size": admin1.stat().st_size, "source": NATURAL_EARTH_ADMIN1, "licence": NATURAL_EARTH_LICENCE}
    level3 = fetch(WGSRPD_LEVEL3, ROOT / "wgsrpd" / "level3.geojson")
    manifest["wgsrpd_level3"] = {"path": str(level3), "size": level3.stat().st_size, "source": WGSRPD_LEVEL3, "licence": WGSRPD_LICENCE}

    for resolution in ISEA3H_RESOLUTIONS:
        grid = fetch(ISEA3H_GRID.format(resolution), ROOT / "grids" / f"ISEA3H-{resolution}.gpkg")
        manifest[f"isea3h_{resolution}"] = {"path": str(grid), "size": grid.stat().st_size}

    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    missing = sorted(k for k, v in manifest.items() if v.get("status") == "missing")
    print(f"\n{len(manifest)} items, {len(missing)} missing: {missing}")


if __name__ == "__main__":
    try:
        main()
    except IncompleteDownload as e:
        sys.exit(f"incomplete download, rerun to resume: {e}")
