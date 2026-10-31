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

    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    missing = sorted(k for k, v in manifest.items() if v.get("status") == "missing")
    print(f"\n{len(manifest)} items, {len(missing)} missing: {missing}")


if __name__ == "__main__":
    try:
        main()
    except IncompleteDownload as e:
        sys.exit(f"incomplete download, rerun to resume: {e}")
