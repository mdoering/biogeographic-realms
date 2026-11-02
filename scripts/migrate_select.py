"""Once, for the runs of M1: move each configuration's select result into its own `select/` directory, as the pipeline now writes it.

`uv run python -I scripts/migrate_select.py [--results results/explore]`

Before, screen and select wrote the same map files into `results/explore/res<r>/<run>/<config>/`, and select, which ran last, overwrote the screen's.
`explore.migrate_select` moves the select files, its stamp and its `log.json` entry into `<config>/select/`, deletes `stage-screen.json` (the screen's map files are
gone, so a screen runs again) and writes `MIGRATED.txt`; `screen.csv` stays. A configuration that was already moved, or that only has a screen, is left alone, so it
can run again. Reads and moves files under `--results` only.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

from realms.explore import migrate_select


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--results", type=Path, default=Path("results/explore"), help="the root of the explore runs, default results/explore")
    args = parser.parse_args()
    if not args.results.is_dir():
        raise SystemExit(f"FAILED {args.results} is not a directory")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S", stream=sys.stdout)
    done = migrate_select(args.results)
    print(json.dumps(done, indent=1))
    print(f"{len(done)} configurations moved, {sum(d['screen_stamp_removed'] for d in done)} screen stamps deleted")


if __name__ == "__main__":
    main()
