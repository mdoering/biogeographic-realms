"""Download Costello et al.'s species-by-cell matrix, resuming until it is complete."""

import time
from pathlib import Path

from realms.fetch import IncompleteDownload, fetch

URL = "https://ndownloader.figshare.com/files/8622547"
DEST = Path("data/costello/obis_2009_csq.csv")
SIZE = 270216630

for attempt in range(30):
    try:
        print("complete:", fetch(URL, DEST, expected_size=SIZE))
        break
    except IncompleteDownload as e:
        print(f"attempt {attempt}: {e}", flush=True)
        time.sleep(10)
else:
    raise SystemExit("gave up after 30 attempts")
