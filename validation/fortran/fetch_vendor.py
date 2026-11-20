#!/usr/bin/env python3
"""Download the CC0 replication files used by the Fortran comparison driver.

The files come from the public Harvard Dataverse deposit
doi:10.7910/DVN/DUCHYS (CC0 1.0) and are written to
validation/fortran/vendor/, which is git-ignored. ``gss_log.txt`` is the
deposit's ``GSS/log.txt``.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENDOR = ROOT / "vendor"
DATASET = "doi:10.7910/DVN/DUCHYS"
API = "https://dataverse.harvard.edu/api/access/datafile/{file_id}"

# Dataverse file id, path under vendor/, md5, label in the deposit.
FILES = (
    (13226319, "alecia_parallel.f", "2ac53ebe4277bb83a078e4e982c902f3", "GSS/alecia_parallel.f"),
    (13226403, "eigrs.f", "700bf8ca7de2e7921fcdcd28e270a4e9", "GSS/eigrs.f"),
    (13226321, "opt.f", "4e7838e1abca1cc06770171693f38f40", "GSS/opt.f"),
    (13226353, "gss_log.txt", "f66c118cec424f54a6f770c9743c1c09", "GSS/log.txt"),
    (13226347, "data/gss2.prn", "8fca7a9e197bbc2f108f07b058eb4c16", "GSS/data/gss2.prn"),
)


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(file_id: int, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = API.format(file_id=file_id)
    temporary = dest.with_name(dest.name + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "attitudedyn-fetch/0.1"})
    with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    temporary.replace(dest)


def ensure_downloaded(vendor: Path | None = None) -> Path:
    """Download any missing or mismatched replication files and return the vendor directory."""
    destination = VENDOR if vendor is None else vendor
    for file_id, relative, expected, label in FILES:
        path = destination / relative
        if path.is_file() and _md5(path) == expected:
            continue
        print(f"Downloading {label} from {DATASET} into {path}", file=sys.stderr)
        _download(file_id, path)
        actual = _md5(path)
        if actual != expected:
            path.unlink(missing_ok=True)
            raise SystemExit(
                f"Checksum mismatch for {label}: expected {expected}, got {actual}."
            )
    return destination


def main() -> None:
    ensure_downloaded()
    print(f"Replication files are in {VENDOR}")


if __name__ == "__main__":
    main()
