#!/usr/bin/env python3
"""Assemble a wheel-build tree that contains the Python project and the C++ core.

``pip install ./python`` finds the core at ``../core``. cibuildwheel packs one
directory into an sdist, so this script copies the Python project and places
that core at ``./core`` inside ``wheel-stage/``. Wheel tests still read the
GSS files from the checkout.
"""

from __future__ import annotations

import shutil
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PYTHON_SOURCE = REPOSITORY_ROOT / "python"
CORE_SOURCE = REPOSITORY_ROOT / "core"
STAGE_ROOT = REPOSITORY_ROOT / "wheel-stage"
SKIPPED_DIRECTORY_NAMES = {
    "__pycache__",
    ".pytest_cache",
    "build",
    "_skbuild",
    "dist",
    "wheel-stage",
}


def skip_python_entries(directory: str, names: list[str]) -> set[str]:
    """Leave build products and caches out of the staged tree."""
    skipped: set[str] = set()
    for name in names:
        if name in SKIPPED_DIRECTORY_NAMES or name.endswith(".egg-info"):
            skipped.add(name)
        elif name.endswith((".pyc", ".so", ".pyd")):
            skipped.add(name)
    return skipped


def main() -> None:
    if STAGE_ROOT.exists():
        shutil.rmtree(STAGE_ROOT)
    shutil.copytree(PYTHON_SOURCE, STAGE_ROOT, ignore=skip_python_entries)
    staged_core = STAGE_ROOT / "core"
    if staged_core.exists():
        shutil.rmtree(staged_core)
    shutil.copytree(
        CORE_SOURCE,
        staged_core,
        ignore=shutil.ignore_patterns("tests", "test_core", "*.o"),
    )


if __name__ == "__main__":
    main()
